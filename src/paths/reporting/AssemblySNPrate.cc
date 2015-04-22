///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* AssemblySNPrate.  Find SNPs in the assembly 
   contigs. AlignReads must have been run first.  
   To work correctly, the bandwidth entries in the read locations have to be set.  
   If we wanted to work around this we could pass a default bandwidth to the 
   alignment code. */

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS


#include "MainTools.h"
#include "Superb.h"
#include "paths/ReadLoc.h"
#include "system/LockedData.h"
#include "system/ParsedArgs.h"
#include "system/SpinLockedData.h"
#include "VecUtilities.h"
#include <omp.h>


int main(int argc, char *argv[]){
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_OrDefault(SUBDIR, "test");
  CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
  CommandArgument_String_OrDefault_Doc(TIGS, "", 
	 "if unspecified, process all contigs;"
	 " otherwise it is one of the following: \n"
	 "(a) a list of contig ids (in ParseIntSet format) or \n"
	 "(b) the letter 's' followed by a list of scaffolds or \n"
	 "(c) s<scaffold id>.<list of indices of contigs in the scaffold"
         "(d) the letter 'p' followed by a percentage of a sum of contig lengths");
  CommandArgument_String_OrDefault_Doc(LIB_TYPES, "frag", 
	 "can be {frag,jump,long_jump} or a subset");
  CommandArgument_Int_OrDefault_Doc(MIN_CVRG, 30,
	 "min base coverage to compute SNPs");
  CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
	 "if specified, file extension is .READLOCS_PREFIX.readlocs "
				       "instead of .readlocs");
  CommandArgument_Bool_OrDefault_Doc( WRITE_ALL, False,
	 "if true, save output of all positions in <ASSEMBLY>.SNPs.out" );
  CommandArgument_Bool_OrDefault_Doc( ARCHIVE, False,
	 "if true, save output as <ASSEMBLY>.SNPrates.out" );
  CommandArgument_Int_OrDefault_Doc( N_THREADS, 16,
	 "number of threads to be used" );
  CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, 133333,
    "Seed value for the random generator." );
  EndCommandArguments;
  
  N_THREADS = configNumThreads(N_THREADS);
  omp_set_num_threads(N_THREADS);

  // Define directories.
  
  String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 
  
  // Output stream.
     
  ofstream archive_out;
  if ( ARCHIVE ) {
    String archive_file =  sub_dir + "/" + ASSEMBLY + ".SNPrates.out";
    cout << "Sending output to " << archive_file << "\n" << endl;
    archive_out.open( archive_file.c_str( ) );
    PrintCommandPretty( archive_out );
  }
  ostream &out = ARCHIVE ? archive_out : * (ostream *) &cout;

  String snps_file =  sub_dir + "/" + ASSEMBLY + ".SNPs.out";
  ofstream snps_out( snps_file.c_str() );
  // find out which classes of reads to use

  vec<Bool> lib_types( 3, False );
  vec<String> lib_types_s;
  ParseStringSet( LIB_TYPES, lib_types_s );
  for ( int i = 0; i < lib_types_s.isize( ); i++ ){    
    if ( lib_types_s[i] == "frag" )
      lib_types[0] = True;
    else if ( lib_types_s[i] == "jump" )
      lib_types[1] = True;
    else if ( lib_types_s[i] == "long_jump" )
      lib_types[2] = True;
    else{    
      out << "Illegal LIB_TYPES." << endl;
      return 1;    
    }    
  }

  // Parse TIGS.
  
  vec<int> tigIds;
  {    
    if (TIGS == ""){    
      int n_tigs = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY
					     + ".contigs.fastb" );
      for ( int j = 0; j < n_tigs; j++ )
	tigIds.push_back(j);    
    }
    else if (TIGS.Contains("s", 0)) {    
      TIGS = TIGS.After("s");
      vec<superb> scaffolds;
      ReadSuperbs(sub_dir + "/" + ASSEMBLY + ".superb", scaffolds);
      if (TIGS.Contains(".")) {    
	int scaffold = TIGS.Before(".").Int();
	ForceAssertLt(scaffold, scaffolds.isize());
	vec<int> spos;
	ParseIntSet(TIGS.After("."), spos);
	for (int j = 0; j < spos.isize(); j++)
	  tigIds.push_back(scaffolds[scaffold].Tig(spos[j]));    
      }
      else{    
	vec<int> s;
	ParseIntSet(TIGS, s);
	for (int i = 0; i < s.isize(); i++){ 
	  int scaffold = s[i];
	  ForceAssertLt(scaffold, scaffolds.isize());
	  for (int j = 0; j < scaffolds[scaffold].Ntigs(); j++)
	    tigIds.push_back(scaffolds[scaffold].Tig(j));
	}    
      }   
    }
    else if (TIGS.Contains("p", 0)) {    
      double perc = TIGS.After("p").Double();
      out << "percent of tig len sum required = " << perc << endl;
      vecbasevector tseqs( sub_dir + "/" + ASSEMBLY + ".contigs.fastb" );
      size_t n_tigs = tseqs.size();
      longlong totLen = 0;
      for ( size_t i = 0; i < n_tigs; i++ )
	totLen += tseqs[i].size();
      
      vec<int> stigids( n_tigs, vec<int>::IDENTITY );
      srand(RANDOM_SEED);
      random_shuffle( stigids.begin(), stigids.end() );
      longlong reqLen = (longlong) (totLen * perc/100.0); 
      longlong currLen = 0;
      
      for (size_t k = 0; k < stigids.size(); k++){
	if ( currLen < reqLen ){
	  tigIds.push_back( stigids[k] );
	  currLen += tseqs[ stigids[k] ].size();
	}
	else break;
      }
    }
    else ParseIntSet(TIGS, tigIds);    
  }

  UniqueSort( tigIds );
  PRINT_TO( out, tigIds.size() );
  //out << "tigIds: "; tigIds.Print(out); out << endl;
  
  
  // Process contigs.
  snps_out << ">TYPE\t" << "TIG" << "\t" << "POS" << "\t" << "CT0\tCT1\tCT2\tCT3\tCT4\tCT5\t" << "BASE\tSNP1\tSNP2" << "\n"; 

  String head = sub_dir + "/" + ASSEMBLY;
  if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
  ForceAssert( IsRegularFile( head + ".readlocs") );
  read_locs_on_disk locs_file( head, run_dir );
  vec<int> tEffLengths( tigIds.size(), 0 );
  vec<int> tnSNPs( tigIds.size(), 0 );
  vec<String> tOutSNPs( tigIds.size() );
  out << Date() << " reading sequences and readlocs for selected tigs" << endl;
  vec< vec<read_loc> > v_locs( tigIds.size() );
  vecbasevector tigs;
  for ( int it = 0; it < tigIds.isize( ); it++ ){
    if ( tigIds.size( ) < 100 || it % ( tigIds.size()/100 ) == 0 ){ 
      out << "." << flush;
    }
    int TIG = tigIds[it];
    vec<read_loc> locs;
    locs_file.LoadContig( TIG, locs );
    v_locs[it] = locs;
  }
  out << endl;
  
  tigs.Read( sub_dir + "/" + ASSEMBLY + ".contigs.fastb", tigIds ); 
  PRINT2_TO( out, tigIds.size(), tigs.size() );

  out << Date() << " computing pileups" << endl;
  SpinLockedData lock;
  #pragma omp parallel for
  for ( int it = 0; it < tigIds.isize( ); it++ ){
    
    if ( tigIds.size( ) < 100 || it % ( tigIds.size()/100 ) == 0 ){
      SpinLocker locker(lock);
      out << "." << flush;
    }
    int TIG = tigIds[it];
    
    vec<dumbcall> calls;
    Pileup( tigs[it], v_locs[it], run_dir, calls, 
	    lib_types[0], lib_types[1], lib_types[2] );

    for ( int j = 0; j < tigs[it].isize( ); j++ ){
      vec<int> order( 6, vec<int>::IDENTITY );
      vec<int> ccalls(6, 0);
      for ( int k = 0; k < 6; k++ )
	ccalls[k] = calls[j].base[k];
      if ( Sum(ccalls) > MIN_CVRG ){
	ReverseSortSync( ccalls, order );
	double sum2   = ccalls[0] + ccalls[1];
	double sumAll = Sum( ccalls );
	if ( sum2/sumAll > 0.95 && order[0] < 4 && order[1] < 4 ){
	  tEffLengths[it]++;
	  if ( ccalls[1] / sum2 > 0.4 ){
	    tnSNPs[it]++;
	    #pragma omp critical
	    snps_out << ">SNP*\t" << TIG << "\t" << j << "\t" << calls[j].base[0] << "\t" << calls[j].base[1] << "\t" << calls[j].base[2] << "\t" << calls[j].base[3] << "\t" << calls[j].base[4] << "\t" << calls[j].base[5] << "\t" << as_base(tigs[it][j]) << "\t" << as_base(order[0]) << "\t" << as_base(order[1]) << endl; 
	  }else if ( WRITE_ALL )
	    #pragma omp critical
	    snps_out << ">HAP \t" << TIG << "\t" << j << "\t" << calls[j].base[0] << "\t" << calls[j].base[1] << "\t" << calls[j].base[2] << "\t" << calls[j].base[3] << "\t" << calls[j].base[4] << "\t" << calls[j].base[5] << "\t" << as_base(tigs[it][j]) << endl;
	}
      }
    }
  }
  out << endl;

  double nSNPs = Sum( tnSNPs );
  double effLength = Sum( tEffLengths );

  double SNPrate = (double) nSNPs / (double) effLength;
  double sigma = sqrt( nSNPs ) / (double) effLength;
  PRINT4_TO( out, SNPrate, 3*sigma, nSNPs, effLength );

  out << Date() << " AssemblySNPrate done!" << endl;
}    

