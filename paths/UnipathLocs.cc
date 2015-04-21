/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------------

// UnipathLocs.  Generate read locations on normal unipaths, as a 
// vec<read_location_short>, where "normal" is defined by MAX_COPY_NUMBER and
// MIN_KMERS.  As compared to using all the unipaths, this saves space.
// Sort by contig, and for fixed contig, by read.  Provide a separate index by 
// read.  Files created:
// reads.unilocs.K[.rindex].

#include "Basevector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "ReadLocation.h"
#include "Vec.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_Int_OrDefault(MAX_COPY_NUMBER, 10);
     CommandArgument_Int_OrDefault(MIN_KMERS, 1);
     CommandArgument_String_OrDefault(READS, "reads");
     CommandArgument_String_OrDefault(PATHS, "paths");
     CommandArgument_String_OrDefault(UNIPATHS, "unipaths");
     CommandArgument_Bool_OrDefault_Doc( DUMP_COVERAGES, False, "Writes a report to RUN/UnipathLocs.coverages.log" );
     EndCommandArguments;

     // Set up directories.

     String run_dir = PRE + "/" + DATA + "/" + RUN;

     // Read in read paths, unipaths, reads, and pairing data.

     vecKmerPath paths( run_dir + "/" + READS + "." + PATHS + ".k" + ToString(K) );
     vecKmerPath paths_rc( run_dir + "/" + READS + "." + PATHS + "_rc.k" + ToString(K) );
     vecKmerPath unipaths( run_dir + "/" + READS + "." + UNIPATHS + ".k" + ToString(K) );
     size_t nreads = paths.size( );
     ForceAssertEq( nreads, paths_rc.size( ) );
     vecKmerPath::size_type nuni = unipaths.size( );
     BREAD2( run_dir + "/" + READS + "." + UNIPATHS + "db.k" + ToString(K), 
          vec<tagged_rpint>, unipathsdb );

     cout << "Reads to place  : " << nreads << endl;
     cout << "Unipath count   : " << nuni << endl;

     // Compute length of each unipath.

     vec<int> ulen(nuni);
     for ( vecKmerPath::size_type i = 0; i < nuni; i++ )
          ulen[i] = unipaths[i].KmerCount( );

     // Read in output of UnipathCoverage.

     VecPdfEntryVec cp;
     cp.ReadAll( (run_dir + "/" + READS + "." + UNIPATHS + ".predicted_count.k" + ToString(K)).c_str() );
     ForceAssertEq( cp.size( ), nuni );


     // Define the normal unipaths.
     vec<int> predicted_copyno(nuni, -1);
     vec<Bool> normal(nuni, False);
     for ( vecKmerPath::size_type i = 0; i < nuni; i++ ) {
       GetMostLikelyValue( predicted_copyno[i], cp[i] );
       if ( predicted_copyno[i] > MAX_COPY_NUMBER ) continue;
       if ( ulen[i] < MIN_KMERS ) continue;
       normal[i] = True;    
     }

     // Find read placements on unipaths.  

     vec<read_location_short> ulocs;
     for ( size_t id = 0; id < nreads; id++ )
     {    for ( int pass = 1; pass <= 2; pass++ )
          {    const KmerPath& p = ( pass == 1 ? paths[id] : paths_rc[id] );
               static vec< pair<int,int> > uo;
               uo.clear( );
               for ( int j = 0; j < p.NSegments( ); j++ )
               {    const KmerPathInterval& I = p.Segment(j);
                    static vec<longlong> locs;
                    Contains( unipathsdb, I, locs );
                    for ( int u = 0; u < locs.isize( ); u++ )
                    {    const tagged_rpint& t = unipathsdb[ locs[u] ];
                         int uid = t.PathId( );
                         if ( !normal[uid] ) continue;
                         longlong offset = t.Start( ) - I.Start( );
                         for ( int r = 0; r < j; r++ )
                              offset += p.Segment(r).Length( );
                         for ( int r = 0; r < t.PathPos( ); r++ )
                              offset -= unipaths[uid].Segment(r).Length( );
                         uo.push_back( make_pair( uid, int(offset) ) );    }    }
               UniqueSort(uo);
               for ( int t = 0; t < uo.isize( ); t++ )
               {    int uid = uo[t].first, offset = uo[t].second;
                    ulocs.push_back( read_location_short( id, uid, -offset, 
                         ( pass == 1 ? ForwardOr : ReverseOr ) ) );    }    }    }
     sort( ulocs.begin( ), ulocs.end( ), cmp_contig_read );

     // Generate output files.

     VecIntVec index;
     index.resize(nreads);
     for ( vec<read_location_short>::size_type ii = 0; ii < ulocs.size(); ++ii )
         index[ ulocs[ii].ReadId() ].push_back(ii);
     for ( VecIntVec::iterator itr(index.begin()), end(index.end());
             itr != end; ++itr )
         sort( itr->begin(), itr->end() );

     String filehead = run_dir + "/" + READS + ".unilocs." + ToString(K);
     BinaryWriter::writeFile( filehead, ulocs );
     index.WriteAll( filehead + ".indexr" );

     // Compute some basic stats
     int unique = 0, multiple = 0, unplaced = 0;
     for (size_t i = 0; i < nreads; i++ ) {
       IntVec::size_type placements = index[i].size();
       if (placements == 0)
	 unplaced++;
       else if (placements <= 2)
	 unique++;
       else
	 multiple++;
     }

     cout << "Normal unipaths : " << normal.CountValue(True)  << endl;
     cout << "Uniquely placed reads : " << unique << endl;
     cout << "Multiply placed reads : " << multiple << endl;
     cout << "Unplaced reads        : " << unplaced << endl;



  // Report on the coverage of each unipath by read locs.
  if ( DUMP_COVERAGES ) {
    
    // Get the genome size.
    String genome_size_file = PRE + "/" + DATA + "/genome.size";
    if ( !IsRegularFile( genome_size_file ) ) {
      cout << "Cannot calculate coverages: File " << genome_size_file
	   << " does not exist.  Ignoring DUMP_COVERAGES." << endl;
      return 0;
    }
    longlong genome_size = FirstLineOfFile( genome_size_file ).Int( );
    double expected_coverage = double( nreads ) / genome_size;
    
    // Open the output file.
    String coverages_file = run_dir + "/UnipathLocs.coverages.out";
    cout << "DUMP_COVERAGES: Writing to " << coverages_file << endl;
    ofstream out( coverages_file.c_str( ) );
    out << "Genome size = " << genome_size << endl;
    
    
    // Find the number of reads with a read_location on each unipath.
    vec<int> locs_per_unipath( nuni, 0 );
    for ( int i = 0; i < ulocs.isize( ); i++ )
      locs_per_unipath[ ulocs[i].Contig( ) ]++;
    

    for ( vecKmerPath::size_type i = 0; i < nuni; i++ ) {
      int CN = predicted_copyno[i];
      
      // Find coverage of unipath by reads.
      double coverage = double( locs_per_unipath[i] ) / ulen[i];
      
      out << "Unipath " << i << ":\tLength " << ulen[i] << "\tCopy # "
	  << CN << "\tN reads " << locs_per_unipath[i]
	  << "\t-> Coverage by read_locations: " << coverage
	  << "\t(expected: " << CN << ")" << endl;
    }
    
    out.close( );
  }
  
  cout << Date( ) << ": Done!" << endl;
  return 0;
}
