///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Read SAM records from standard input, generate files
// OUT_HEAD.{fastb,qualb,qltout,pairto,pairtob,pairto_index,names}.
// Currently only set up for one library.

#include "Basevector.h"
#include "FastaFileset.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "VecString.h"
#include "lookup/LibInfo.h"
#include "lookup/LookAlign.h"
#include "lookup/SAM2CRD.h"

int main( int argc, char** argv )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(SAM,"/dev/stdin",
       "input file, defaults to stdin");
     CommandArgument_String(OUT_HEAD);
     CommandArgument_String_OrDefault_Doc(LIBINFO,
       "", "text file with <lib name> <insert size> <std dev>");
     CommandArgument_Int_OrDefault_Doc(SEP, 0,
       "mean separation for reads from unknown or unspecified library");
     CommandArgument_Int_OrDefault_Doc(DEV, 0,
       "std dev of separation for reads from unknown or unspecified library");
     CommandArgument_String_OrDefault_Doc(LIB_NAME, "",
       "library name for reads from unknown library");
     CommandArgument_Bool_OrDefault_Doc(REQUIRE_LIBINFO, False,
       "library names must be present in LIBINFO file - unknown librariess not allowed");
     CommandArgument_Int_OrDefault_Doc(TIG, -1, "force target_id to be this");
     CommandArgument_Bool_OrDefault_Doc( FLIP, False,
       "reverse-complement reads and reverse quals.");
     CommandArgument_Bool_OrDefault( PRESERVE_ORIENTATION, False );
     CommandArgument_Bool_OrDefault(KEEP_DUPS, True);
     CommandArgument_Bool_OrDefault(MAPPED_PAIRS_ONLY, False);
     CommandArgument_Bool_OrDefault(PFONLY, False);
     CommandArgument_Bool_OrDefault(LOG_TO_CERR, True);
     CommandArgument_Int_OrDefault_Doc(MIN_MAPQ, 0, "minimum mapping quality");
     CommandArgument_Bool_OrDefault_Doc(WRITE_ALIGNS, True,
       "Output <OUT_HEAD>.{qltout,mpaq} (alignment information)");
     CommandArgument_Bool_OrDefault_Doc(WRITE_PAIRS, True,
       "Output <OUT_HEAD>.pairs (new-format pairs file)");
     CommandArgument_Bool_OrDefault_Doc(WRITE_NAMES, True,
       "Output <OUT_HEAD>.{names} (original read names)");
     CommandArgument_Bool_OrDefault_Doc(NAMES_PLUS, False,
       "read names include .1 or .2 depending on first read in pair");
     CommandArgument_Bool_OrDefault_Doc(APPLY_CLIPPING, False,
        "Experimental feature to clip output reads" );
     CommandArgument_Bool_OrDefault_Doc(USE_OQ, False,
       "Use the original quality scores if the OQ tag exists "
       "(otherwise fall back to the regular quality scores).");
     CommandArgument_UnsignedInt_OrDefault_Doc(NOMINAL_READ_LEN, 0,
        "Subtract twice this value from insert size to get an approximate value for pair separation.");
     CommandArgument_Bool_OrDefault_Doc(PRINT_VISUAL, False,
                                        "print visualizations of the "
                                        "alignments in the QLTOUT file.");
     CommandArgument_String_OrDefault_Doc(REF_FASTA, "", "the alignment "
                                          "reference, required for "
                                          "PRINT_VISUAL=True, when present "
                                          "allows mismatches to be correctly "
                                          "noted in parseable QLTOUT lines");
     EndCommandArguments;

     vecbasevector seqs;
     vecbasevector refs;
     vec<size_t> ref_indices;
     vecqualvector quals;
     vec<look_align_x> alns;
     vec<pairinfo> pairinfoVec;
     vecString namesv;
     vecString libNames;
     Logger const& logger = LOG_TO_CERR ? Logger(std::cerr) : Logger::nullLogger();
     SAM::SAMFile* samfile;
     if (SAM == "/dev/stdin" || !SAM.EndsWith(".bam"))
     {
         samfile = new SAM::SAMFile(SAM,logger);
     }
     else
     {
        samfile = new SAM::BAMFile(SAM,"",false,logger);
     }

     if (PRINT_VISUAL && REF_FASTA.empty())
     {
         std::cerr << "PRINT_VISUAL=True requires REF_FASTA to be set "
                   << "to the relevant reference." << std::endl;
         CRD::exit(EXIT_FAILURE);
     }

     if (!REF_FASTA.empty())
     {
         // load the reference chromosomes from REF_FASTA and build a map
         // relating chromosome names to indices
         vecString ref_names;
         FastFetchReads(refs, &ref_names, REF_FASTA);
         validateReferenceDictionary(*samfile,
                                 vec<String>(ref_names.begin(),ref_names.end()),
                                 refs);
     }
     
     vec<Bool> first_in_pair;
     SAM2CRD( *samfile, seqs, quals, alns, pairinfoVec, namesv, first_in_pair,
              libNames, MAPPED_PAIRS_ONLY, KEEP_DUPS, USE_OQ, PFONLY, 
              APPLY_CLIPPING, NAMES_PLUS );
     if (NAMES_PLUS)
     {    for ( int64_t i = 0; i < (int64_t) namesv.size( ); i++ )
               namesv[i] += ( first_in_pair[i] ? ".1" : ".2" );    }

     delete samfile;

     if ( TIG >= 0 )
     {    for ( int i = 0; i < alns.isize( ); i++ )
               alns[i].target_id = TIG;    }
     

     if ( FLIP ){
       cout << Date( ) << ": rc-ing query bases" << endl;
       for (ulonglong i=0; i < seqs.size(); i++)
	 seqs[i].ReverseComplement();
       for (ulonglong i=0; i < quals.size(); i++ )
	 reverse( quals[i].begin( ), quals[i].end( ) );
     }

     if (WRITE_ALIGNS) {
       Ofstream( aout, OUT_HEAD + ".qltout" );
       Ofstream( mapout, OUT_HEAD + ".mapq" );
       for ( int j = 0; j < alns.isize( ); j++ )
	 {    
	   if ( (int) alns[j].mapQ( ) < MIN_MAPQ ) continue;
	   if (REF_FASTA.empty())
	     {
	       alns[j].PrintParseable(aout);
	     }
	   else
	     {
	       alns[j].PrintParseable(aout, &seqs[alns[j].QueryId()],
				      &refs[ref_indices[alns[j].TargetId()]]);
	     }
	   if (PRINT_VISUAL)
	     {
	       alns[j].PrintVisual(aout, seqs[alns[j].QueryId()],
				   refs[ref_indices[alns[j].TargetId()]], false);
	     }         
	   mapout << (int) alns[j].mapQ( ) << "\n";    }
     }

     if (WRITE_NAMES)
       namesv.WriteAll( OUT_HEAD + ".names" );

     PairsManager pairs( seqs.size() );
     if ( WRITE_PAIRS ) {

       if ( LIBINFO == "" )
       {     for ( vecString::iterator itr(libNames.begin()),
	 	   end(libNames.end()); itr != end; ++itr ) 
             {     String sam_lib_name = *itr;
	           if (sam_lib_name == "unknown")
	                sam_lib_name = (LIB_NAME != "" ? LIB_NAME : "unknown");
	           pairs.addLibrary(SEP, DEV, sam_lib_name);    }    }

       else if (IsRegularFile( LIBINFO )) {
	 LibInfoDB libInfoDB(LIBINFO);
	 for ( vecString::iterator itr(libNames.begin()),
		 end(libNames.end()); itr != end; ++itr ) {
	   LibInfo const* pInfo = libInfoDB.getInfo(*itr);

	   if ( pInfo ) {
	     pairs.addLibrary(pInfo->mMean - 2 * NOMINAL_READ_LEN,
			      pInfo->mStdDev, *itr);
	   } else {
	     std::cout << (REQUIRE_LIBINFO ? "ERROR: " : "Warning: ")
		       << "Cannot find entry for library: " << (*itr) 
		       << " in library information file: " << std::endl
		       << LIBINFO << std::endl;
	     
	     if (REQUIRE_LIBINFO) // We must have library information
		 exit(1);
	     
	     std::cout << "Using default SEP and DEV instead" << std::endl;
	     
	     String sam_lib_name = *itr;
	     if (sam_lib_name == "unknown")
	       sam_lib_name = (LIB_NAME != "" ? LIB_NAME : "unknown");
	     pairs.addLibrary(SEP, DEV, sam_lib_name);
	   }
	 }
	 
       } else if ( REQUIRE_LIBINFO ) {
	 std::cout << "ERROR: Cannot find library information file: " << std::endl
		   << LIBINFO << std::endl;
	 exit(1);

       } else {
	 std::cout << "Warning: No library information file, using default SEP and DEV "
		   << "for all libraries instead" << std::endl;
	 for ( vecString::iterator itr(libNames.begin()),
		 end(libNames.end()); itr != end; ++itr ) {
	   String sam_lib_name = (LIB_NAME == "" ? *itr : LIB_NAME);
	   pairs.addLibrary(SEP, DEV, sam_lib_name);
	 }
       }
       
       // Create a PairsManager with pairing information.
       typedef vec<pairinfo>::iterator Itr;
       for ( Itr itr(pairinfoVec.begin()), end(pairinfoVec.end());
               itr != end;
               ++itr )
           pairs.addPairToLib(itr->readID1,itr->readID2,itr->libraryID,false);

       // Write to the pairs file and/or the pairto file.
       if ( WRITE_PAIRS )
	 pairs.Write( OUT_HEAD + ".pairs" );
      
    }

     if ( PRESERVE_ORIENTATION )
     {    vec<Bool> flipped( seqs.size( ), False );
          for ( int i = 0; i < alns.isize( ); i++ )
          {    int64_t id1 = alns[i].query_id;
               int64_t id2 = pairs.getPartnerID(id1);
               if ( alns[i].Fw1( ) && id2 >= 0 && !flipped[id2] )
               {    seqs[id2].ReverseComplement( );
                    quals[id2].ReverseMe( );
                    flipped[id2] = True;    }
               if ( alns[i].Rc1( ) && !flipped[id1] )
               {    seqs[id1].ReverseComplement( );
                    quals[id1].ReverseMe( );
                    flipped[id1] = True;    }    }    }

     seqs.WriteAll( OUT_HEAD + ".fastb" );
     quals.WriteAll( OUT_HEAD + ".qualb" );
}
