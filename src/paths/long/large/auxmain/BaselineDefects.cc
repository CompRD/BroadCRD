///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BaselineDefects.  Find defects in assisted assemblies of 100 Fosmid regions.
// This might be faster if we followed the procedure in RunSome, in 
// paths/long/varcomp/CompareVarsTools.cc.

// results, r48308:
// totals: gaps = 4, endgaps = 0, indels = 69, subs = 131

// Dubious cases:
// fid = 9: gaps = {-76}, endgaps = {}, indels = {11,1}, subs = 23
// fid = 11: gaps = {}, endgaps = {}, indels = {23,16,13,4,3,1}, subs = 23
// fid = 28: gaps = {}, endgaps = {}, indels = {4,4,1}, subs = 1
// fid = 40: gaps = {-78}, endgaps = {}, indels = {12,2}, subs = 10
// fid = 99: gaps = {-49,-49}, endgaps = {}, indels = {1,1,1}, subs = 10
// fid = 106: gaps = {}, endgaps = {}, indels = {364,181}, subs = 3

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency LongProto

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "paths/long/fosmid/Fosmids.h"

int main( )
{
     vec<int> fids = AllFosmids( );
     vec<String> GAPS, ENDGAPS, INDELS;
     int SUBS = 0;
     #pragma omp parallel for
     for ( int i = 0; i < fids.isize( ); i++ )
     {    int fid = fids[i];
          SystemSucceed( "LongProto SAMPLE=human.hpool2 READS=#picard "
               "TMP=tmp.xxx/" + ToString(fid)
               + " IN_GENOME=/wga/dev/references/Homo_sapiens/"
               + "NA12878_Fosmid_Pool.regions.fin/fos." + ToString(fid) + ".fasta "
               + "X=" + ToString(fid) + " LOGGING=EVAL_ASSEMBLY=True "
               + "> xxx" + ToString(fid) );
          vec<String> gaps, endgaps, indels;
          int subs = -1;
          fast_ifstream in( "xxx" + ToString(fid) );
          String line;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "SUMMARY GAPS", 0 ) )
                    Tokenize( line.After( "SUMMARY GAPS: " ), ',', gaps );
               if ( line.Contains( "SUMMARY ENDGAPS", 0 ) )
                    Tokenize( line.After( "SUMMARY ENDGAPS: " ), ',', endgaps );
               if ( line.Contains( "SUMMARY INDELS", 0 ) )
                    Tokenize( line.After( "SUMMARY INDELS: " ), ',', indels );
               if ( line.Contains( "SUMMARY SUBS", 0 ) )
                    subs = line.After( "SUMMARY SUBS: " ).Int( );    }
          #pragma omp critical
          {    cout << "fid = " << fid << ": ";
               cout << "gaps = {" << printSeq(gaps) << "}" << ", ";
               cout << "endgaps = {" << printSeq(endgaps) << "}" << ", ";
               cout << "indels = {" << printSeq(indels) << "}" << ", ";
               cout << "subs = " << subs << endl;
               GAPS.append(gaps), ENDGAPS.append(endgaps), INDELS.append(indels);
               SUBS += subs;    }    }
     cout << "\ntotals: ";
     cout << "gaps = " << GAPS.size( ) << ", endgaps = " << ENDGAPS.size( )
          << ", indels = " << INDELS.size( ) << ", subs = " << SUBS << endl;    }
