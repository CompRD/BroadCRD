///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// CallRefTrace.  Call RefTrace on a given assembly and reference.  To generalize
// as needed.

#include <omp.h>

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/RefTrace.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "util/NullOStream.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(ASSEMBLY, "", "fastb file for assembly");
     CommandArgument_String_OrDefault_Doc(HBV_ASSEMBLY, "", 
          ".hbv file for assembly");
     CommandArgument_String_Doc(REF, 
          "fasta or HyperBasevector file for reference sequence");
     CommandArgument_Int_OrDefault_Doc(VERBOSITY, 0, "logging verbosity");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
          "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault_Doc(BEST_GLOBAL_OUT, "", 
          "fastb file for best global alignment");
     CommandArgument_Bool_OrDefault(FIX_BUG, False);

     // Heuristics:

     CommandArgument_String_OrDefault(MAX_OFFSET_DIFF, "");
     CommandArgument_String_OrDefault(MAX_ERROR_RATE, "");
     CommandArgument_String_OrDefault(OFFSET_ADD, "");
     CommandArgument_String_OrDefault(MAX_TWIDDLE, "");
     CommandArgument_String_OrDefault(MIN_GROUP_FRAC, "");
     CommandArgument_String_OrDefault(MIN_GROUP_SAVE, "");

     EndCommandArguments;

     // Parse heuristics.

     vec<int> max_offset_diff;
     vec<double> max_error_rate;
     vec<int> offset_add;
     vec<int> max_twiddle;
     vec<double> min_group_frac;
     vec<int> min_group_save;
     if ( MAX_OFFSET_DIFF.Contains( "," ) && !MAX_OFFSET_DIFF.Contains( "{" ) )
          MAX_OFFSET_DIFF = "{" + MAX_OFFSET_DIFF + "}";
     if ( MAX_ERROR_RATE.Contains( "," ) && !MAX_ERROR_RATE.Contains( "{" ) )
          MAX_ERROR_RATE = "{" + MAX_ERROR_RATE + "}";
     if ( OFFSET_ADD.Contains( "," ) && !OFFSET_ADD.Contains( "{" ) )
          OFFSET_ADD = "{" + OFFSET_ADD + "}";
     if ( MAX_TWIDDLE.Contains( "," ) && !MAX_TWIDDLE.Contains( "{" ) )
          MAX_TWIDDLE = "{" + MAX_TWIDDLE + "}";
     if ( MIN_GROUP_FRAC.Contains( "," ) && !MIN_GROUP_FRAC.Contains( "{" ) )
          MIN_GROUP_FRAC = "{" + MIN_GROUP_FRAC + "}";
     if ( MIN_GROUP_SAVE.Contains( "," ) && !MIN_GROUP_SAVE.Contains( "{" ) )
          MIN_GROUP_SAVE = "{" + MIN_GROUP_SAVE + "}";
     ParseIntSet( MAX_OFFSET_DIFF, max_offset_diff, false );
     ParseDoubleSet( MAX_ERROR_RATE, max_error_rate, false );
     ParseIntSet( OFFSET_ADD, offset_add, false );
     ParseIntSet( MAX_TWIDDLE, max_twiddle, false );
     ParseDoubleSet( MIN_GROUP_FRAC, min_group_frac, false );
     ParseIntSet( MIN_GROUP_SAVE, min_group_save, false );

     // Check arguments.

     if ( ASSEMBLY != "" && HBV_ASSEMBLY != "" )
	  FatalErr("choose ONE of ASSEMBLY and HBV_ASSEMBLY");

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Set up reference.

     HyperBasevector ref;
     vec<HyperBasevector> GH;
     if ( REF.Contains( ".hbv", -1 ) ) BinaryReader::readFile( REF, &ref );
     else if ( REF.Contains( ".fasta", -1 ) )
     {    vecbasevector r;
          FetchReads( r, 0, REF );
          vec<basevector> rv( r.size( ) );
          for ( int i = 0; i < rv.isize( ); i++ )
               rv[i] = r[i];
          const int K = 80;
          ref = HyperBasevector( K, rv );    }
     else
     {    cout << "Unknown file type for REF." << endl;
          exit(1);    }
     GH.push_back(ref);
     vec<bool> is_circular( 1, false );

     // Set up assembly.

     HyperBasevector hb;
     if ( ASSEMBLY != "" ) {
	  vecbasevector A;
	  A.ReadAll(ASSEMBLY);
	  const int K = 80;
	  vec<basevector> a;
	  for (int i = 0; i < (int) A.size(); i++)
	       a.push_back(A[i]);
          // If the following is slow, we'll call RefTrace separately for each case.
	  hb = HyperBasevector(K,a);
     } else if ( HBV_ASSEMBLY != "" ) {
	  BinaryReader::readFile(HBV_ASSEMBLY, &hb);
     } else
	  FatalErr("you must supply either ASSEMBLY or HBV_ASSEMBLY");
     vec<int> inv( hb.EdgeObjectCount( ), -1 );

     SupportedHyperBasevector shb;
     {   // create a dummy supported hyperbasevector
         vec< vec<int> > usu;
         vec<fix64_6> count_fw, count_rc;
         vec< pair< vec<int>, vec<int> > > pairs;
         vec< vec<pair_point> > pair_data;
         vec< vec< pair<fix64_6,int64_t> > > weights_fw_origin, weights_rc_origin;
         IntDistribution read_length_dist;
         double fudge_mult = 1.0;
         shb = SupportedHyperBasevector( hb, inv, usu, count_fw, count_rc, 
                 weights_fw_origin, weights_rc_origin, pairs, pair_data, 
                 0, read_length_dist, fudge_mult );
     }

     // Call RefTrace.

     long_logging logc( "" );
     logc.REFTRACE_VARIANTS_SUMMARY=True;
     int n = std::max( {1, max_offset_diff.isize( ), max_error_rate.isize( ),
          offset_add.isize( ), max_twiddle.isize( ), min_group_frac.isize( ), 
          min_group_save.isize( ) } );
     RefTraceHeuristics rth_def;
     vec<String> results(n);
     vec< pair<int,int> > gp( n, make_pair(-1,-1) ); 
     vec<int> ids( n, vec<int>::IDENTITY );
     for ( int i = 0; i < n; i++ )
     {    int MAX_OFFSET_DIFF, OFFSET_ADD, MAX_TWIDDLE, MIN_GROUP_SAVE;
          double MAX_ERROR_RATE, MIN_GROUP_FRAC;

          if ( max_offset_diff.empty( ) ) MAX_OFFSET_DIFF = rth_def.max_offset_diff;
          else if ( i > max_offset_diff.isize( ) - 1 )
          {    MAX_OFFSET_DIFF = max_offset_diff.back( );    }
          else MAX_OFFSET_DIFF = max_offset_diff[i];

          if ( max_error_rate.empty( ) ) MAX_ERROR_RATE = rth_def.max_error_rate;
          else if ( i > max_error_rate.isize( ) - 1 )
          {    MAX_ERROR_RATE = max_error_rate.back( );    }
          else MAX_ERROR_RATE = max_error_rate[i];

          if ( offset_add.empty( ) ) OFFSET_ADD = rth_def.offset_add;
          else if ( i > offset_add.isize( ) - 1 )
          {    OFFSET_ADD = offset_add.back( );    }
          else OFFSET_ADD = offset_add[i];

          if ( max_twiddle.empty( ) ) MAX_TWIDDLE = rth_def.max_twiddle;
          else if ( i > max_twiddle.isize( ) - 1 )
          {    MAX_TWIDDLE = max_twiddle.back( );    }
          else MAX_TWIDDLE = max_twiddle[i];

          if ( min_group_frac.empty( ) ) MIN_GROUP_FRAC = rth_def.min_group_frac;
          else if ( i > min_group_frac.isize( ) - 1 )
          {    MIN_GROUP_FRAC = min_group_frac.back( );    }
          else MIN_GROUP_FRAC = min_group_frac[i];

          if ( min_group_save.empty( ) ) MIN_GROUP_SAVE = rth_def.min_group_save;
          else if ( i > min_group_save.isize( ) - 1 )
          {    MIN_GROUP_SAVE = min_group_save.back( );    }
          else MIN_GROUP_SAVE = min_group_save[i];

          RefTraceHeuristics rth( MAX_OFFSET_DIFF, MAX_ERROR_RATE, OFFSET_ADD,
               MAX_TWIDDLE, MIN_GROUP_FRAC, MIN_GROUP_SAVE );

          ostringstream out;
          String BEST_GLOBAL_OUT_i = BEST_GLOBAL_OUT;
          if ( BEST_GLOBAL_OUT != "" ) BEST_GLOBAL_OUT_i += "." + ToString(i);

          ref_data ref;
          ref.GH = GH;
          ref.is_circular= is_circular;
          RefTrace( ref, shb, inv, VERBOSITY, logc, out, rth, 
                  BEST_GLOBAL_OUT_i, FIX_BUG);
          results[i] = out.str( );
          String line;
          istringstream in( results[i] );
          int gaps = -1, penalty = -1;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( " total gaps = " ) )
                    gaps = line.Between( " total gaps = ", "," ).Int( );
               if ( line.Contains( " penalty" ) )
                    penalty = line.Before( " penalty" ).Int( );    }
          gp[i] = make_pair( gaps, penalty );    }
     SortSync( gp, ids );
     if ( gp[0].first < 0 )
     {    cout << "\nSomething went wrong, RefTrace must have failed.\n";
          exit(1);    }
     cout << results[ ids[0] ];
     if( BEST_GLOBAL_OUT!="") Cp2( BEST_GLOBAL_OUT + "." + ToString( ids[0] ), BEST_GLOBAL_OUT );
     cout << endl;    }
