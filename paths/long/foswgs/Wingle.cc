#include "MainTools.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"

int main( )
{    RunTime( );

     // Load shb.

     SupportedHyperBasevector shb;
     BinaryReader::readFile( "alpha.shbv", &shb );

     // Set up logging etc.

     long_logging logc( "" );
     long_heuristics heur( "" );
     ref_data ref;
     vec<ref_loc> readlocs;
     String VERB = "", OUT_INT_HEAD = "";
     long_logging_control log_control( ref, &readlocs, OUT_INT_HEAD, VERB );

     // Heuristics.

     const double min_weight_split = 10.0;
     const double junk_ratio = 10.0;
     const int max_del = 100;

     PRINT( shb.EdgeObjectCount( ) );

     for ( int pass = 1; pass <= 3; pass++ )
     {    shb.RemoveSmallMostlyAcyclicComponents(logc);
          shb.DeleteReverseComplementComponents(logc);
          shb.DeleteLowCoverage( heur, log_control, logc );    
               }

     shb.PullApart( min_weight_split, logc );
     shb.PullApart2( min_weight_split, logc );

     shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
     shb.RemoveSmallMostlyAcyclicComponents(logc);

     shb.Gulp(log_control,logc);
     shb.Ungulp(logc);

     PRINT( shb.EdgeObjectCount( ) );

     logc.USE_GENOME_FOR_DUMP = False;
     shb.DumpFiles( "beta", log_control, logc );    }

