///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ErrorCorrectionEfficiency.  For a "micro"-type assembly, determine the fraction 
// of reads that are corrected as a function of the position on the genome.
//
// Output = lines of the form:
//
// tig pos corrected_count all_count frac
//
// where (tig,pos) give coordinates on the reference genome, all_count is the
// count of all reads at that position, corrected_count is the count of reads
// that survive error correction, and frac = corrected_count/all_count.
//
// Note three modes: PRINT_COV, PLOT_COV, PRINT_UNUSED_READS.
//
// In its present form this code won't work if the assembly is from more than one
// scaffold.
//
// THIS CODE MAY BE BUGGY!

#include "Basevector.h"
#include "FeudalMimic.h"
#include "Intvector.h"
#include "MainTools.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "util/ReadTracker.h"

// MakeDepend: dependency PlotPoints

int main( int argc, char *argv[] )
{
  
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_Int_OrDefault(START, -1);
     CommandArgument_Int_OrDefault(STOP, -1);
     CommandArgument_Bool_OrDefault(PRINT_COV, True);
     CommandArgument_String_OrDefault_Doc(PLOT_COV, "", 
          "name of optional plot file, which must end in .png");
     CommandArgument_Bool_Abbr_OrDefault(PRINT_UNCORRECTED_READS, PUR, False);
     EndCommandArguments;

     // Define directories.
  
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     // Define region.

     String region = FirstLineOfFile( data_dir + "/../region" );
     int tig_id = region.Before( ":" ).Int( );
     int origin = region.Between( ":", "-" ).Int( );

     // Define reads files.

     String reads_orig_head = run_dir + "/frag_reads_orig";
     String reads_filt_head = run_dir + "/frag_reads_filt";
     size_t nreads_orig = MastervecFileObjectCount( reads_orig_head + ".fastb" );
     size_t nreads_filt = MastervecFileObjectCount( reads_filt_head + ".fastb" );

     // Define mapping from _orig to _filt reads.

     vec<int> id_filt( nreads_orig, -1 );
     ReadTracker rt, rtnext;
     rt.Load(reads_filt_head);
     for ( size_t r = 0; r < nreads_filt; r++ )
     {    String source = rt.GetReadSource(r);
          rtnext.Load(source);
          uint32_t index = rt.GetReadIndex(r);
          id_filt[index] = r;    }

     // Load alignments.

     vec<look_align> aligns;
     LoadLookAligns( data_dir + "/frag_reads_orig.qltout", aligns );

     // Load keep.

     vec<Bool> keep;
     BinaryReader::readFile( run_dir + "/frag_reads_edit.keep", &keep );

     // Set up coverage tracking.

     vecbasevector genome( data_dir + "/genome.fastb" );
     VecIntVec cov, ecov;
     Mimic( genome, cov ), Mimic( genome, ecov );

     // Define start and stop positions.

     int start = 0, stop = genome[0].size( );
     if ( START >= 0 ) start = START;
     if ( STOP >= 0 ) stop = STOP;

     // Print uncorrected reads.

     if (PRINT_UNCORRECTED_READS)
     {    for ( int i = 0; i < aligns.isize( ); i++ )
          {    const look_align& la = aligns[i];
               int id = id_filt[la.query_id];
               if ( id < 0 || keep[id] ) continue;
               if ( IntervalOverlap( start, stop, la.pos2( ) - origin,
                    la.Pos2( ) - origin ) > 0 )
               {    cout << id << "\n";    }    }
          if ( !PRINT_COV && PLOT_COV == "" ) return 0;    }

     // Compute coverage.

     for ( int i = 0; i < aligns.isize( ); i++ )
     {    const look_align& la = aligns[i];
          int id = id_filt[la.query_id];
          if ( id < 0 ) continue;
          int t = la.target_id - tig_id;
          for ( int j = Max( 0, la.pos2( ) - origin ); 
               j < Min( genome[t].isize( ), la.Pos2( ) - origin ); j++ )
          {    cov[t][j] = cov[t][j] + 1;
               if ( keep[id] ) ecov[t][j] = ecov[t][j] + 1;    }    }

     // Report results.

     if ( PRINT_COV || PLOT_COV != "" )
     {    ofstream* outp = 0;
          String PLOT_IN;
          if ( PLOT_COV != "" )
          {    outp = new ofstream;
               PLOT_IN = PLOT_COV.Before( ".png" ) + ".txt";
               outp->open( PLOT_IN.c_str( ), ios::app | ios::out );    }
          for ( size_t t = 0; t < genome.size( ); t++ )
          {    for ( int pos = start; pos < stop; pos++ )
               {    double efrac = double( ecov[t][pos] ) / double( cov[t][pos] );
                    if (PRINT_COV)
                    {    cout << t << " " << pos << " " << ecov[t][pos] << " "
                              << cov[t][pos] << " " << efrac << "\n";    }
                    if ( PLOT_COV != "" && cov[t][pos] > 0 )
                         (*outp) << pos << " " << efrac << "\n";    }    }
          if ( PLOT_COV != "" )
          {    outp->close( );
               SystemSucceed( "PlotPoints " + ARG(IN, PLOT_IN) + ARG(OUT, PLOT_COV)
                    + ARG(POINTSIZE, 0.1) + ARG(CONNECT, True) + ARG(PR, 1)
                    + ARG(X_AXIS_EXTEND, 0) + ARG(Y_AXIS_EXTEND, 0) + ARG(MAX_Y, 1)
                    + ARG(PR, 1) );
               Remove(PLOT_IN);    }    }    }
     
