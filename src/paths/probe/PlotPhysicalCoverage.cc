///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PlotPhysicalCoverage.  Make a physical coverage plot from read locations.

#include "MainTools.h"
#include "Superb.h"
#include "math/Functions.h"
#include "paths/ReadLoc.h"

int main(int argc, char *argv[])
{
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
     CommandArgument_String_Doc(REGION, 
          "of form s<n>, where n is the id of a scaffold");
     CommandArgument_Int_OrDefault(START, -1);
     CommandArgument_Int_OrDefault(STOP, -1);
     CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
          "if specified, file extension is .READLOCS_PREFIX.readlocs "
          "instead of .readlocs");
     CommandArgument_String_Doc(OUT, "name of png file to generate");
     EndCommandArguments;

     // Define directories.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

     // Define scaffold.

     int s = REGION.After( "s" ).Int( );

     // Load scaffolds.

     int ntigs 
          = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY + ".contigs.fastb" );
     vec<superb> scaffolds;
     ReadSuperbs( sub_dir + "/" + ASSEMBLY + ".superb", scaffolds );
     vec<int> to_super( ntigs, -1 ), to_super_pos( ntigs, -1 );
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    int n = scaffolds[i].Ntigs( );
          for ( int j = 0; j < n; j++ )
          {    to_super[ scaffolds[i].Tig(j) ] = i;
               to_super_pos[ scaffolds[i].Tig(j) ] = j;    }    }

     // Define coordinates on scaffold.

     const superb& S = scaffolds[s];
     const int min_gap = 100;
     int N = 0;
     vec<int> gap( S.Ntigs( ) - 1 );
     for ( int i = 0; i < S.Ntigs( ); i++ )
     {    N += S.Len(i);
          if ( i < S.Ntigs( ) - 1 ) 
          {    gap[i] = Max( min_gap, S.Gap(i) );
               N += gap[i];    }    }

     // Set up.

     vec< vec<int> > cov(3);
     for ( int j = 0; j < cov.isize( ); j++ )
          cov[j].resize( N, 0 );
     String head = sub_dir + "/" + ASSEMBLY;
     if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
     read_locs_on_disk locs_file( head, run_dir );
     String IN = OUT.Before( ".png" ) + ".points";
     ForceAssert( OUT.Contains( ".png", -1 ) );

     // Go through the contigs.

     const int dev_mult = 4.0;
     for ( int i = 0; i < S.Ntigs( ); i++ )
     {    int tig = S.Tig(i);
          vec<read_loc> locs;
          locs_file.LoadContig( tig, locs );
          for ( int z = 0; z < locs.isize( ); z++ )
          {    const read_loc& rl = locs[z];
               int cl = rl.ReadClass( );
               if ( !rl.PartnerPlaced( ) ) continue;
               if ( !rl.Fw( ) || !rl.PartnerRc( ) ) continue;
               int m1 = rl.ContigId( ), m2 = rl.PartnerContigId( );
               int s1 = to_super[m1], s2 = to_super[m2];
               int p1 = to_super_pos[m1], p2 = to_super_pos[m2];
               if ( s1 != s2 ) continue;

               // Find the coordinates of the read and its partner on the scaffold.

               int start1 = rl.Start( ), start2 = rl.PartnerStart( );
               for ( int j = 0; j < p1; j++ )
                    start1 += S.Len(j) + gap[j];
               int stop1 = start1 + rl.ReadLength( );
               for ( int j = 0; j < p2; j++ )
                    start2 += S.Len(j) + gap[j];
               int stop2 = start2 + rl.PartnerReadLength( );

               // Handle the forward case.

               if ( start1 <= stop2 )
               {    int sep = start2 - stop1;
                    if ( Abs( sep - rl.Sep( ) ) <= dev_mult * rl.Dev( ) )
                    {    for ( int j = start1; j < stop2; j++ )
                         {    if ( j >= 0 && j < N ) cov[cl][j]++;    }    }    }

               // Handle the reverse case -- assuming a circular scaffold.

               else
               {    int sep = ( N - stop1 ) + start2;
                    if ( Abs( sep - rl.Sep( ) ) <= dev_mult * rl.Dev( ) )
                    {    for ( int j = start1; j < N; j++ )
                              if ( j >= 0 ) cov[cl][j]++;
                         for ( int j = 0; j < stop2; j++ )
                              if ( j < N ) cov[cl][j]++;    }    }    }    }

     // Plot.

     {    int npasses = 3;
          if ( BigSum( cov[2] ) == 0 ) npasses = 2;
          Ofstream( out, IN );
          for ( int j = 0; j < npasses; j++ )
          {    int start = ( START < 0 ? 0 : START );
               int stop = ( STOP < 0 ? cov[j].isize( ) : STOP );
               for ( int i = start; i < stop; i++ )
               {    out << i << " " << cov[j][i];
                    if ( j == 0 ) out << " 0 0 1\n";
                    if ( j == 1 ) out << " 0 1 0\n";
                    if ( j == 2 ) out << " 1 0 0\n";    }    }    }
     String title = "Scaffold " + ToString(s) + ", physical coverage by frag [blue]"
          + ", jump [green], long jump [red]";
     SystemSucceed( "PlotPoints" + ARGC(IN) + ARGC(OUT) + ARG(CONNECT, True) 
          + " X_AXIS_EXTEND=0 Y_AXIS_EXTEND=0 X_AXIS_OFFSET=0 TITLE=\"" + title
          + "\" TITLE_FONTSIZE=12 COLOR_BY_POINT=True" );    }
