// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// ReadGC.  Analyze GC content of reads.  This code was written with the goal of
// identifying and removing bacterial contamination from a particular project.

// Output: for each supercontig, we give the histogram of GC contents observed
// in the reads (one data point per reads).  We also give the cumulative total.
// There are command-line parameters to control which reads and which supers are
// processed.

#include "math/Arith.h"
#include "Basevector.h"
#include "MainTools.h"
#include "ReadLocation.h"
#include "Superb.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_UnsignedInt_OrDefault(HG, 5); // histogram granularity
     CommandArgument_UnsignedInt_OrDefault(MIN_SUPER, 0);
     CommandArgument_UnsignedInt_OrDefault(MAX_SUPER, 2000000000);
     CommandArgument_UnsignedInt_OrDefault(MIN_READ_LENGTH, 1);
     CommandArgument_Double_OrDefault(MIN_MEAN_GC, 0);
     EndCommandArguments;

     // Check arguments.

     if ( 100/HG * HG != 100 )
     {    cout << "HG must divide 100\n";
          exit(1);    }

     // Set up directories.

     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String source_dir = run_dir +
          ( IsDirectory( run_dir + "/preFindSeeds" ) ? "/preFindSeeds" : "" );
     String subdir = run_dir + "/" + SUBDIR;

     // Load assembly data.

     vecbasevector reads( source_dir + "/reads.fastb" );
     int ntigs = MastervecFileObjectCount( subdir + "/mergedcontigs.fastb" );
     READ( subdir + "/mergedcontigs_orig.locs", vec<read_location>, locs );
     int N = reads.size( );
     vec< vec<int> > locs_by_contig(ntigs), locs_by_read(N);
     for ( int i = 0; i < locs.isize( ); i++ )
     {    locs_by_contig[ locs[i].Contig( ) ].push_back(i);
          locs_by_read[ locs[i].ReadId( ) ].push_back(i);    }
     vec<superb> supers;
     ReadSuperbs( subdir + "/mergedcontigs.superb", supers );

     // Look at reads.

     int hg = HG, hi = 100/HG;
     vec<int> histall(hi, 0);
     for ( int i = 0; i < supers.isize( ); i++ )
     {    if ( supers[i].FullLength( ) < (int) MIN_SUPER ) continue;
          if ( supers[i].FullLength( ) > (int) MAX_SUPER ) continue;
          cout << "supercontig " << i << ", length = "
               << supers[i].FullLength( ) << "\n";
          const superb& sup = supers[i];
          vec<int> hist(hi, 0);
          static vec<Float> gcvals;
          gcvals.clear( );
          for ( int j = 0; j < sup.Ntigs( ); j++ )
          {    int m = sup.Tig(j);
               for ( int k = 0; k < locs_by_contig[m].isize( ); k++ )
               {    int id1 = locs[ locs_by_contig[m][k] ].ReadId( );
                    if ( reads[id1].size( ) < MIN_READ_LENGTH ) continue;
                    Float gcp = GcPercent(reads[id1]);
                    ++hist[ Min( hi - 1, int(floor(gcp)) / hg ) ];
                    gcvals.push_back(gcp);    }    }
          Float sum = 0;
          for ( int j = 0; j < gcvals.isize( ); j++ )
               sum += gcvals[j];
          if ( sum / Float( gcvals.size( ) ) < Float(MIN_MEAN_GC) ) continue;
          for ( int j = 0; j < hi; j++ )
          {    histall[j] += hist[j];
               if ( hist[j] == 0 ) continue;
               cout << "GC " << hg * j << "-" << hg * (j+1) << "%: "
                    << hist[j] << " reads\n";    }    }
     cout << "\nSUM FOR ABOVE:\n";
     for ( int j = 0; j < hi; j++ )
     {    if ( histall[j] == 0 ) continue;
          cout << "GC " << hg * j << "-" << hg * (j+1) << "%: "
               << histall[j] << " reads\n";    }    }
