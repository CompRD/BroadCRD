//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/Simulation.h"
#include "simulation/ReadSimulatorSimpleCore.h"
#include "random/Shuffle.h"

void MakeSimReads( long_sim sim, const ref_data& ref,
     vecbasevector& reads, vec<ref_loc>& readlocs )
{    const int num_threads = 32;

     const vecbasevector& G = ref.G;
     const vec<double>& ploidy = ref.ploidy;
     const vec<bool>& is_bv_circular = ref.is_circular;

     vec<double> ploidy_u(ploidy);
     UniqueSort(ploidy_u);
     if ( !ploidy_u.solo( ) )
     {    cout << "WARNING: ploidy appears to be mixed.  At present "
               << "the simulator is not correctly set up to handle this."
               << "endl." << endl;    }
     else sim.COVERAGE *= ploidy_u[0];

     int LEN_SIG = int( round( double(sim.LEN) * sim.LEN_SIG_PC/100.0 ) );
     RandomLength len_rnd_gen( sim.LEN, LEN_SIG, sim.LEN/100, sim.RANDOM_SEED );
     vec<double> COPY_NUMBERS( G.size( ), 1.0 );
     vec<size_t> sim_gid, sim_start, sim_stop;
     vec<Bool> rc;
     for ( size_t ibv = 0; ibv < G.size(); ibv++ ) 
     {    float coverage = COPY_NUMBERS[ibv] * sim.COVERAGE; 
          reads_size_and_position_compute( G[ibv], is_bv_circular[ibv], ibv,
               coverage, len_rnd_gen, sim.RANDOM_SEED, False, 0.0, &reads,
               &sim_gid, &sim_start, &rc );    }
     sim_stop.resize( reads.size( ) );
     for ( size_t id = 0; id < reads.size( ); id++ )
          sim_stop[id] = sim_start[id] + reads[id].size( );
//     rc.resize( reads.size( ) );
     #pragma omp parallel for
     for ( int i_thread = 0; i_thread < num_threads; i_thread++ )
     {    reads_sample_parallel( G, is_bv_circular, sim_gid, sim_start, rc, &reads, 
               sim.RANDOM_SEED, i_thread, num_threads );    }
     #pragma omp parallel for
     for ( int i_thread = 0; i_thread < num_threads; i_thread++ )
     {    reads_perturb_parallel( &reads, 0, sim.ERR_DEL, sim.ERR_INS, sim.ERR_SUB,
               sim.RANDOM_SEED, i_thread, num_threads );    }
     readlocs.resize( sim_gid.size( ) );
     for ( int i = 0; i < sim_gid.isize( ); i++ )
          readlocs[i] = ref_loc( sim_gid[i], sim_start[i], sim_stop[i], rc[i] );

     // Randomize order of reads.

     vec<int> perm;
     Shuffle( reads.size( ), perm );
     PermuteSwappableVec( reads, perm );
     PermuteVec( readlocs, perm );    }

// Note below: should be filtering readlocs too!!

void FilterSimReads( const long_sim& sim,
     const vecbasevector& G, const vec<ref_loc>& readlocs, String& RID )
{
     if ( sim.SIM_FILTER != "" || sim.SIM_END_EXCLUDE > 0 )
     {    vec<int> ids;
          int N = readlocs.size( );
          if ( sim.SIM_FILTER != "" )
          {    vec<String> sf;
               ParseStringSet( sim.SIM_FILTER, sf );
               vec< triple<int,int,int> > sfs;
               for ( int j = 0; j < sf.isize( ); j++ )
               {    int g = sf[j].Before( "." ).Int( );
                    int start = sf[j].Between( ".", "-" ).Int( );
                    int stop = sf[j].After( "-" ).Int( );
                    sfs.push( g, start, stop );    }
               for ( int id = 0; id < N; id++ )
               {    int gstart = readlocs[id].start, gstop = readlocs[id].stop;
                    Bool good = False;
                    for ( int j = 0; j < sfs.isize( ); j++ )
                    {    int g = sfs[j].first;
                         int start = sfs[j].second, stop = sfs[j].third;
                         if ( (int) readlocs[id].id == g 
                              && IntervalOverlap(gstart,gstop,start,stop) > 0 ) 
                         {    good = True;    }    }
                    if (good) ids.push_back(id);    }    }
          else
          {    for ( int id = 0; id < N; id++ )
               {    int gstart = readlocs[id].start, gstop = readlocs[id].stop;
                    if ( gstart < G[ readlocs[id].id ].isize( ) - sim.SIM_END_EXCLUDE
                         && gstop > sim.SIM_END_EXCLUDE )
                    {    ids.push_back(id);    }    }    }
          vec<int> rid;
          if ( RID != "" ) ParseIntSet( RID, rid, true, 0, N );
          else rid = vec<int>( N, vec<int>::IDENTITY );
          ids = Intersection( ids, rid );
          cout << Date( ) << ": after filtering there are "
               << ids.size( ) << " reads" << endl;
          RID = "{";
          for ( int j = 0; j < ids.isize( ); j++ )
               RID += ( j > 0 ? "," : "" ) + ToString( ids[j] );
          RID += "}";    }    }
