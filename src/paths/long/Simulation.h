///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONGPROTO_SIMULATION_H
#define LONGPROTO_SIMULATION_H

#include "CoreTools.h"

#define SimulationBool( SIM, DEFAULT, EXPLANATION )                           \
     if ( String(#DEFAULT) == "True" ) SIM = True;                            \
     else if ( String(#DEFAULT) == "False" ) SIM = False;                     \
     else FatalErr( "Illegal value for " << #SIM << "." );                    \
     for ( int i = 0; i < sim.isize( ); i++ )                                 \
     {    const String& x = sim[i];                                           \
          if ( !x.Contains( "=" ) )                                           \
          {    FatalErr( "Illegal simulation parameter " << x << "." );    }  \
          if ( x.Before( "=" ) == #SIM )                                      \
          {    used[i] = True;                                                \
               if ( x.After( "=" ) == "True" ) SIM = True;                    \
               else if ( x.After( "=" ) == "False" ) SIM = False;             \
               else                                                           \
               {    FatalErr( "Illegal value for simulation parameter "       \
                         << x.After( "=" ) << "." );    }                     \
               break;    }    }

#define SimulationInt( SIM, DEFAULT, EXPLANATION )                             \
     if ( String(#DEFAULT).IsInt( ) ) SIM = String(#DEFAULT).Int( );           \
     else FatalErr( "Illegal value for " << #SIM << "." );                     \
     for ( int i = 0; i < sim.isize( ); i++ )                                  \
     {    const String& x = sim[i];                                            \
          if ( !x.Contains( "=" ) )                                            \
          {    FatalErr( "Illegal simulation parameter " << x << "." );    }   \
          if ( x.Before( "=" ) == #SIM )                                       \
          {    used[i] = True;                                                 \
               if ( x.After( "=" ).IsInt( ) )                                  \
                    SIM = x.After( "=" ).Int( );                               \
               else                                                            \
               {    FatalErr( "Illegal value for simulation parameter "        \
                         << x.After( "=" ) << "." );    }                      \
               break;    }    }

#define SimulationDouble( SIM, DEFAULT, EXPLANATION )                          \
     if ( String(#DEFAULT).IsDouble( ) ) SIM = String(#DEFAULT).Double( );     \
     else FatalErr( "Illegal value for " << #SIM << "." );                     \
     for ( int i = 0; i < sim.isize( ); i++ )                                  \
     {    const String& x = sim[i];                                            \
          if ( !x.Contains( "=" ) )                                            \
          {    FatalErr( "Illegal simulation parameter " << x << "." );    }   \
          if ( x.Before( "=" ) == #SIM )                                       \
          {    used[i] = True;                                                 \
               if ( x.After( "=" ).IsDouble( ) )                               \
                    SIM = x.After( "=" ).Double( );                            \
               else                                                            \
               {    FatalErr( "Illegal value for simulation parameter "        \
                         << x.After( "=" ) << "." );    }                      \
               break;    }    }

#define SimulationString( SIM, DEFAULT, EXPLANATION )                          \
     SIM = DEFAULT;                                                            \
     for ( int i = 0; i < sim.isize( ); i++ )                                  \
     {    const String& x = sim[i];                                            \
          if ( !x.Contains( "=" ) )                                            \
          {    FatalErr( "Illegal simulation parameter " << x << "." );    }   \
          if ( x.Before( "=" ) == #SIM )                                       \
          {    used[i] = True;                                                 \
               if ( x.After( "=" ).IsInt( ) )                                  \
                    SIM = x.After( "=" ).Int( );                               \
               else                                                            \
               {    FatalErr( "Illegal value for simulation parameter "        \
                         << x.After( "=" ) << "." );    }                      \
               break;    }    }

class long_sim {

     public:

     long_sim( const String& SIM )
     {    vec<String> sim;
          ParseStringSet( SIM, sim );
          vec<Bool> used( sim.size( ), False );

// Definition of simulation parameters:

SimulationInt( RANDOM_SEED, 1252723, "random seed" );

SimulationInt( LEN, 3000, "simulated read length" );

SimulationDouble( LEN_SIG_PC, 10, "simulated read length dev percent" );

SimulationDouble( COVERAGE, 60.0, "simulated read coverage" );

SimulationDouble( ERR_DEL, 0.03, "deletion rate" );

SimulationDouble( ERR_INS, 0.002, "insertion rate" );

SimulationDouble( ERR_SUB, 0.008, "substitution rate" );

SimulationBool( PERF, False, "shorthand for ERR_DEL = ERR_INS = ERR_SUB = 0." );

// Stuff related to simulation:

SimulationBool( CHEAT_MODEL, False, "use simulation "
     "parameters instead of computing model parameters from data" );

SimulationString( SIM_FILTER, "",
     "if specified, generate the full set of simulated reads, but correct "
     "only those touching the given regions {g.start-stop}; coordinates are "
     "relative to the genome as defined by X, etc." );

SimulationInt( SIM_END_EXCLUDE, 0,
     "if specified, generate the full set of simulated reads, but correct "
     "only those not lying within the given number of end bases of a reference "
     "contig" );

     Bool fail = False;
     for ( int i = 0; i < used.isize( ); i++ )
     {    if ( !used[i] )
          {    cout << "\nIllegal heuristic " << sim[i] << "." << endl;
               fail = True;    }    }
     if (fail)
     {    cout << "Abort." << endl;
          _exit(1);    }

     if (PERF) ERR_DEL = ERR_INS = ERR_SUB = 0;

}

     int RANDOM_SEED;
     int LEN;
     double LEN_SIG_PC;
     double COVERAGE;
     double ERR_DEL;
     double ERR_INS;
     double ERR_SUB;
     Bool PERF;
     Bool CHEAT_MODEL;
     String SIM_FILTER;
     int SIM_END_EXCLUDE;

};

#endif
