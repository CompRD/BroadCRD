///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// DiscoMerge.  Merge DISCOVAR de novo assemblies.
//
// Deliberately does not create a data directory.

// TO DO:
// - try running computing/loading in parallel
// - can bubbles be merged?

// Shares a lot with AddNewStuffTest in GapToyTools4.cc.
// Copies a bunch of command-line args from GapToyCore.cc.

// INPUT          OPTIONS                N50     MIS   EDGES  SAVED  HOURS
//                                       CONTIG        /KB    AS
// ---------------------------------------------------------------------------
// 50946.F1,51472.F2,51476.F3             72.6   20   14.7    B5     26
// 50946.F1,51472.F2,51476.F3,50946.CEPH,51452.YRI -- DIED    B6

//
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "kmers/BigKPather.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/FinalFiles.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/ImprovePath.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/MakeGaps.h"
#include "paths/long/large/Simplify.h"
#include "system/HostName.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(INSTANCES, "comma-separated list of instances");
     CommandArgument_String_Doc(INST_OUT, "output instance");
     CommandArgument_Bool_OrDefault_Doc(NAME2, False,
          "name assemblies using INSTANCE.After(.)");
     CommandArgument_String_OrDefault_Doc(NAMES, "",
          "optional comma-separated list of names for samples; only makes "
          "sense if each assembly has a single sample");
     CommandArgument_String_OrDefault_Doc(IN_DIR, "a.final",
          "could be a.patched or a.fin instead");
     CommandArgument_Bool_OrDefault_Doc(CREATE_DATA_DIR, False,
          "create data directory; uses lots of memory and disk");
     CommandArgument_String_OrDefault_Doc(USER, "",
          "use this instead of Getenv(USER)");
     EndCommandArguments;

     // Mess with args.

     vec<String> names;
     if ( NAMES != "" ) ParseStringSet( "{" + NAMES + "}", names );

     // Make output directory.

     String ROOT = "/wga/scr4";
     String user = ( USER != "" ? USER : Getenv("USER","noname") );
     String out_work_dir = ROOT + "/" + user + "/GapToy/" + INST_OUT;
     Mkdir777(out_work_dir);

     // Echo command.

     {    string hostname = getHostName();
          hostname = hostname.substr( 0, hostname.find('.') );
          Ofstream( out, out_work_dir + "/the_command" );
          out << "\n" << hostname << ": " << command.TheCommand( ) << endl;    }

     // Find instances.

     vec<String> instances;
     ParseStringSet( "{" + INSTANCES + "}", instances );
     ForceAssert( instances.nonempty( ) );
     if (NAME2)
     {    for ( int i = 0; i < instances.isize( ); i++ )
          {    if ( !instances[i].Contains( "." ) )
               {    cout << "To use NAME2, each instance must contain a dot."
                         << endl;
                    Scram(1);    }    }    }
     if ( instances.nonempty( ) && names.size( ) != instances.size( ) )
     {    cout << "INSTANCES and NAMES must have the same length." << endl;
          Scram(1);    }

     // Load assemblies.

     String dir0 = ROOT + "/" + user + "/GapToy/" + instances[0];
     if ( IsRegularFile( dir0 + "/genome.fastb" ) )
          Cp2( dir0 + "/genome.fastb", out_work_dir );
     if ( IsRegularFile( dir0 + "/genome.names" ) )
          Cp2( dir0 + "/genome.names", out_work_dir );
     const int K = 200;
     HyperBasevector hb(K);
     ReadPathVec paths2;
     {    vec<HyperBasevector> HB( instances.size( ) );
          for ( int i = 0; i < instances.isize( ); i++ )
          {    cout << Date( ) << ": loading assembly " << i+1 << " of "
                    << instances.size( ) << endl;
               String work_dir = ROOT + "/" + user + "/GapToy/" + instances[i];
               BinaryReader::readFile( work_dir + "/" + IN_DIR + "/a.hbv", &HB[i] );
               ForceAssertEq( K, HB[i].K( ) );    }
          cout << Date( ) << ": merging assemblies" << endl;
          hb.SetToDisjointUnionOf(HB);    }

     // Merge HyperBasevectors.

     {    cout << Date( ) << ": building allx, peak mem = " 
               << PeakMemUsageGBString( ) << endl;
          vecbasevector allx;
          BuildAll( allx, hb );
          HyperBasevector hb3;
          {    ReadPathVec allx_paths;
               const int coverage = 4;
               cout << Date( ) << ": building big hbv, peak mem = " 
                    << PeakMemUsageGBString( ) << endl;
               buildBigKHBVFromReads( K, allx, coverage, &hb3, &allx_paths );    }
          hb = hb3;    }

     // Get the involution.

     cout << Date( ) << ": finding involution, peak mem = "
          << PeakMemUsageGBString( ) << endl;
     vec<int> inv2;
     hb.Involution(inv2);
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // Create data directory.

     if (CREATE_DATA_DIR) 
     {    Mkdir777( out_work_dir + "/data" );
          vecbasevector bases;
          VecPQVec quals;
          for ( int i = 0; i < instances.isize( ); i++ )
          {    String work_dir = ROOT + "/" + user + "/GapToy/" + instances[i];
               bases.ReadAll( work_dir + "/data/frag_reads_orig.fastb", True );
               quals.ReadAll( work_dir + "/data/frag_reads_orig.qualp", True );    }
          bases.WriteAll( out_work_dir + "/data/frag_reads_orig.fastb" );
          quals.WriteAll( out_work_dir + "/data/frag_reads_orig.qualp" );    }

     // Compute paths from scratch.

     {    cout << Date( ) << ": recomputing paths" << endl;
          const int64_t batch = 100000000;
          int64_t total = 0;
          double lclock = 0, cclock = 0;
          int64_t total_bytes = 0;
          for ( int i = 0; i < instances.isize( ); i++ )
          {    cout << Date( ) << ": processing instance " << i+1 
                    << " of " << instances.size( ) << endl;
               String work_dir = ROOT + "/" + user + "/GapToy/" + instances[i];
               int64_t N = MastervecFileObjectCount(
                    work_dir + "/data/frag_reads_orig.fastb" );
               total_bytes += FileSize( work_dir + "/data/frag_reads_orig.fastb" );
               total_bytes += FileSize( work_dir + "/data/frag_reads_orig.qualp" );
               for ( int64_t start = 0; start < N; start += batch )
               {    int64_t stop = Min( start + batch, N );
                    vecbasevector bases;
                    VecPQVec quals;

                    // The vast majority of the time in what follows is spent 
                    // reading from files, so not much could be saved by running
                    // it and the path extension in parallel.  (Old comment, perhaps
                    // still valid.)

                    lclock -= WallClockTime( );
                    bases.ReadRange( 
                         work_dir + "/data/frag_reads_orig.fastb", start, stop );
                    quals.ReadRange( 
                         work_dir + "/data/frag_reads_orig.qualp", start, stop );
                    lclock += WallClockTime( );
                    cclock -= WallClockTime( );
                    ReadPathVec paths2i( bases.size( ) );
                    path_improver pimp;
                    vec<int64_t> ids;
                    ImprovePaths( paths2i, hb, inv2, bases, quals, ids,
                         pimp, False, False );
                    const int EXT_MODE = 1;
                    const int MIN_GAIN = 5;
                    Bool extend_paths_verbose = False;
                    #pragma omp parallel for
                    for ( int64_t id = 0; id < (int64_t) paths2i.size( ); id++ )
                    {    if ( paths2i[id].size( ) > 0 ) paths2i[id].resize(1);
                         ExtendPath( paths2i[id], id, hb, to_right, bases[id], 
                              quals[id], MIN_GAIN, extend_paths_verbose, 
                              EXT_MODE );    }
                    paths2.Append(paths2i);
                    cclock += WallClockTime( );    }
               total += N;    }
          cout << ConvertTime(cclock) << " used computing" << endl;
          cout << ConvertTime(lclock) << " used loading" << endl;
          cout << "read rate = " << setprecision(3) 
               << ( total_bytes/1000000.0 ) / lclock << " MB/second" << endl;
          cout << Date( ) << ": validating" << endl;
          Validate( hb, inv2, paths2 );    }

     // Compute subsam names and subsam_starts;

     cout << Date( ) << ": computing subsams" << endl;
     vec<String> subsam_names;
     for ( int i = 0; i < instances.isize( ); i++ )
     {    String work_dir = ROOT + "/" + user + "/GapToy/" + instances[i];
          vec<String> sn;
          if ( NAMES != "" ) sn.push_back( names[i] );
          else
          {    BinaryReader::readFile( work_dir + "/subsam.names", &sn );
               for ( int j = 0; j < sn.isize( ); j++ )
               {    String head = ToString(i+1);
                    if (NAME2) head = instances[i].After( "." );
                    sn[j] = head + "." + sn[j];    }    }
          subsam_names.append(sn);    }
     BinaryWriter::writeFile( out_work_dir + "/subsam.names", subsam_names );
     vec<int64_t> subsam_starts;
     int64_t all_reads = 0;
     for ( int i = 0; i < instances.isize( ); i++ )
     {    String work_dir = ROOT + "/" + user + "/GapToy/" + instances[i];
          vec<int64_t> s;
          BinaryReader::readFile( work_dir + "/subsam.starts", &s );
          for ( int j = 0; j < s.isize( ); j++ )
               s[j] += all_reads;
          subsam_starts.append(s);
          all_reads += MastervecFileObjectCount(
               work_dir + "/data/frag_reads_orig.fastb" );    }
     cout << "subsam_names: " << printSeq(subsam_names) << endl;
     cout << "subsam_starts: " << printSeq(subsam_starts) << endl;
     BinaryWriter::writeFile( out_work_dir + "/subsam.starts", subsam_starts );

     // Do some minimal simplification.

     RemoveHangs( hb, inv2, paths2, 700 );
     RemoveSmallComponents3(hb);
     Cleanup( hb, inv2, paths2 );
     Validate( hb, inv2, paths2 );

     // Directories.

     String fin_dir = out_work_dir + "/a.fin";
     Mkdir777(fin_dir);

     // OLD SIMPLIFICATION -- WOULD REQUIRE ALL BASES AND QUALS.

     /*
     {
          // Run PartnersToEnds.

          cout << Date( ) << ": running PartnersToEnds" << endl;
          PartnersToEnds( hb, paths2, bases, quals );
          TestInvolution( hb, inv2 );
          Validate( hb, inv2, paths2 );

          // Simplify the assembly.

          if (SIMPLIFY)
          {    cout << Date( ) << ": simplify" << endl;
               const int MAX_SUPP_DEL = 0;
               const Bool TAMP_EARLY = True;
               const int MIN_RATIO2 = 8;
               const int MAX_DEL2 = 200;
               const Bool PLACE_PARTNERS = False;
               const Bool DEGLOOP = True;
               const Bool EXT_FINAL = True;
               const int EXT_FINAL_MODE = 1;
               const int DEGLOOP_MODE = 1;
               const double DEGLOOP_MIN_DIST = 2.5;
               const Bool IMPROVE_PATHS = True;
               const Bool IMPROVE_PATHS_LARGE = False;
               const Bool FINAL_TINY = True;
               const Bool UNWIND3 = True;
               Simplify( fin_dir, hb, inv2, paths2, bases, quals, MAX_SUPP_DEL, 
                    TAMP_EARLY, MIN_RATIO2, MAX_DEL2, PLACE_PARTNERS, False, "", 
                    DEGLOOP, EXT_FINAL, EXT_FINAL_MODE, False, vec<int>( ), 
                    DEGLOOP_MODE, DEGLOOP_MIN_DIST, IMPROVE_PATHS, 
                    IMPROVE_PATHS_LARGE, FINAL_TINY, UNWIND3 );
               TestInvolution( hb, inv2 );    }    }
     */

     // Find lines and write files.

     cout << Date( ) << ": finding lines" << endl;
     const int MAX_CELL_PATHS = 50;
     const int MAX_DEPTH = 10;
     vec<vec<vec<vec<int>>>> lines;
     FindLines( hb, inv2, lines, MAX_CELL_PATHS, MAX_DEPTH );
     BinaryWriter::writeFile( fin_dir + "/a.lines", lines );
     {    vec<int> llens, npairs;
	  GetLineLengths( hb, lines, llens );
	  GetLineNpairs( hb, inv2, paths2, lines, npairs );
	  BinaryWriter::writeFile( fin_dir + "/a.lines.npairs", npairs );
	  vec<vec<covcount>> covs;
	  ComputeCoverage( hb, inv2, paths2, lines, subsam_starts, covs );
	  WriteLineStats( fin_dir + "/a", lines, llens, npairs, covs ); 
	  double cn_frac_good = CNIntegerFraction(hb, covs);
          cout << "CN fraction good = " << cn_frac_good << endl;    }

     // Scaffold.

     const int MIN_LINE = 5000;
     const int MIN_LINK_COUNT = 3;
     const Bool GAP_CLEANUP = True;
     const String final_dir = out_work_dir + "/a.final";
     Mkdir777(final_dir);
     {    VecULongVec invPaths;
          invert( paths2, invPaths, hb.EdgeObjectCount( ) );
          MakeGaps( hb, inv2, paths2, invPaths, MIN_LINE, MIN_LINK_COUNT, 
               out_work_dir, "fin", False, GAP_CLEANUP );    }

     // Carry out final analyses and write final assembly files.

     const String EVALUATE = "False";
     const Bool EVALUATE_VERBOSE = False;
     const Bool ALIGN_TO_GENOME = True;
     const String X;
     map<String,GapToyResults> res;
     const String SAMPLE;
     const String species;
     const vec<int> fosmids;
     const vecbasevector G;
     Bool SAVE_FASTA = False; // ???????????????????????????????????????????????????
     FinalFiles( hb, inv2, paths2, subsam_names, subsam_starts, out_work_dir, 
          final_dir, MAX_CELL_PATHS, MAX_DEPTH, ALIGN_TO_GENOME, EVALUATE, 
          EVALUATE_VERBOSE, X, res, SAMPLE, species, fosmids, G, SAVE_FASTA );

     // Done.

     cout << Date( ) << ": done" << endl;
     cout << "peak mem usage = " << PeakMemUsageGB( ) << " GB" << endl;
     Scram(0);    }
