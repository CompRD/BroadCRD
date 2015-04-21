///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: private
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency SAM2CRDDump

#include <functional>
#include <algorithm>

#include "Basevector.h"
#include "CoreTools.h"
#include "FastIfstream.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "TokenizeString.h"
#include "efasta/EfastaTools.h"
#include "lookup/LookAlign.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/DataSpec.h"
#include "paths/long/DiscovarTools.h"
#include "paths/long/LoadAndCorrect.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/Satellite.h"
#include "paths/long/fosmid/FosmidPool.h"
#include "paths/long/fosmid/FosmidVector.h"
#include "paths/long/fosmid/Fosmids.h"
#include "paths/long/local/Setup.h"
#include "util/Fastq.h"

void LoadPacBio( const String& X, const String& TMP, vecbasevector& pb,
     const long_logging_control& log_control, const long_logging& logc,
     const Bool PACBIO_MORE )
{    double clock = WallClockTime( );
     String pbdir = "/wga/scr4/picard/pacbio";
     vec<String> runs;
     if ( !PACBIO_MORE ) runs = { "019914", "019915", "019916" };
     else
     {    runs = { "019892", "019891", "019899", "019900",
               "019914", "019915", "019916" };    }
     for ( int r = 0; r < runs.isize( ); r++ )
     {    String region;
          vec<String> regions;
          vec< vec< pair<String,String> > > junctions, breaks, edits;
          ParseFosmidPoolMetainfo( regions, junctions, breaks, edits );
          vec<int> xs;
          int status;
          ParseIntSet( X, xs, status );
          for ( int j = 0; j < xs.isize( ); j++ )
          {    if ( j > 0 ) region += " ";
               region += regions[ xs[j] ];    }
          String getsam = "samtools view -h " + pbdir + "/" + runs[r]
               + "/comprd-data/filtered_subreads.bam " + region;
          String dump_tail = " SEP=-15 DEV=12 NOMINAL_READ_LEN=251 USE_OQ=True "
               "NH=True LOG_TO_CERR=False WRITE_ALIGNS=False WRITE_NAMES=False "
               "LIBINFO=/wga/scr4/dexter/libinfo/dexter_libs "
               " > " + TMP + "/pb.SAM2CRDDump.log";
          status = System( getsam + " | SAM2CRDDump OUT_HEAD=" + TMP + "/pb."
               + ToString(r) + dump_tail );
          if ( status != 0 ) DiscovarTools::ExitSamtoolsFailed( );
          pb.ReadAll( TMP + "/pb." + ToString(r) + ".fastb", True );     }
     vec<int> lens;
     for ( int i = 0; i < (int) pb.size( ); i++ )
          lens.push_back( pb[i].size( ) );
     Sort(lens);
     int nref = (*(log_control.G))[0].size( );
     cout << Date( ) << ": " << pb.size( )
          << " PacBio reads of median length " << lens[ lens.size( )/2 ] 
          << " = " << setiosflags(ios::fixed) << setprecision(1)
          << Sum(lens) / double(nref) << resetiosflags(ios::fixed)
          << "X coverage" << endl;
     REPORT_TIME( clock, "used loading PacBio data" );    }

void LoadIlluminaReads( const String& SAMPLE, const vec<picard_spec>& specs,
     const String& region, const String& TMP, const long_logging& logc, 
     const Bool USE_PF_ONLY, const Bool KEEP_NAMES )
{
     // For SAMPLE = hpool2 or hpool3, load clone end sequences.

     vec<String> endseqs, endseq_ids;
     if ( ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) && region != "" )
     {    // Note that the following requires that the program EndSeq has been run.
          fast_ifstream in( "/wga/scr4/human_data/NA12878_A2925/endseqs.sam" );
          String line;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               endseqs.push_back(line);
               endseq_ids.push_back( line.Before( "\t" ) );    }    }
     UniqueSort(endseq_ids);

     // For SAMPLE = hpool2 or hpool3, determine which end sequences are present in
     // the aligned reads.

     vec<Bool> have_end( endseq_ids.size( ), False );
     vec<String> posts3;
     vec<char> sep;
     sep.push_back( '\t' );
     if ( ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) && region != "" )
     {    basevector s1( "ACACGACGCTCTTCCGATCTAGTTGCTT" );
          basevector s2( "GACGTGTGCTCTTCCGATCTAGTTGCTT" );
          String s1f = s1.ToString( ), s1r, s2f = s2.ToString( ), s2r;
          StringReverseComplement( s1f, s1r );
          StringReverseComplement( s2f, s2r );
          vec<String> posts;
          for ( int i = 0; i < specs.isize( ); i++ )
          {    const picard_spec& s = specs[i];
               String getsam;
               if ( s.special != "" )
                    getsam = "samtools view -h " + s.special + " " + region;
               else
               {    String picard = ( !s.direct ? "/wga/scr4/picard" 
                         : "/seq/picard" );
                    getsam = "samtools view -h " + picard + "/" + s.flowcell + "/"
                         + s.picard_run + "/" + "/" + ToString(s.lane) + "/" + s.lib
                         + "/" + s.flowcell + "." + ToString(s.lane)
                         + ".aligned.duplicates_marked.bam " + region;    }
               fast_pipe_ifstream in(getsam);
               String line;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    String id = line.Before( "\t" );
                    int p = BinPosition( endseq_ids, id );
                    if ( p >= 0 )
                    {    have_end[p] = True;
                         vec<String> fields;
                         TokenizeStrictly( line, sep, fields );
                         const String& s = fields[9];
                         if ( line.Contains(s1f) )
                              posts.push_back( s1f + s.After(s1f) );
                         if ( line.Contains(s2f) )
                              posts.push_back( s2f + s.After(s2f) );
                         if ( line.Contains(s1r) )
                         {    String x = s.Before(s1r), y;
                              StringReverseComplement( x, y );
                              posts.push_back( s1f + y );    }
                         if ( line.Contains(s2r) )
                         {    String x = s.Before(s2r), y;
                              StringReverseComplement( x, y );
                              posts.push_back( s2f + y );    }    }    }    }
          vec<String> posts2;
          const int mlen = 48;
          for ( String& x : posts )
          {    if ( x.isize( ) < mlen ) continue;
               x.resize(mlen);
               posts2.push_back(x);    }
          Sort(posts2);
          vec< pair<int,String> > px;
          for ( int i = 0; i < posts2.isize( ); i++ )
          {    int j = posts2.NextDiff(i);
               int n = j - i;
               px.push( n, posts2[i] );
               i = j - 1;    }
          ReverseSort(px);
          for ( int i = 0; i < Min( 2, px.isize( ) ); i++ )
          {    posts3.push_back( px[i].second );
               String s;
               StringReverseComplement( px[i].second, s );
               posts3.push_back(s);    }    }

     // Now determine which end sequences should be included.

     vec<Bool> extra( endseq_ids.size( ), False );
     if ( ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) && region != "" )
     {    for ( const String& x : endseqs )
          {    vec<String> fields;
               TokenizeStrictly( x, sep, fields );
               Bool found = False;
               for ( const String& y : posts3 )
               {    if ( fields[9].Contains(y) )
                    {    String id = x.Before( "\t" );
                         int p = BinPosition( endseq_ids, id );
                         if ( !have_end[p] ) extra[p] = True;    }    }    }    }

     // Now do the main work.

     // Define primary datasets.
     vec<String> heads;

     if( specs.size() > 0 && specs.front().special.Contains(".fastq",String::npos)){
         const char q_offset=33;
         if(specs.size()==1){
             fastq::Fastq2FastbQualbPair(specs.front().special,q_offset,TMP+"/0");
              heads.push_back( "0");
         }
         else{
             ForceAssert( specs.size()%2==0);
             for(size_t ii=0;ii*2<specs.size();++ii){
                 const String& r1 = specs[2*ii].special;
                 const String& r2 = specs[2*ii+1].special;
                 ForceAssert( r1.Contains(".fastq",String::npos));
                 ForceAssert( r2.Contains(".fastq",String::npos));
                 fastq::Fastq2FastbQualbPair_Interleave(r1,r2,q_offset,TMP+"/"+ToString(ii));
                 heads.push_back( ToString(ii) );
             }
         }
     }
     else{
         for ( int i = 0; i < specs.isize( ); i++ )
         {    const picard_spec& s = specs[i];
              String getsam;
              if ( s.special != "" )
              {    String bam = s.special;
                   if (logc.PRINT_BAMS) cout << Date( ) << ": loading " << bam << endl;
                   getsam = "samtools view -h " + bam + " " + region;    }
              else
              {    String picard = ( !s.direct ? "/wga/scr4/picard" : "/seq/picard" );
                   String bam = picard + "/" + s.flowcell + "/"
                        + s.picard_run + "/" + "/" + ToString(s.lane) + "/" + s.lib + "/"
                        + s.flowcell + "." + ToString(s.lane)
                        + ".aligned.duplicates_marked.bam";
                   if (logc.PRINT_BAMS) cout << Date( ) << ": loading " << bam << endl;
                   getsam = "samtools view -h " + bam + " " + region;    }
              SamIAm( i, getsam, TMP, logc.KEEP_LOCS,
                   "/wga/scr4/dexter/libinfo/dexter_libs", USE_PF_ONLY,
                   KEEP_NAMES );
              heads.push_back( ToString(i) );
         }
     }



     // Harvest the extra reads.

     if ( ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) && region != "" )
     {    vec<String> sam;
          for ( int i = 0; i < specs.isize( ); i++ )
          {    const picard_spec& s = specs[i];

               String getsam_header;
               if ( s.special != "" )
                    getsam_header = "samtools view -H " + s.special;
               else
               {    String picard = ( !s.direct ? "/wga/scr4/picard" : "/seq/picard" );
                    getsam_header = "samtools view -H " + picard + "/" + s.flowcell
                         + "/" + s.picard_run + "/" + "/" + ToString(s.lane) + "/"
                         + s.lib + "/" + s.flowcell + "." + ToString(s.lane)
                         + ".aligned.duplicates_marked.bam";    }
               fast_pipe_ifstream in(getsam_header);
               String line;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    sam.push_back(line);    }    }
          for ( const String& s : endseqs )
          {    int p = BinPosition( endseq_ids, s.Before( "\t" ) );
               if ( extra[p] ) sam.push_back(s);    }
          {    Ofstream( out, TMP + "/ends.sam" );
               for ( const String& s : sam )
                    out << s << "\n";    }
          String dump_tail =
                  " SEP=-15 DEV=12 NOMINAL_READ_LEN=251 USE_OQ=True NH=True "
                  "LOG_TO_CERR=False WRITE_NAMES=False"
                  " LIBINFO=/wga/scr4/dexter/libinfo/dexter_libs "
                  "WRITE_ALIGNS="
               + String( logc.KEEP_LOCS ? "True" : "False" )
               + " > " + TMP + "/SAM2CRDDump.log";
          int status = System( "cat " + TMP + "/ends.sam | SAM2CRDDump "
               "OUT_HEAD=" + TMP + "/ends" + dump_tail );
          if ( status != 0 )
          {    cout << "\nSystem call to samtools/SAM2CRDDump failed.\n\n";
               cout << "Here is the output:\n";
               CpAppend( TMP + "/SAM2CRDDump.log", cout );
               cout << "\nAborting.\n" << endl;
               TracebackThisProcess( );    }
          heads.push_back( "ends" );    }

     // Merge readsets.

     MergeReadSets( heads, TMP, logc, KEEP_NAMES );    }

void LoadAndCorrectIllumina( const String& SAMPLE, const String& READS,
     const String& DATASET, const String& X, const String& TMP,
     const long_data_spec& spec, VecEFasta& corrected, vec<int>& cid,
     vec<pairing_info>& cpartner, const Bool HAVE_PICARD,
     const long_heuristics& heur, const long_logging_control& log_control,
     const long_logging& logc, const int NUM_THREADS, const String& EXIT,
     const double clock, bool useOldLRPMethod, const Bool KEEP_NAMES )
{
     // General setup.

     double clock1 = WallClockTime( );
     double SELECT_FRAC = spec.SELECT_FRAC;
     double COVERAGE_CAP = spec.COVERAGE_CAP;
     String HUMAN_CONTROLS = spec.HUMAN_CONTROLS;
     if ( READS == "#picard" )
     {    if ( SAMPLE != "rhody" && SAMPLE != "ecoli" && SAMPLE != "ecoli11" 
               && SAMPLE != "plasmo" && SAMPLE != "human" && SAMPLE != "hpool1" 
               && SAMPLE != "hpool2" && SAMPLE != "human.hpool2"
               && SAMPLE != "hpool3" && SAMPLE != "human.hpool3" 
               && SAMPLE != "rhino"
               && SAMPLE != "entero" && SAMPLE != "bifi" && SAMPLE != "scardovia"
               && SAMPLE != "tb148" && SAMPLE != "ecoli12" && SAMPLE != "ecoli_scs"
               && SAMPLE != "bcereus" && SAMPLE != "tb" && SAMPLE != "arabidopsis" )
          {    FAIL_MSG( "If READS=#picard, at present, SAMPLE must be rhody, "
                    "bcereus, ecoli, ecoli11, ecoli12, ecoli_scs, entero, tb, "
                    "tb148, plasmo, arabidopsis, "
                    "bifi, scardovia, rhino, human, hpool1, hpool2, hpool3, "
                    "human.hpool2 or human.hpool3." );    }    }
     Mkdir777(TMP);
     vec<picard_spec> specs;
     String region;
     vec<String> chrlist, regions;
     ParseRegions( X, regions );

     // If human controls are to be used, append to 'regions' and remove overlaps.

     if ( HUMAN_CONTROLS != "" )
     {    if ( !isdigit( HUMAN_CONTROLS[0] ) )
          {    ostringstream out;
               out << printSeq( expand_fosmids(HUMAN_CONTROLS) );
               HUMAN_CONTROLS = out.str( );    }
          vec<String> regions0;
          vec< vec< pair<String,String> > > junctions, breaks, edits;
          ParseFosmidPoolMetainfo( regions0, junctions, breaks, edits );
          vec<String> all = GetFinishedFosmidFiles( );
          vec<int> hc;
          if ( HUMAN_CONTROLS != "all" ) ParseIntSet( HUMAN_CONTROLS, hc );
          vec<Bool> to_delete0( regions0.size( ), True );
          for ( const String& f : all )
          {    int id = f.Between( ".", "." ).Int( );
               if ( HUMAN_CONTROLS != "all" && !BinMember( hc, id ) ) continue;
               to_delete0[ f.Between( ".", "." ).Int( ) ] = False;    }
          EraseIf( regions0, to_delete0 );
          for ( int i = 0; i < regions0.isize( ); i++ )
          {    String id = regions0[i].Before( ":" );
               String range = regions0[i].After( ":" );
               if ( range.Contains( " " ) ) range = range.Before( " " );
               int start = range.Before( "-" ).Int( ) - ExtendRegionBy( );
               int stop = range.After( "-" ).Int( ) + ExtendRegionBy( );
               range = ToString(start) + "-" + ToString(stop);
               regions.push_back( id + ":" + range );    }
          vec<Bool> to_delete( regions.size( ), False );
          combine_regions:
          for ( int i1 = 0; i1 < regions.isize( ); i1++ )
          for ( int i2 = 0; i2 < regions.isize( ); i2++ )
          {    if ( i2 == i1 || to_delete[i1] || to_delete[i2] ) continue;
               String &r1 = regions[i1], &r2 = regions[i2];
               if ( r1.Before( ":" ) != r2.Before( ":" ) ) continue;
               int start1 = r1.Between( ":", "-" ).Int( );
               int stop1 = r1.After( "-" ).Int( );
               int start2 = r2.Between( ":", "-" ).Int( );
               int stop2 = r2.After( "-" ).Int( );
               if ( IntervalOverlap( start1, stop1, start2, stop2 ) > 0 )
               {    int start = Min( start1, start2 ), stop = Max( stop1, stop2 );
                    r1 = r1.Before( ":" ) + ":" + ToString(start)
                         + "-" + ToString(stop);
                    to_delete[i2] = True;
                    goto combine_regions;    }    }
          EraseIf( regions, to_delete );    }

     vec<int> trace_ids, precorrect_seq;
     ParseIntSet( logc.TRACE_IDS, trace_ids );

     vec<int> trace_pids;
     ParseIntSet( logc.TRACE_PIDS, trace_pids );
     for ( int i = 0; i < trace_pids.isize( ); i++ )
     {    int64_t pid = trace_pids[i];
          trace_ids.push_back( 2*pid, 2*pid + 1 );    }

     ParseIntSet( "{" + heur.PRECORRECT_SEQ + "}", precorrect_seq, false );

     if(   READS.Contains(".bam",String::npos) || READS.Contains(".fastq",String::npos)){
         vec<String> files;
         ParseStringSet( "{" + READS + "}", files );
         bool bHasBAM=false;
         bool bHasFASTQ=false;
         for(const auto& file: files){
             bHasBAM = bHasBAM || file.Contains(".bam",String::npos);
             bHasFASTQ = bHasFASTQ || file.Contains(".fastq",String::npos);
             picard_spec s; s.special = file;
             specs.push_back(s);
         }
         if( bHasFASTQ ) {
             if(specs.size()%2==1 && specs.size()!=1){
                 FatalErr("By design, a FASTQ file list can take either one (interleaved) file, or an even number of files (to interleave the reads from file 2i and 2i+1)");
             }
             if( bHasBAM ){ FatalErr( "Read files contain FASTQ plus BAM/FASTB files" ); }
             if( regions.nonempty() ){
                 if(logc.STATUS_LOGGING)
                 {    cout << Date( ) << ": region doesn't mean anything for FASTQ extraction"
                           << endl;    }
             }
         }

         if(bHasBAM)
         {   if ( regions.nonempty( ) )
             {    if ( SAMPLE == "human" )
                  {    for ( int c = 1; c <= 22; c++ )
                            chrlist.push_back( ToString(c) );
                       chrlist.push_back( "X", "Y", "MT" );    }
                  else
                  {    FatalErr( "For READS = list of bam files, I'm not sure what to "
                            << "do if the SAMPLE is not human." );    }    }
             for ( int i = 0; i < regions.isize( ); i++ )
             {    if ( i > 0 ) region += " ";
                  if ( regions[i].IsInt( ) ) region += chrlist[ regions[i].Int( ) ];
                  else
                  {    Bool bad = False;
                       String chr = regions[i].Before( ":" );
                       String range = regions[i].After( ":" );
                       if ( !range.Contains( "-" ) ) bad = True;
                       if ( SAMPLE == "human" )
                       {    if ( !Member( chrlist, chr ) ) bad = True;    }
                       else if ( !chr.IsInt( ) ) bad = True;
                       String start = range.Before( "-" ), stop = range.After( "-" );
                       if ( !start.IsInt( ) || !stop.IsInt( ) ) bad = True;
                       if (bad) FAIL_MSG( "Failed to parse X argument." );
                       if ( SAMPLE != "human" ) chr = chrlist[ chr.Int( ) ];
                       region += chr + ":" + ToString( start.Int( ) ) + "-"
                        + ToString( stop.Int( ) );    }    }    }
     }
     else
     {    SetupIlluminaData( spec, SAMPLE, DATASET, regions, COVERAGE_CAP,
               spec.HPOOL_ORIG, specs, chrlist, region, SELECT_FRAC, logc );    }
     REPORT_TIME( clock1, "used in initial load" );

     // Get the reads.

     if ( !HAVE_PICARD )
     {    double vclock = WallClockTime( );
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": begin read import, see log files in TMP = " 
                    << TMP << endl << Date( ) << ": loading reads from picard cache" 
                    << endl;    }
          LoadIlluminaReads( SAMPLE, specs, region, TMP, logc, heur.USE_PF_ONLY,
               KEEP_NAMES );
          if ( MastervecFileObjectCount( TMP + "/frag_reads_orig.fastb" ) == 0 )
          {    cout << "\nThere are no reads, nothing to do." << endl;
               Scram( );    }

          // Cut coverage and randomize order.

          if ( SELECT_FRAC < 0 && COVERAGE_CAP < 0 )
          {    if ( SAMPLE == "hpool2" || SAMPLE == "hpool3" )
                    COVERAGE_CAP = spec.HPOOL_COVERAGE_CAP;
               else SELECT_FRAC = 1.0;    }
          if ( COVERAGE_CAP >= 0 )
          {    vecbasevector bases( TMP + "/frag_reads_orig.fastb" );
               double true_size = 0;
               for ( int g = 0; g < (int) (*log_control.G).size( ); g++ )
               {    true_size +=
                         (*log_control.G)[g].size( ) * (*log_control.ploidy)[g];    }
               double nominal_cov = bases.SizeSum( ) / true_size;
               if ( nominal_cov > COVERAGE_CAP )
               {    cout << Date( ) << ": lowering nominal coverage to cap" << endl;
                    SELECT_FRAC = COVERAGE_CAP / nominal_cov;    }
               else SELECT_FRAC = 1.0;    }
          SelectRandom( TMP, SELECT_FRAC, logc, spec );

          // Truncate hpool1 reads to length 250.  Note that this doesn't adjust
          // read locations, so they are off.

          if ( SAMPLE == "hpool1" )
          {    vecbasevector bases( TMP + "/frag_reads_orig.fastb" );
               vecqualvector quals( TMP + "/frag_reads_orig.qualb" );
               for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
               {    bases[id].resize(250);
                    quals[id].resize(250);    }
               bases.WriteAll( TMP + "/frag_reads_orig.fastb" );
               quals.WriteAll( TMP + "/frag_reads_orig.qualb" );    }
          REPORT_TIME( vclock, "used fetching data" );    }

     // If requested, delete reads carrying human satellite DNA.

     if (heur.DELETE_SATELLITE) 
          CleanSatelliteReads( TMP, heur.SATELLITE_TARGETS, heur.MAX_ALPHA_SCORE );

     // If requested, delete Q2 tails.

     if (heur.DELETE_Q2)
     {    vecbasevector bases( TMP + "/frag_reads_orig.fastb" );
          vecqualvector quals( TMP + "/frag_reads_orig.qualb" );
          for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
          {    while( bases[id].size( ) > 0 && quals[id].back( ) == 2 )
               {    bases[id].resize( bases[id].size( ) - 1 );
                    quals[id].resize( quals[id].size( ) - 1 );    }    }
          bases.WriteAll( TMP + "/frag_reads_orig.fastb" );
          quals.WriteAll( TMP + "/frag_reads_orig.qualb" );    }

     // For hpool2 and hpool3, if the entire dataset is to be assembled, clean junk
     // out of the reads.  Also clean out reads from 'finished' Fosmids if requested.

     if ( ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) && X == "" )
     {    CleanFosmidReads(TMP);
          if (spec.REMOVE_FINISHED) RemoveFinished( SAMPLE, TMP );    }

     // Get read locations.

     if ( logc.KEEP_LOCS )
     {    int N = MastervecFileObjectCount( TMP + "/frag_reads_orig.fastb" );
          (log_control.readlocs)->resize(N);
          vec<look_align> aligns;
          LoadLookAligns( TMP + "/frag_reads_orig.qltout", aligns );
          for ( int i = 0; i < aligns.isize( ); i++ )
          {    const look_align& la = aligns[i];
               (*(log_control.readlocs))[ la.query_id ] = ref_loc(
                    la.target_id, la.pos2( ), la.Pos2( ), la.rc1 );    }    }

     // Do the error correction.

     if ( EXIT == "LOAD" ) Done(clock);
     vecbasevector creads;
     CorrectionSuite( TMP, heur, logc, log_control, creads, corrected,
          cid, cpartner, NUM_THREADS, EXIT, clock, useOldLRPMethod );

     // Remove vector reads.

     vec<Bool> to_delete( corrected.size( ), False );
     if ( SAMPLE == "hpool1" && heur.REMOVE_VECTOR )
          FlagVector( corrected, to_delete, logc );

     // Define pairing info.  Note that for now we set all the library ids to 0.

     DefinePairingInfo(
          TMP, creads, to_delete, cid, corrected, cpartner, logc );    }
