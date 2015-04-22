///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ScaffoldInfoProber.  Experimental code to tinker with scaffolding possibilities.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Fastavector.h"
#include "FindGaps.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Superb.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "random/Random.h"
#include <omp.h>

// A tiglet represents part or all of a contig.

class tiglet {

     public:

     tiglet( ) { }
     tiglet( const int tig, const int start, const int stop )
          : tig(tig), start(start), stop(stop) { }

     int tig;
     int start, stop;

     int Len( ) const { return stop - start; }

};

// A tlink is a link between two contigs.  It is always forward on the first 
// contig, reverse on the second.

class tlink {

     public:

     tlink( ) { }
     tlink( const int tig1, const int tig2, const int pos1, const int Pos1, 
          const int pos2, const int Pos2, const int sep, const int dev )
          : tig1(tig1), tig2(tig2), pos1(pos1), Pos1(Pos1), pos2(pos2),
          Pos2(Pos2), sep(sep), dev(dev) { }

     int tig1, tig2;
     int pos1, Pos1;
     int pos2, Pos2;
     int sep, dev;

     friend Bool operator<( const tlink& x, const tlink& y )
     {    if ( x.tig1 < y.tig1 ) return True;
          if ( x.tig1 > y.tig1 ) return False;
          if ( x.tig2 < y.tig2 ) return True;
          if ( x.tig2 > y.tig2 ) return False;
          if ( x.pos1 < y.pos1 ) return True;
          if ( x.pos1 > y.pos1 ) return False;
          if ( x.Pos1 < y.Pos1 ) return True;
          if ( x.Pos1 > y.Pos1 ) return False;
          if ( x.pos2 < y.pos2 ) return True;
          if ( x.pos2 > y.pos2 ) return False;
          if ( x.Pos2 < y.Pos2 ) return True;
          if ( x.Pos2 > y.Pos2 ) return False;
          if ( x.sep < y.sep ) return True;
          if ( x.sep > y.sep ) return False;
          if ( x.dev < y.dev ) return True;
          return False;    }

};

TRIVIALLY_SERIALIZABLE(tlink);

// A linkster represents a link between tiglets.  It is defined with reference to
// the tlinks data structure.

class linkster {

     public:

     linkster( ) { }
     linkster( const int link_start, const int link_stop, const int sep, 
          const int dev ) 
          : link_start(link_start), link_stop(link_stop), sep(sep), dev(dev) { }

     int link_start, link_stop; // [link_start, link_stop) on tlinks
     int sep, dev;

     int NLinks( ) const { return link_stop - link_start; }

     void Print( ostream& out, const vec<tlink>& links, 
          const int m1, const int m2, const int ntigs0 ) const
     {    for ( int k = link_start; k < link_stop; k++ )
          {    const tlink& T = links[k];
               out << ( m1 < ntigs0 ? m1 : m1 - ntigs0 )
                    << ( m1 < ntigs0 ? "fw" : "rc" )
                    << "[" << T.pos1 << "," << T.Pos1 << "] --> " 
                    << ( m2 < ntigs0 ? m2 : m2 - ntigs0 )
                    << ( m2 < ntigs0 ? "fw" : "rc" )
                    << "[" << T.pos2 << "," << T.Pos2 << "]: "
                    << T.sep << " +/- " << T.dev << "\n";    }    }
     void PrintBrief( ostream& out, 
          const int m1, const int m2, const int ntigs0 ) const
     {    out << ( m1 < ntigs0 ? m1 : m1 - ntigs0 )
               << ( m1 < ntigs0 ? "fw" : "rc" ) << " --> "
               << ( m2 < ntigs0 ? m2 : m2 - ntigs0 ) << ( m2 < ntigs0 ? "fw" : "rc" )
               << ": " << sep << " +/- " << dev << "\n";    }

};

TRIVIALLY_SERIALIZABLE(linkster);

// A tigat (pronounced "tig at") is a contig placed in some coordinate system.
// For development purposes, we include genomic coordinates.

class tigat {

     public:

     tigat( ) { }
     tigat( const int t, const int len, const int start, const int dev,
          const int nlinks, const int support, const int gmult, const int gtig, 
          Bool grc, const int gstart, const int gstop ) : t(t), len(len), 
          start(start), dev(dev), nlinks(nlinks), support(support), gmult(gmult), 
          gtig(gtig), grc(grc), gstart(gstart), gstop(gstop) { }

     int t;
     int len;
     int start, dev;
     int nlinks;
     int support;
     int gmult; // number of places on genome
     int gtig;
     Bool grc;
     int gstart, gstop;

     friend Bool operator<( const tigat& T1, const tigat& T2 )
     {    if ( T1.start < T2.start ) return True;
          if ( T1.start > T2.start ) return False;
          if ( T1.dev < T2.dev ) return True;
          if ( T1.dev > T2.dev ) return False;
          if ( T1.t < T2.t ) return True;
          if ( T1.t > T2.t ) return False;
          if ( T1.len < T2.len ) return True;
          if ( T1.len > T2.len ) return False;
          if ( T1.nlinks < T2.nlinks ) return True;
          if ( T1.nlinks > T2.nlinks ) return False;
          if ( T1.support < T2.support ) return True;
          return False;    }

     // Relativize.  This translates genomic coordinates to a scale that mimics
     // the start-stop estimates.

     friend void Relativize( vec<tigat>& X, const int ntigs0 )
     {    int N = X.size( ), next = -1, gtig = 0;
          vec<Bool> changed( N, False ), rc(N);
          for ( int j = 0; j < N; j++ )
          {    if ( X[j].start == 0 ) next = j;
               rc[j] = X[j].grc;
               if ( X[j].gmult == 0 ) changed[j] = True;    }
          ForceAssertGe( next, 0 );
          while(1)
          {    Bool this_rc = rc[next];
               int this_gtig = X[next].gtig; 
               int this_gstart = X[next].gstart, this_gstop = X[next].gstop;
               for ( int k = 0; k < N; k++ )
               {    if ( changed[k] ) continue;
                    if ( rc[k] != this_rc || X[k].gtig != this_gtig ) continue;
                    int glen = X[k].gstop - X[k].gstart;
                    if ( !this_rc ) 
                         X[k].gstart = X[next].start + X[k].gstart - this_gstart;
                    else X[k].gstart = X[next].start + this_gstop - X[k].gstop;
                    X[k].gstop = X[k].gstart + glen;
                    X[k].gtig = gtig, X[k].grc = False;
                    changed[k] = True;    }
               gtig++;
               int closest = -1, closest_score = 1000000000;
               for ( int j = 0; j < N; j++ )
               {    if ( !changed[j] && Abs( X[j].start ) < closest_score )
                    {    closest = j, closest_score = Abs( X[j].start );    }    }
               if ( closest < 0 ) break;
               next = closest;    }    }

     friend void Print( const vec<tigat>& X, const int ntigs0, Bool show_gor )
     {    vec< vec<String> > rows;
          vec<String> row0;
          row0.push_back( "#", "id", "or", "len", "start", "stop", "dev" );
          row0.push_back( "nlinks", "supp", "m", "gtig" );
          if (show_gor) row0.push_back( "gor" );
          row0.push_back( "gstart", "gstop" );
          rows.push_back(row0);
          for ( int i = 0; i < X.isize( ); i++ )
          {    vec<String> row;
               const tigat& x = X[i];
               int tx = ( x.t < (int) ntigs0 ? x.t : x.t - ntigs0 );
               row.push_back( ToString(i), ToString(tx), 
                    ( x.t < (int) ntigs0 ? "fw" : "rc" ), ToString(x.len),
                    ToString(x.start), ToString(x.start+x.len), ToString(x.dev) );
               row.push_back( ToString(x.nlinks), ToString(x.support), 
                    ToString(x.gmult), 
                    ( x.gmult > 0 ? ToString(x.gtig) : String("-") ) );
               if (show_gor) row.push_back(
                    ( x.gmult > 0 ? ( x.grc ? "rc" : "fw" ) : String("--") ) );
               row.push_back(
                    ( x.gmult > 0 ? ToString(x.gstart) : String("-----") ), 
                    ( x.gmult > 0 ? ToString(x.gstop) : String("-----") ) );
               rows.push_back(row);    }
          PrintTabular( cout, rows, 2, "rrrrrrrrrrrrrr" );    }

};

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     EndCommandArguments;

     // Get started.

     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     cout << Date( ) << ": " << run_dir << endl;

     // Set up checkpointing.

     Bool CHECKPOINT = True;
     String checkpoint_dir = sub_dir + "/check";
     Mkdir777(checkpoint_dir);
     String checkpoint_head = checkpoint_dir + "/scaff.";

     // Check for input files.

     String READS = "scaffold_reads", ALIGNS = "scaffold_reads";
     String reads_file = run_dir + "/" + READS + ".fastb";
     String pairs_file = run_dir + "/" + READS + ".pairs"; 
     String contigs_in_file = sub_dir + "/" + "hyper_plus_edited.fasta";
     String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
     String aligns_index_file = sub_dir + "/" + ALIGNS + ".qltoutlet.index";
     vec<String> required;
     required.push_back( reads_file, pairs_file, aligns_file, aligns_index_file );
     required.push_back( contigs_in_file );
     for ( int i = 0; i < required.isize( ); i++ )
     {    if ( !IsRegularFile( required[i] ) )
          {    cout << "Can't find " << required[i] << endl << "Abort." << endl;
               exit(1);    }    }

     // Load contigs.

     cout << Date( ) << ": loading contigs" << endl;
     vec<fastavector> contigs;
     LoadFromFastaFile( contigs_in_file, contigs );
     size_t ntigs0 = contigs.size( );

     // Build alignments to genome reference sequence.

     vec<look_align> galigns;
     if ( !CHECKPOINT || !IsRegularFile( checkpoint_head + "GALIGNS" ) )
     {    size_t nbatches = 4 * omp_get_max_threads( );
          cout << Date( ) << ": aligning contigs to genome reference" << endl;
          #pragma omp parallel for
          for ( size_t i = 0; i < nbatches; i++ )
          {    int start = ( i * ntigs0 ) / nbatches;
               int stop = ( (i+1) * ntigs0 ) / nbatches;
               SystemSucceed( "QueryLookupTable" + ARG(K, 12) + ARG(MM, 12) 
                    + ARG(MC, 0.15) + ARG(SEQS, sub_dir + "/hyper_plus_edited.fastb")
                    + ARG(PARSEABLE, True) + ARG(L, data_dir + "/genome.lookup")
                    + ARG(SEQS_TO_PROCESS, 
                         "\"[" + ToString(start) + "," + ToString(stop) + ")\"")
                    + " > " + checkpoint_head + "GALIGNS." + ToString(i) );    }
          Remove( checkpoint_head + "GALIGNS" );
          for ( size_t i = 0; i < nbatches; i++ )
          {    System( "cat " + checkpoint_head + "GALIGNS." + ToString(i)
                    + " >> " + checkpoint_head + "GALIGNS" );    
               Remove( checkpoint_head + "GALIGNS." + ToString(i) );    }    }
     cout << Date( ) << ": loading genome alignments" << endl;
     vec< vec<int> > galigns_index;
     LoadLookAligns( checkpoint_head + "GALIGNS", galigns, galigns_index, ntigs0 );

     // Load reads and get read lengths.  (Perhaps not needed!)

     cout << Date( ) << ": computing read lengths" << endl;
     vec<uint16_t> read_len;
     String READLEN_file = checkpoint_head + "READLEN";
     if ( CHECKPOINT && IsRegularFile(READLEN_file) )
          BinaryReader::readFile( READLEN_file, &read_len );
     else
     {    vecbasevector reads(reads_file);
          read_len.resize( reads.size( ) );
          for ( size_t i = 0; i < reads.size( ); i++ )
          {    ForceAssertLt( reads[i].size( ), 65536u );
               read_len[i] = reads[i].size( );    }
          if (CHECKPOINT) BinaryWriter::writeFile( READLEN_file, read_len );    }

     // Build links.

     vec<tlink> links;
     if ( CHECKPOINT && IsRegularFile( checkpoint_head + "LINKS" ) )
     {    cout << Date( ) << ": reading links" << endl;
          BinaryReader::readFile( (checkpoint_head + "LINKS").c_str(),
                                  &links );    }
     else
     {
          // Load read pairs.

          cout << Date( ) << ": loading pairs" << endl;
          size_t nreads = MastervecFileObjectCount(reads_file);
          PairsManager pairs( pairs_file );

          // Load alignments of reads to contigs.
          // aligns0_index works as this: say aligns0_index[id] = ii, then:
          //  ii = -2 :   no look_aligns for this read
          //  ii = -1 :   more than one look_align for this read (id)
          //  ii >= 0 :   aligns0[ii] is the only look_align for this read
          
          cout << Date( ) << ": loading aligns" << flush;
          vec<alignlet> aligns0;
          BinaryReader::readFile( aligns_file, &aligns0 );
          vec<int> aligns0_index;
          BinaryReader::readFile( aligns_index_file, &aligns0_index );
          cout << " (" << ToStringAddCommas( aligns0.size( ) ) 
               << " aligns found)" << endl;

          // Add alignments for the rcs of contigs.

          vec<alignlet> aligns;
          vec<int> aligns_index( aligns0_index );
          {    for ( size_t i = 0; i < nreads; i++ )
                    aligns_index[i] = 2 * aligns_index[i];
               aligns.reserve( 2 * aligns0.size( ) );
               for ( size_t i = 0; i < aligns0.size( ); i++ )
               {    alignlet a = aligns0[i];
                    aligns.push_back(a);
                    a.Reverse( );
                    a.SetTargetId( ntigs0 + a.TargetId( ) );
                    aligns.push_back(a);    }
               Destroy(aligns0), Destroy(aligns0_index);    }

          // Build a data structure that defines the links.

          cout << Date( ) << ": defining links" << endl;
          links.reserve( aligns.size( ) / 2 );
          for ( size_t pi = 0; pi < pairs.nPairs( ); pi++ )
          {    int64_t id1 = pairs.ID1(pi), id2 = pairs.ID2(pi);
               int sep = pairs.sep(pi), dev = pairs.sd(pi);
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 ) swap( id1, id2 );
                    if ( aligns_index[id1] < 0 || aligns_index[id2] < 0 ) continue;
                    for ( int j1 = 0; j1 < 2; j1++ )
                    {    const alignlet& a1 = aligns[ aligns_index[id1] + j1 ];
                         for ( int j2 = 0; j2 < 2; j2++ )
                         {    const alignlet& a2 = aligns[ aligns_index[id2] + j2 ];
                              if ( !a1.Fw1( ) || a2.Fw1( ) ) continue;
                              links.push( a1.TargetId( ), a2.TargetId( ), 
                                   a1.pos2( ), a1.Pos2( ), a2.pos2( ), a2.Pos2( ),
                                   sep - ( a1.TargetLength( ) - a1.Pos2( ) ) 
                                   - a2.pos2( ), dev );    }    }    }    }
          cout << Date( ) << ": sorting " << ToStringAddCommas( links.size( ) ) 
               << " links" << endl;
          ParallelSort(links);
          if (CHECKPOINT)
          {    cout << Date( ) << ": writing links" << endl;
               BinaryWriter::writeFile( (checkpoint_head + "LINKS").c_str(),
                                        links );    }    }

     // Combine contigs and their rcs.

     cout << Date( ) << ": combining contigs with their rcs" << endl;
     {    contigs.reserve( 2 * contigs.size( ) );
          for ( size_t i = 0; i < ntigs0; i++ )
          {    contigs.push_back( contigs[i] );
               contigs.back( ).ReverseComplement( );    }    }
     size_t ntigs = contigs.size( );

     // Condense links.

     vec<linkster> linksters;
     vec< pair<int,int> > m1m2;
     if ( CHECKPOINT && IsRegularFile( checkpoint_head + "LINKSTERS" ) )
     {    cout << Date( ) << ": reading linksters" << endl;
          BinaryReader r( (checkpoint_head + "LINKSTERS").c_str() );
          r.read( &linksters );
          r.read( &m1m2 );    }
     else
     {    cout << Date( ) << ": condensing links" << endl;
          size_t nbatches = 10 * omp_get_max_threads( );
          vec< vec<linkster> > linksters_b(nbatches);
          vec< vec< pair<int,int> > > m1m2_b(nbatches);
          #pragma omp parallel for
          for ( size_t bi = 0; bi < nbatches; bi++ )
          {    size_t start = ( bi * links.size( ) ) / nbatches;
               size_t stop = ( (bi+1) * links.size( ) ) / nbatches;
               for ( size_t i = start; i < stop; i++ )
               {    size_t j;
                    for ( j = i + 1; j < links.size( ); j++ )
                    {    if ( links[j].tig1 != links[i].tig1 ) break;
                         if ( links[j].tig2 != links[i].tig2 ) break;    }
                    vec<int> seps, devs;
                    for ( size_t k = i; k < j; k++ )
                    {    seps.push_back( links[k].sep ); 
                         devs.push_back( links[k].dev );    }
                    int sep, dev;
                    GapStats( seps, devs, sep, dev );
                    linksters_b[bi].push( i, j, sep, dev );
                    m1m2_b[bi].push( links[i].tig1, links[i].tig2 );
                    i = j - 1;    }    }
          cout << Date( ) << ": appending" << endl;
          for ( size_t bi = 0; bi < nbatches; bi++ )
          {    linksters.append( linksters_b[bi] );
               m1m2.append( m1m2_b[bi] );    }
          if (CHECKPOINT)
          {    cout << Date( ) << ": writing linksters" << endl;
               BinaryWriter w( (checkpoint_head + "LINKSTERS").c_str() );
               w.write(linksters);
               w.write(m1m2);
               w.close();    }    }
     cout << Date( ) << ": found " << ToStringAddCommas( linksters.size( ) )
           << " linksters" << endl;

     // Define the "scaffold graph" structure.  It is a digraph whose vertices are
     // base ranges on contigs (tiglets), and whose edges are links between them
     // (linksters).  Initially we use entire contigs, but during the scaffolding
     // process, we sometimes split them.

     cout << Date( ) << ": defining scaffold graph" << endl;
     vec<tiglet> tiglets;
     digraphE<linkster> G;
     {    tiglets.reserve(ntigs);
          for ( size_t t = 0; t < ntigs; t++ )
               tiglets.push( t, 0, contigs[t].size( ) );
          size_t N = tiglets.size( );
          vec< vec<int> > &from = G.FromMutable( ), &to = G.ToMutable( );
          vec< vec<int> >& from_edge_obj = G.FromEdgeObjMutable( );
          vec< vec<int> >& to_edge_obj = G.ToEdgeObjMutable( );
          from.resize(N), to.resize(N);
          from_edge_obj.resize(N), to_edge_obj.resize(N);
          cout << Date( ) << ": filling in from, etc." << endl;
          for ( size_t i = 0; i < linksters.size( ); i++ )
          {    int m1 = m1m2[i].first, m2 = m1m2[i].second;
               from[m1].push_back(m2), to[m2].push_back(m1);
               from_edge_obj[m1].push_back(i), to_edge_obj[m2].push_back(i);    }
          cout << Date( ) << ": sorting" << endl;
          #pragma omp parallel for
          for ( size_t i = 0; i < N; i++ )
          {    SortSync( from[i], from_edge_obj[i] );
               SortSync( to[i], to_edge_obj[i] );    }
          G.EdgesMutable( ) = linksters;
          Destroy(linksters), Destroy(m1m2);    }
     cout << Date( ) << ": done initializing" << endl;

     // Start the process of building neighborhoods.

     for ( int attempt = 1; attempt <= 5; attempt++ )
     {    int t = 0;
          while(1)
          {    t = randomx( ) % G.N( );
               if ( tiglets[t].Len( ) >= 5000 ) break;    }
          int t0 = ( t < (int) ntigs0 ? t : t - (int) ntigs0 );
          cout << "\nATTEMPT " << attempt << ", NEIGHBORHOOD OF CONTIG "
               << t0 << " " << ( t < (int) ntigs0 ? "fw" : "rc" ) << endl;
          vec<tigat> X;
          int gmult = galigns_index[t0].size( ), gtig = 0, gstart = 0, gstop = 0;
          Bool grc = False;
          if ( gmult > 0 )
          {    const look_align& la = galigns[ galigns_index[t0][0] ];
               gtig = la.target_id;
               grc = la.Rc1( );
               gstart = la.pos2( ), gstop = la.Pos2( );    }
          X.push( t, tiglets[t].Len( ), 0, 0, 0, tiglets[t].Len( ), gmult, gtig,
               grc ^ ( t >= (int) ntigs0 ), gstart, gstop );

          for ( int pass = 1; pass <= 2; pass++ )
          {    
               PRINT( X.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

               // Find contigs that are linked to.

               vec<int> have, to, tox;
               for ( int i = 0; i < X.isize( ); i++ )
               {    have.push_back( X[i].t );
                    if ( X[i].t < (int) ntigs0 ) have.push_back( X[i].t + ntigs0 );
                    else have.push_back( X[i].t - ntigs0 );    }
               for ( int i = 0; i < X.isize( ); i++ )
               {    to.append( G.From( X[i].t ) ), to.append( G.To( X[i].t ) );    }
               UniqueSort(have), UniqueSort(to);
               for ( int i = 0; i < to.isize( ); i++ )
                    if ( !BinMember( have, to[i] ) ) tox.push_back( to[i] );
               to = tox;

               // Find and condense links to each.

               int nX = X.size( );
               for ( int l = 0; l < to.isize( ); l++ )
               {    int tx = to[l], nlinks = 0;
                    vec<int> supports, starts, devs;
                    for ( int i = 0; i < nX; i++ )
                    {    int t2 = X[i].t;
                         for ( int j = 0; j < G.From(t2).isize( ); j++ )
                         {    if ( G.From(t2)[j] != tx ) continue;
                              const linkster& L = G.EdgeObjectByIndexFrom( t2, j );
                              nlinks += L.NLinks( );
                              vec<int> pos1, pos2;
                              // L.Print( cout, links, t2, tx, ntigs0 ); // XXXXXXXX
                              // L.PrintBrief( cout, t2, tx, ntigs0 ); // XXXXXXXXXX
                              for ( int k = L.link_start; k < L.link_stop; k++ )
                              {    const tlink& T = links[k];
                                   pos1.push_back(T.pos1); 
                                   pos2.push_back(T.pos2);    }
                              Sort(pos1), Sort(pos2);
                              supports.push_back( std::min( {pos1.back( ) - pos1.front( ),
                                   pos2.back( ) - pos2.front( ), X[i].support} ) );
                              starts.push_back( 
                                   X[i].start + tiglets[t2].Len( ) + L.sep );
                              devs.push_back( int( round( sqrt( double( X[i].dev 
                                   * X[i].dev + L.dev * L.dev ) ) ) ) );    }
                         for ( int j = 0; j < G.To(t2).isize( ); j++ )
                         {    if ( G.To(t2)[j] != tx ) continue;
                              const linkster& L = G.EdgeObjectByIndexTo( t2, j );
                              nlinks += L.NLinks( );
                              vec<int> pos1, pos2;
                              // L.Print( cout, links, tx, t2, ntigs0 ); // XXXXXXXX
                              // L.PrintBrief( cout, tx, t2, ntigs0 ); // XXXXXXXXXX
                              for ( int k = L.link_start; k < L.link_stop; k++ )
                              {    const tlink& T = links[k];
                                   pos1.push_back(T.pos1); 
                                   pos2.push_back(T.pos2);    }
                              Sort(pos1), Sort(pos2);
                              supports.push_back( std::min( {pos1.back( ) - pos1.front( ),
                                   pos2.back( ) - pos2.front( ), X[i].support} ) );
                              starts.push_back( 
                                   X[i].start - tiglets[tx].Len( ) - L.sep );
                              devs.push_back( int( round( sqrt( double( X[i].dev 
                                   * X[i].dev + L.dev * L.dev ) ) ) ) );    }    }
                    int support = Max(supports);
                    const int min_support = 3;
                    const int min_links_unsupported = 3;
                    if ( support < min_support && nlinks < min_links_unsupported )
                         continue;
                    if ( supports.size( ) > 1 ) support = Max( 1000, support );
                    if ( support < min_support ) continue;
                    const int radius = 40000;
                    int start, dev;
                    GapStats( starts, devs, start, dev );
                    int tx0 = ( tx < (int) ntigs0 ? tx : tx - (int) ntigs0 );
                    int gmult = galigns_index[tx0].size( ); 
                    int gtig = 0, gstart = 0, gstop = 0;
                    if ( gmult > 0 )
                    {    const look_align& la = galigns[ galigns_index[tx0][0] ];
                         gtig = la.target_id;
                         grc = la.Rc1( );
                         gstart = la.pos2( ), gstop = la.Pos2( );    }
                    int len = tiglets[tx].Len( );
                    if ( IntervalOverlap( start, start+len, -radius, radius ) == 0 )
                         continue;
                    X.push( tx, len, start, dev, nlinks, support, gmult, gtig, 
                         grc ^ ( tx >= (int) ntigs0 ), gstart, gstop );    }    }
          Sort(X);
          Relativize( X, ntigs0 );
          Print( X, ntigs0, False );    }

     // Finish.

     cout << "\n" << Date( ) << ": done, time used = " << TimeSince(clock) 
          << endl;    }
