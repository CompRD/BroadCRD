///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Find tumor-only simple bubbles.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/cancer/SomaticNormal.h"
#include "polymorphism/Edit.h"

class var {

     public:

     var( ) { }
     var( const String& chr, const int start, const String& normal, 
          const String& tumor )
     {    chr_ = chr;
          start_ = start;
          normal_ = normal;
          tumor_ = tumor;    }

     String chr_;
     int start_;
     String normal_;
     String tumor_;

     friend Bool operator<( const var& v1, const var& v2 )
     {    if ( v1.chr_ < v2.chr_ ) return True;
          if ( v1.chr_ > v2.chr_ ) return False;
          if ( v1.start_ < v2.start_ ) return True;
          if ( v1.start_ > v2.start_ ) return False;
          if ( v1.normal_ < v2.normal_ ) return True;
          if ( v1.normal_ > v2.normal_ ) return False;
          if ( v1.tumor_ < v2.tumor_ ) return True;
          if ( v1.tumor_ > v2.tumor_ ) return False;
          return False;    }

};

int main( )
{    RunTime( );

     // Define hardcoded paths.

     // String work_dir = "/wga/scr4/jaffe/GapToy/50935.HCC1954+BL";
     String work_dir = "/wga/scr4/jaffe/GapToy/50921.HCC1143+BL";
     String DIR = work_dir + "/a.final";
     String fin_dir = DIR;
     Bool ALTREF = True;
     String suffix = ( ALTREF ? "_alt" : "" );

     // Load genome.

     vecbasevector genome( DIR + "/../genome.fastb" + suffix );
     int ng = genome.size( );
     vecbasevector genome2(genome);
     {    vecbasevector genome_rc(genome);
          for ( int g = 0; g < (int) genome.size( ); g++ )
               genome_rc[g].ReverseComplement( );
          genome2.Append(genome_rc);    }
     String gnf = DIR + "/../genome.names" + suffix ;
     fast_ifstream in(gnf);
     String line;
     vec<String> genome_names;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          genome_names.push_back(line);    }

     // Parse known variants.

     vec<var> vars;
     fast_ifstream sin( "/cga/meyerson/home/marcin/Projects/Jumps/files/"
          "HCC1143-wgs-pcrfree.snp.maf" );
     vec<String> x;
     while(1)
     {    getline( sin, line );
          if ( sin.fail( ) ) break;
          if ( !line.Contains( "KEEP" ) ) continue;
          Tokenize( line, '\t', x );
          String chr = x[1];
          int start = x[2].Int( ), stop = x[3].Int( );
          String ar = x[4], at1 = x[5], at2 = x[6];
          ForceAssert( at1 == ar );
          ForceAssert( at2 != ar );
          String normal = ar, tumor = at2;
          vars.push( chr, start, normal, tumor );
          int g = Position( genome_names, chr );
          if ( g < 0 )
          {    PRINT(chr);
               ForceAssert( 0 == 1 );    }
          ForceAssertEq( as_base( genome[g][start-1] ), normal[0] );    }
     /*
     fast_ifstream iin( "/cga/meyerson/home/marcin/Projects/Jumps/files/"
          "HCC1143-wgs-pcrfree.indel.maf" );
     */
     fast_ifstream iin( "/wga/scr4/cancer/"
          "HCC1143-wgs_pcrfree_coclean.strelka.pass.somatic.indels.maf" );
     while(1)
     {    getline( iin, line );
          if ( iin.fail( ) ) break;
          if ( !line.Contains( "KEEP" ) ) continue;
          Tokenize( line, '\t', x );
          String chr = x[1];
          int start = x[2].Int( ), stop = x[3].Int( );
          String ar = x[4], at1 = x[5], at2 = x[6];
          String normal, tumor;
          if ( ar != at2 && ar == at1 )
          {    normal = ar, tumor = at2;    }
          else if ( ar != at1 && ar == at2 )
          {    normal = ar, tumor = at1;    }
          else if ( ar != at1 && at1 == at2 )
          {    normal = ar, tumor = at1;    }
          else ForceAssert( 0 == 1 );
          vars.push( chr, start, normal, tumor );
          int g = Position( genome_names, chr );
          if ( g < 0 )
          {    PRINT(chr);
               ForceAssert( 0 == 1 );    }
          if ( normal != "-" )
          {    for ( int j = 0; j < normal.isize( ); j++ )
               {    if ( as_base( genome[g][ start-1+j ] ) != normal[j] )
                    {    PRINT5( chr, start, stop, normal, tumor );    
                         cout << "confused" << endl;
                         Scram(0);    }    }    }    }
     Sort(vars);
     Ofstream( out, "varlist" );
     for ( int i = 0; i < vars.isize( ); i++ )
     {    String chr = vars[i].chr_;
          int start = vars[i].start_;
          String tumor = vars[i].tumor_, normal = vars[i].normal_;
          PRINT4_TO( out, chr, start, tumor, normal );    }

     // Load assembly.

     HyperBasevectorX hb;
     BinaryReader::readFile( DIR + "/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( DIR + "/a.inv", &inv );
     vec<vec<int>> count;
     BinaryReader::readFile( DIR + "/a.countsb", &count );
     vec< vec< pair<int,int> > > hits;
     BinaryReader::readFile( fin_dir + "/a.aligns" + suffix, &hits );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( DIR + "/a.lines", &lines );
     vec<int> tol;
     GetTol( hb, lines, tol );

     // Find hits in bad regions.  Use hg19.

     vec<Bool> bad( hb.E( ), False );
     FindBadEdges( hb, inv, lines, hits, bad, ALTREF );

     // Load new alignments.

     vec< vec< pair<int,align> > > alignsx;
     vec< pair<int,int> > bad_pairs;
     BinaryReader::readFile( DIR + "/a.alignsx" + suffix, &alignsx );
     BinaryReader::readFile( DIR + "/a.bad_pairs" + suffix, &bad_pairs );

     // Go through bubbles.

     int ecount = 0;
     for ( int v = 0; v < hb.N( ); v++ )
     {    if ( hb.To(v).size( ) != 1 || hb.From(v).size( ) != 2 ) continue;
          if ( hb.From(v)[0] != hb.From(v)[1] ) continue;
          int w = hb.From(v)[0];
          if ( hb.From(w).size( ) != 1 ) continue;
          int x = hb.To(v)[0], y = hb.From(w)[0];
          vec<int> s = { x, v, w, y };
          UniqueSort(s);
          if ( s.size( ) != 4 ) continue;
          int e1 = hb.IFrom( v, 0 ), e2 = hb.IFrom( v, 1 );
          if ( bad[e1] || bad[e2] ) continue;
          if ( count[1][e1] > 0 ) swap( e1, e2 );
          if ( count[1][e1] > 0 ) continue;
          if ( count[0][e1] < 10 ) continue;

          Bool bad1 = False, bad2 = False;
          for ( int j = 0; j < bad_pairs.isize( ); j++ )
          {    if ( bad_pairs[j].first == e1 || bad_pairs[j].second == e1 )
                    bad1 = True;
               if ( bad_pairs[j].first == e2 || bad_pairs[j].second == e2 )
                    bad2 = True;    }

          // Test for stuff not to show.

          if ( alignsx[e1].solo( ) && alignsx[e2].solo( )
               && alignsx[e1][0].first >= ng && alignsx[e2][0].first >= ng )
          {    continue;    }

          cout << "\n<" << ++ecount << "> ";
          PRINT2( e1, e2 );    
          ostringstream out;
          out << "e1 has " << alignsx[e1].size( ) << " alignments" << endl;
          out << "e2 has " << alignsx[e2].size( ) << " alignments" << endl;
          for ( int m = 1; m <= 2; m++ )
          {    int e = ( m == 1 ? e1 : e2 );
               for ( int j = 0; j < alignsx[e].isize( ); j++ )
               {    int g = alignsx[e][j].first;
                    Bool fw = True;
                    if ( g >= (int) genome.size( ) )
                    {    fw = False;
                         g -= genome.size( );    }
                    out << "[" << j+1 << "] e" << ToString(m) << " aligned at " 
                         << ( fw ? "+" : "-" ) << g << "." 
                         << alignsx[e][j].second.Extent2( ) << endl;    }    }
          if (bad1) out << "don't trust alignment of e1" << endl;
          if (bad2) out << "don't trust alignment of e2" << endl;    

          if ( alignsx[e1].solo( ) && alignsx[e2].solo( ) && !bad1 && !bad2 )
          {    vec< pair<int,edit0> > edits1, edits2;
               int g1 = alignsx[e1][0].first, g2 = alignsx[e2][0].first;
               edits1 = AlignToEdits( 
                    alignsx[e1][0].second, hb.EdgeObject(e1), genome2[g1] );
               edits2 = AlignToEdits( 
                    alignsx[e2][0].second, hb.EdgeObject(e2), genome2[g2] );
               if ( g1 == g2 && Subset( edits2, edits1 ) )
               {    cout << "somatic mutation:\n";
                    vec< pair<int,edit0> > add;
                    for ( int i = 0; i < edits1.isize( ); i++ )
                    {    if ( !Member( edits2, edits1[i] ) )
                              add.push_back( edits1[i] );    }
                    for ( int i = 0; i < add.isize( ); i++ )
                    {    int g = g1, gpos = add[i].first;
                         cout << i+1 << ". "
                              << genome_names[g] << "." << gpos + 1 << ": "
                              << add[i].second;
                         String from, to;
                         int start;

                         if ( add[i].second.etype == SUBSTITUTION )
                         {    from = String( as_base( genome[g][gpos] ) );
                              to = add[i].second.seq;
                              start = gpos;    }

                         else if ( add[i].second.etype == DELETION )
                         {    int n = add[i].second.n;
                              /* vcf version
                              for ( int j = -1; j < n; j++ )
                                   from += as_base( genome[g][gpos+j] );
                              to += as_base( genome[g][gpos-1] );
                              start = gpos - 1;    
                              */
                              for ( int j = 0; j < n; j++ )
                                   from += as_base( genome[g][gpos+j] );
                              start = gpos;    }

                         else // if ( add[i].second.etype == INSERTION )
                         {    /* vcf version
                              from += as_base( genome[g][gpos-1] );
                              to += as_base( genome[g][gpos-1] );
                              to += add[i].second.seq;
                              start = gpos - 1;    
                              */
                              to = add[i].second.seq;
                              start = gpos;    }

                         if ( from == "" ) from = "-";
                         if ( to == "" ) to = "-";

                         var v( genome_names[g], start + 1, from, to );
                         Bool known = BinMember( vars, v );
                         cout << " (" << ( known ? "known" : "UNKNOWN" ) << ")";
                         cout << endl;    }    }
               else
               {    cout << out.str( );
                    cout << "edits1:\n";
                    for ( int i = 0; i < edits1.isize( ); i++ )
                    {    cout << i+1 << ". " << g1 << "." << edits1[i].first << ": "
                              << edits1[i].second << endl;    }
                    cout << "edits2:\n";
                    for ( int i = 0; i < edits2.isize( ); i++ )
                    {    cout << i+1 << ". " << g2 << "." << edits2[i].first << ": "
                              << edits2[i].second << endl;    }    }    }
          else cout << out.str( );    }    }
