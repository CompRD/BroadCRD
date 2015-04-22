/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <sstream>

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "bias/DinukeBias.h"
#include "bias/StartBias.h"
#include "bias/UniformBias.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"

String StartBias( const vecbasevector& reads, const vecbasevector& ref,
     const vecbasevector& rref, const vec<look_align>& aligns, const Bool brief,
     String marked_reference_file, const Bool summary_form, const Bool silent )
{
     // Compute read start points.

     ForceAssertGt( ref.size( ), 0u );
     vec< vec<int> > starts( 2 * ref.size( ) );
     for ( int i = 0; i < starts.isize( ); i++ )
          starts[i].resize( ref[i/2].size( ), 0 );
     for ( int j = 0; j < aligns.isize( ); j++ )
     {    const look_align& la = aligns[j];
          int id = la.query_id, id2 = la.target_id;
          align a = la.a;
          const basevector& rd1 = reads[la.query_id];
          const basevector& rd2 = ( !la.rc1 ? ref[id2] : rref[id2] );
          if (la.rc1) a.ReverseThis( rd1.size( ), rd2.size( ) );
          ForceAssertGe( a.pos2( ), 0 );
          ForceAssertLe( static_cast<unsigned>(a.Pos2()), rd2.size( ) );
          ++starts[ (2*id2) + (la.rc1 ? 1 : 0) ][ a.pos2( ) ];    }

     // Compute start point bias.

     String bias = StartBias( reads, ref, rref, starts, brief, 
          marked_reference_file, summary_form, silent );
     if (brief) return bias;

     // Generate dinucleotide stats.

     cout << "\n";
     DinukeBias( reads, ref, rref, aligns );

     // Generate stats for GC content of windows starting at the read start point.
     // For the moment the window size is hardwired.

     unsigned int window = 50;
     vec<int> GC( window+1, 0 ), GC_ref( window+1, 0 );
     for ( int i = 0; i < starts.isize( ); i++ )
     {    const basevector& R = ( i % 2 == 0 ? ref[i/2] : rref[i/2] );
          for ( unsigned int p = 0; p + window <= R.size( ); p++ )
          {    int gc = 0;
               for ( unsigned int k = 0; k < window; k++ )
                    if ( R[p+k] == 1 || R[p+k] == 2 ) ++gc;
               GC[gc] += starts[i][p];
               GC_ref[gc] += 1;    }    }
     cout << "\nGC content statistics:\n\n";
     cout << "\nDistribution of number of GCs in " << window << "-based window "
          << "starting at read start point:\n\n";
     cout << "Field 1. GC content of 50-base window\n";
     cout << "Field 2. Percent of read starts whose window has this GC content\n";
     cout << "Field 3. Bias multiplier for given GC content\n";
     cout << "Field 4. Number of read starts for this GC content\n";
     cout << "Field 5. Number of windows on reference for this GC content\n\n";
     longlong total = BigSum(GC);
     double mult = double(BigSum(GC_ref)) / double(BigSum(GC));
     vec< vec<String> > rows;
     for ( unsigned int i = 0; i <= window; i++ )
     {    vec<String> row;
          ostringstream out1, out2;
          out1 << setiosflags(ios::fixed) << setprecision(4)
               << 100.0 * double( GC[i] ) / double(total) << "%";
          if ( GC_ref[i] == 0 ) out2 << "N/A";
          else 
          {    out2 << setiosflags(ios::fixed) << setprecision(4) 
                    << mult * double( GC[i] ) / double( GC_ref[i] );    }
          row.push_back( ToString(2*i), out1.str( ), out2.str( ) );
          row.push_back( ToString( GC[i] ), ToString( GC_ref[i] ) );
          rows.push_back(row);    }
     PrintTabular( cout, rows, 2, "rrrrr" );

     // Return answer.

     return bias;    }

String StartBias( const vecbasevector& reads, const vecbasevector& ref,
     const vecbasevector& rref, const vec<placement_mark>& places, const Bool brief,
     String marked_reference_file, const Bool summary_form, const Bool silent )
{
     // Compute read start points.

     ForceAssertGt( ref.size( ), 0u );
     vec< vec<int> > starts( 2 * ref.size( ) );
     for ( int i = 0; i < starts.isize( ); i++ )
          starts[i].resize( ref[i/2].size( ), 0 );
     for ( int i = 0; i < places.isize( ); i++ )
     {    const placement_mark& p = places[i];
          unsigned int pos = p.Pos( );
          if ( p.Fw1( ) ) pos = ref[ p.Tig( ) ].size( ) - pos - 1;
          ++starts[ 2*p.Tig( ) + ( p.Fw1( ) ? 0 : 1 ) ][pos];    }

     return StartBias( reads, ref, rref, starts, brief, marked_reference_file, 
          summary_form, silent );    }

String StartBiasRelRaw( const vecbasevector& ref, const vecbasevector& rref, 
     const vec<placement_mark>& places1, const vec<placement_mark>& places2 )
{    int nref = ref.size( );
     ForceAssertGt( nref, 0 );
     vec< vec<int> > starts1( 2 * nref ), starts2( 2 * nref );
     for ( int i = 0; i < 2 * nref; i++ )
     {    starts1[i].resize( ref[i/2].size( ), 0 );
          starts2[i].resize( ref[i/2].size( ), 0 );    }
     for ( int pass = 1; pass <= 2; pass++ )
     {    const vec<placement_mark>& places = ( pass == 1 ? places1 : places2 );
          vec< vec<int> >& starts = ( pass == 1 ? starts1 : starts2 );
          for ( int i = 0; i < places.isize( ); i++ )
          {    const placement_mark& p = places[i];
               unsigned int pos = p.Pos( );
               if ( p.Fw1( ) ) pos = ref[ p.Tig( ) ].size( ) - pos - 1;
               ++starts[ 2*p.Tig( ) + ( p.Fw1( ) ? 0 : 1 ) ][pos];    }    }
     vec<int> NN;
     for ( size_t i = 0; i < ref.size( ); i++ )
          NN.push_back( ref[i].size( ), ref[i].size( ) );
     ostringstream out;
     RightPrecisionOut( out, HowBiasedRelRaw( starts1, starts2, NN ), 1 );
     return out.str( );    }

String StartBias( const vecbasevector& reads, const vecbasevector& ref,
     const vecbasevector& rref, const vec< vec<int> >& starts, const Bool brief,
     String marked_reference_file, const Bool summary_form, const Bool silent )
{
     // Set up.

     int nreads = reads.size( );
     ForceAssertGt( ref.size( ), 0u );

     // Show coverage by reference contig.

     if ( !brief )
     {    cout << "Number of starts per base by reference contig:\n";
          for ( size_t i = 0; i < ref.size( ); i++ )
          {    cout << i << ": " << setprecision(3) 
                    << double( Sum( starts[ 2*i ] ) + Sum( starts[ 2*i + 1 ] ) )
                         / double( ref[i].size( ) ) << "\n";    }
          cout << "\n";    }

     // Create marked reference file.

     if ( marked_reference_file != "" )
     {    Ofstream( mout, marked_reference_file );
          for ( int pass = 0; pass <  2; pass++ )
          {    for ( size_t i = 0; i < ref.size( ); i++ )
               {    const basevector& R = ( pass == 0 ? ref[i] : rref[i] );
                    mout << ">" << i << ( pass == 0 ? "fw" : "rc" ) << "\n";
                    int count = 0;
                    for ( unsigned int j = 0; j < R.size( ); j++ )
                    {    int m = starts[ 2*i + pass ][j];
                         for ( int k = 0; k < m; k++ )
                         {    if ( count > 0 && count % 80 == 0 ) mout << "\n";
                              count++;
                              mout << '^';    }
                         if ( count > 0 && count % 80 == 0 ) mout << "\n";
                         count++;
                         mout << as_base( R[j] );    }
                    mout << "\n";    }    }    }

     // Compute duplicate start points.

     int nbases = 0, ndups = 0;
     for ( int j = 0; j < starts.isize( ); j++ )
     {    nbases += starts[j].size( );
          for ( int u = 0; u < starts[j].isize( ); u++ )
               if ( starts[j][u] > 1 ) ndups += starts[j][u] - 1;    }
     if ( !brief ) 
     {    PRINT3( nbases, nreads, ndups );
          cout << "\n";    }

     // Compute start point bias.

     ostringstream out;
     vec<int> NN;
     for ( size_t i = 0; i < ref.size( ); i++ )
          NN.push_back( ref[i].size( ), ref[i].size( ) );
     out << ( summary_form ? "{{" : "r = 1: " );
     RightPrecisionOut( out, HowBiased( starts, NN ), 1 );
     out << ( summary_form ? "," : " +/- " );
     RightPrecisionOut( out, HowBiasedDev( starts, NN ), 1 );
     out << ( summary_form ? "}," : "\n" );
     vec< vec<int> > startsx;
     int r = 1;
     for ( int j = 0; j < 3; j++ )
     {    r *= 10;
          Smooth( starts, r, startsx );
          if ( summary_form && j > 0 ) out << ",";
          out << ( summary_form ? "{" : "r = " + ToString(r) + ": " );
          RightPrecisionOut( out, HowBiased( startsx, NN ), 1 );
          out << ( summary_form ? "," : " +/- " );
          RightPrecisionOut( out, HowBiasedDev( startsx, NN ), 1 );
          out << ( summary_form ? "}" : "\n" );    }
     if (summary_form) out << "}";
     if ( !silent ) cout << out.str( );
     if ( summary_form && !silent ) cout << "\n";
     if ( !silent ) flush(cout);
     return out.str( );    }
