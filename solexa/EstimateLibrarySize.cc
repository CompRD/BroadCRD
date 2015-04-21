/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// EstimateLibrarySize.  Estimate the number of distinct molecules in an
// Illumina library.  Also compute some other interesting stuff.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "solexa/SolexaMetrics.h"
#include "math/Functions.h"

// MakeDepend: dependency ComputeUniqueMolecules

template<int k> class bytes {

     public:

     unsigned char x[k];

     friend Bool operator==( const bytes& b1, const bytes& b2 )
     {    for ( int i = 0; i < k; i++ )
               if ( b1.x[i] != b2.x[i] ) return False;
          return True;    }


     friend Bool operator<( const bytes& b1, const bytes& b2 )
     {    for ( int i = 0; i < k; i++ )
          {    if ( b1.x[i] < b2.x[i] ) return True;
               if ( b1.x[i] > b2.x[i] ) return False;    }
          return False;    }

};

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault(LANES, "");
     CommandArgument_String_OrDefault(LIBRARY, "");
     CommandArgument_String_OrDefault(PIPELINE, Getenv("SOLEXA_PIPELINE_DIR"));
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Double_OrDefault(MAX_ERRS_HEAD, 0.5);
     CommandArgument_Double_OrDefault(MAX_ERRS_ALL, -1.0);
     EndCommandArguments;

     // Get list of lanes.

     if ( !( LANES == "" ^ LIBRARY == "" ) )
     {    cout << "Exactly one of LANES or LIBRARY must be specified." << endl;
          exit(1);    }
     vec<String> lanesx, lanes0, dates0, lanes, dates, libraries;
     vec<String> all = AllFiles(PIPELINE);
     vec<String> lanes_all;
     if ( LANES != "" ) ParseStringSet( LANES, lanesx );
     for ( int i = 0; i < all.isize( ); i++ )
     {    if ( !all[i].Contains( "_" ) ) continue;
          if ( all[i].size( ) > 12 ) continue;
          String date = all[i].Before( "_" ), FC = all[i].After( "_" );
          String dir = PIPELINE + "/" + all[i];
          for ( int L = 1; L <= 8; L++ )
          {    String FCL = FC + "." + ToString(L);
               if ( LANES != "" && !BinMember( lanesx, FCL ) ) continue;
               String HEAD = dir + "/" + FCL;
               if ( !IsRegularFile( HEAD + ".metrics" ) ) continue;
               solexa_metric_db db( HEAD + ".metrics" );
               if ( !db.Defined( "library_name" ) ) continue;
               if ( LIBRARY != "" && db.Value( "library_name" ) != LIBRARY )
                    continue;
               lanes_all.push_back(FCL);
               libraries.push_back( db.Value( "library_name" ) );
               lanes0.push_back(FCL), dates0.push_back(date);    }    }
     UniqueSort(lanes_all);
     if ( lanes0.empty( ) )
     {    cout << "No lanes found!" << endl;
          exit(1);    }
     UniqueSort(libraries);
     if ( !libraries.solo( ) )
     {    cout << libraries.size( ) << " libraries, should have exactly one"
               << endl;
          exit(1);    }
     if (VERBOSE) cout << lanes0.size( )/2 << " lanes requested\n" << endl;
     SortSync(lanes0, dates0);
     for ( int i = 0; i < lanes0.isize( ); i++ )
     {    int j = lanes0.NextDiff(i);
          if ( j - i != 2 ) continue;
          String head1 = PIPELINE + "/" + dates0[i] + "_"
               + lanes0[i].Before( "." ) + "/" + lanes0[i];
          String head2 = PIPELINE + "/" + dates0[i+1] + "_"
               + lanes0[i+1].Before( "." ) + "/" + lanes0[i+1];
          if ( IsRegularFile( head1 + ".fastb" )
               && IsRegularFile( head1 + ".qualb" )
               && IsRegularFile( head2 + ".fastb" )
               && IsRegularFile( head2 + ".qualb" )
               && MastervecFileObjectCount( head1 + ".fastb" )
                    == MastervecFileObjectCount( head2 + ".fastb" ) )
          {    lanes.push_back( lanes0[i], lanes0[i] );
               dates.push_back( dates0[i], dates0[i+1] );
               if (VERBOSE) cout << "using " << lanes0[i] << endl;    }
          i = j - 1;    }
     int nlanes = lanes.size( ) / 2;
     if ( nlanes == 0 )
     {    cout << "No lanes found." << endl;
          exit(0);    }

     // Map quality scores to probabilities.

     vec<double> prob(256);
     for ( int i = 0; i < 256; i++ )
          prob[i] = pow( 10.0, double(-i)/10.0 );

     // Define key size.

     const int nbases = 20;
     ForceAssert( nbases % 2 == 0 );
     const int nbytes = nbases/4;

     // Go through the lanes.

     vec< vec< bytes<nbytes> > > X(nlanes);
     longlong total_pairs = 0;
     longlong total_PF_pairs = 0;
     for ( int l = 0; l < lanes.isize( ); l += 2 )
     {
          // Load bases and quals.

          String head1 = PIPELINE + "/" + dates[l] + "_"
               + lanes[l].Before( "." ) + "/" + lanes[l];
          String head2 = PIPELINE + "/" + dates[l+1] + "_"
               + lanes[l+1].Before( "." ) + "/" + lanes[l+1];
          vecbasevector bases1( head1 + ".fastb" ), bases2( head2 + ".fastb" );
          vecqualvector quals1( head1 + ".qualb" ), quals2( head2 + ".qualb" );
          ForceAssertEq( bases1.size( ), quals1.size( ) );
          ForceAssertEq( bases2.size( ), quals2.size( ) );
          int N = bases1.size( );
          total_pairs += N;

          // Figure out which reads are PF.

          vec<int> pf;
          fast_ifstream in( head1 + ".filter.vector" );
          String line;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               pf.push_back( line.Int( ) );    }
          total_PF_pairs += Sum(pf);

          // Go through the reads.

          for ( int id = 0; id < N; id++ )
          {
               // Test for quality.

               double errs_head = 0.0, errs_all = 0.0;;
               for ( int pass = 1; pass <= 2; pass++ )
               {    const qualvector& q = ( pass == 1 ? quals1[id] : quals2[id] );
                    for ( int j = 0; j < nbases; j++ )
                    {    errs_head += prob[ q[j] ];
                         if ( MAX_ERRS_ALL > 0 )
                              errs_all += prob[ q[j] ];    }    }
               if ( errs_head > MAX_ERRS_HEAD ) continue;
               if ( MAX_ERRS_ALL > 0 )
               {    for ( int pass = 1; pass <= 2; pass++ )
                    {    const qualvector& q
                              = ( pass == 1 ? quals1[id] : quals2[id] );
                         for ( qvec::size_type j = nbases; j < q.size( ); j++ )
                              errs_all += prob[ q[j] ];    }
                    if ( errs_all > MAX_ERRS_ALL ) continue;    }

               // Test for noise.

               static basevector b1, b2;
               b1.SetToSubOf( bases1[id], 0, nbases );
               b2.SetToSubOf( bases2[id], 0, nbases );
               static basevector A( "AAAAAAAAAAAAAAAAAAAA" );
               if ( b1 == A && b2 == A ) continue;

               // Create key.

               static basevector key;
               if ( b1 < b2 ) key = Cat(b1, b2);
               else key = Cat(b2, b1);
               static bytes<nbytes> K;
               for ( int i = 0; i < nbytes; i++ )
                    K.x[i] = key.extractKmer(4*i,4);
               X[l/2].push_back(K);    }     }

     // Handle multi-lane case.

     int duplicates = 0;
     longlong yield = 0;
     if ( nlanes > 1 )
     {
          // Find cross set.

          vec<int> id( nlanes, vec<int>::IDENTITY ), sizes(nlanes);
          for ( int i = 0; i < nlanes; i++ )
               sizes[i] = X[i].size( );
          ReverseSortSync( sizes, id );
          int nuse;

          if (VERBOSE)
          {    cout << "\naccepted reads:\n";
               for ( int i = 0; i < nlanes; i++ )
               {    cout << "[" << i+1 << "] " << lanes[2*i] << " --> "
                         << X[ id[i] ].size( ) << endl;    }    }

          yield = X[ id[0] ].size( );
          for ( nuse = 2; nuse <= nlanes; nuse++ )
          {    longlong new_yield = nuse * X[ id[nuse-1] ].size( );
               if ( new_yield < yield && nuse > 2 ) break;
               yield = new_yield;    }
          --nuse;
          if (VERBOSE) cout << "\n" << nuse << " lanes used\n" << endl;
          for ( int i = 0; i < nuse; i++ )
          {    RandomShuffle( X[ id[i] ] );
               X[ id[i] ].resize( X[ id[nuse-1] ].size( ) );
               Sort( X[ id[i] ] );    }

          // Find duplicates in cross set.

          for ( int i1 = 0; i1 < nuse; i1++ )
          {    for ( int j = 0; j < X[ id[i1] ].isize( ); j++ )
               {    Bool dup = False;
                    const bytes<nbytes>& K = X[ id[i1] ][j];
                    for ( int i2 = i1+1; i2 < nuse; i2++ )
                         if ( BinMember( X[ id[i2] ], K ) ) dup = True;
                    if (dup)
                    {    duplicates++;
                         if (VERBOSE)
                         {    basevector key;
                              key.assignBaseBits(nbases,K.x);
                              for ( int i = 0; i < nbases; i++ )
                                   cout << as_base( key[i] );
                              cout << " ";
                              for ( int i = nbases; i < 2*nbases; i++ )
                                   cout << as_base( key[i] );
                              cout << "\n";    }    }    }    }
          duplicates = ( duplicates * nuse ) / ( nuse - 1 );
          if (VERBOSE)
          {    cout << "estimated duplicate pairs = " << duplicates << endl;
               cout << "total pairs = " << yield << endl;
               cout << "unique fraction = "
                    << PERCENT_RATIO( 4, yield-duplicates, yield ) << endl;    }    }

     // Find total number of distinct molecules.

     vec< bytes<nbytes> > Xall;
     for ( int i = 0; i < nlanes; i++ )
          Xall.append( X[i] );
     vec< vec<String> > rows;
     vec<String> rowm, row0, rowe, row1, row2, row3, row4, row5;
     rowm.push_back( ToString( lanes_all.size( ) ), "lanes" );
     row0.push_back( ToString(nlanes), "fully processed lanes" );
     row1.push_back( ToStringAddCommas(total_pairs), "total pairs" );
     row2.push_back( ToStringAddCommas(total_PF_pairs), "total PF pairs" );
     longlong good_pairs = Xall.size( );
     row3.push_back( ToStringAddCommas(good_pairs),  "good pairs" );
     UniqueSort(Xall);
     longlong distinct_good_pairs = Xall.size( );
     row4.push_back( ToStringAddCommas(distinct_good_pairs), "distinct good pairs" );
     longlong total, distinct;
     if ( nlanes > 1 )
     {    total = yield;
          distinct = yield - duplicates;    }
     else
     {    total = good_pairs;
          distinct = distinct_good_pairs;    }
     longlong libsize = StringOfOutput( "ComputeUniqueMolecules" + ARG(n,total)
          + ARG(c,distinct) + ARG(NH,True) + ARG(QUIET,True) ).Int( );
     row5.push_back( ToStringAddCommas(libsize), "estimated library size" );
     rows.push_back( rowm, row0, rowe, row1, row2, row3, row4, row5 );
     cout << "\n";
     PrintTabular( cout, rows, 2, "rl" );
     cout << "\n";    }
