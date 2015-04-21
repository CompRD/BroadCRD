///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Get background rate of placements.
// Computation supposed to mimic computation in FindTenxNhood.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/tenx/TenxDirs.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_Doc(N, "1 or 2 or 3 or 4 or 5 or 6");
     EndCommandArguments;

     // Hardcoded directories.

     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Load counts and assembly.

     cout << "\n" << Date( ) << ": loading assembly" << endl;
     HyperBasevectorX hb;
     vec<int> inv, lens;
     vec<vec<vec<vec<int>>>> lines;
     #pragma omp parallel sections
     {    
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.hbx", &hb );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.inv", &inv );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.lines", &lines );    }    }
     GetLineLengths( hb, lines, lens );

     // Load alignments and barcodes.

     cout << Date( ) << ": loading alignments" << endl;
     vec< pair<int,int> > places;
     BinaryReader::readFile( odir + "/10X.aligns", &places );
     cout << Date( ) << ": loading barcodes" << endl;
     vec<uint32_t> bcs;
     BinaryReader::readFile( tdir + "/10X.bc", &bcs );
     vec<int64_t> bci;
     BinaryReader::readFile( tdir + "/10X.bci", &bci );
     int nbc_total = bci.size( ) - 1;

     // Create index to lines.

     cout << Date( ) << ": making tol4" << endl;
     vec< quad<int,int,int,int> > tol4( hb.E( ) );
     for ( int i = 0; i < lines.isize( ); i++ )
     for ( int j = 0; j < lines[i].isize( ); j++ )
     for ( int k = 0; k < lines[i][j].isize( ); k++ )
     for ( int l = 0; l < lines[i][j][k].isize( ); l++ )
          tol4[ lines[i][j][k][l] ] = make_quad( i, j, k, l );

     // Precompute medians.

     cout << Date( ) << ": precomputing medians" << endl;
     vec<vec<int>> medians( lines.size( ) );
     for ( int lid = 0; lid < lines.isize( ); lid++ )
     {    const vec<vec<vec<int>>>& L = lines[lid];
          medians[lid].resize( L.size( ), 0 );
          for ( int m = 0; m < L.isize( ); m++ )
          {    if ( L[m].empty( ) ) continue;
               vec<int> lens;
               for ( int r = 0; r < L[m].isize( ); r++ )
               {    int n = 0;
                    for ( int u = 0; u < L[m][r].isize( ); u++ )
                         n += hb.Kmers( L[m][r][u] );
                    lens.push_back(n);    }
               Sort(lens);
               medians[lid][m] = Median(lens);    }    }

     // Mark edges that are deep within a line.  

     cout << Date( ) << ": looking deep" << endl;
     vec<Bool> deep( hb.E( ), False );
     for ( int lid = 0; lid < lines.isize( ); lid++ )
     {    const vec<vec<vec<int>>>& L = lines[lid];
          int i1, i2, d = 0;
          for ( i1 = 0; i1 < L.isize( ); i1++ )
          {    d += medians[lid][i1];
               if ( d >= 20000 ) break;    }
          d = 0;
          for ( i2 = L.isize( ) - 1; i2 >= 0; i2-- )
          {    d += medians[lid][i2];
               if ( d >= 20000 ) break;    }
          for ( int m = i1 + 1; m < i2; m++ )
          {    for ( int r = 0; r < L[m].isize( ); r++ )
               for ( int u = 0; u < L[m][r].isize( ); u++ )
                    deep[ L[m][r][u] ] = True;    }    }
     PRINT( Sum(deep) );

     // Precompute lefts and rights.

     cout << Date( ) << ": precomputing lefts and rights" << endl;
     vec<int> lefts( hb.E( ), 0 ), rights( hb.E( ), 0 );
     for ( int eid = 0; eid < hb.E( ); eid++ )
     {    if ( deep[eid] ) continue;
          int lid = tol4[eid].first;
          const vec<vec<vec<int>>>& L = lines[lid];
          int j = tol4[eid].second, k = tol4[eid].third, l = tol4[eid].fourth;
          for ( int m = 0; m < j; m++ )
               lefts[eid] += medians[lid][m];
          for ( int u = 0; u < l; u++ )
               lefts[eid] += hb.Kmers( L[j][k][u] );
          for ( int m = j+1; m < L.isize( ); m++ )
               rights[eid] += medians[lid][m];
          for ( int u = l+1; u < L[j][k].isize( ); u++ )
               rights[eid] += hb.Kmers( L[j][k][u] );    }

     // Sort.  For unclear reasons, the ParallelSortSync below is very slow.
     // Therefore there is other code below it to replace it.

     vec<int64_t> ids( places.size( ), vec<int64_t>::IDENTITY );
     // ParallelSortSync( places, ids );
     cout << Date( ) << ": building placesz" << endl;
     vec< triple<int,int,int64_t> > placesz( places.size( ) );
     #pragma omp parallel for
     for ( int64_t i = 0; i < places.jsize( ); i++ )
     {    placesz[i].first = places[i].first;
          placesz[i].second = places[i].second;
          placesz[i].third = i;    }
     cout << Date( ) << ": sorting" << endl;
     ParallelSort(placesz);
     cout << Date( ) << ": building places and ids" << endl;
     for ( int64_t i = 0; i < places.jsize( ); i++ )
     {    places[i].first = placesz[i].first;
          places[i].second = placesz[i].second;
          ids[i] = placesz[i].third;    }

     // Find line hits.

     cout << Date( ) << ": finding line hits" << endl;
     vec<int> lhits( lines.size( ), 0 );
     vec< vec<int> > lhitsb( lines.size( ) );
     int64_t count = 0;
     const int64_t batch = 100000;
     int64_t batches = ( places.jsize( ) + batch - 1 ) / batch;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t b = 0; b < batches; b++ )
     {    
          #pragma omp critical
          {    if ( count % 50 == 0 ) DPRINT2( count, batches );
               count++;    }
          int64_t rid0 = b * batch;
          vec<int> bz;
          for ( int64_t xrid = rid0; xrid < Min( rid0 + batch, places.jsize( ) );
               xrid++ )
          {    int eid = places[xrid].first;
               if ( eid < 0 || deep[eid] ) continue;
               int lid = tol4[eid].first;
               int64_t z;
               for ( z = xrid + 1; z < Min( rid0 + batch, places.jsize( ) ); z++ )
                    if ( tol4[places[z].first].first != lid ) break;

               // Find the barcodes that are within 20 kb of an end.

               bz.clear( );
               for ( int64_t ridx = xrid; ridx < z; ridx++ )
               {    int eidx = places[ridx].first;
                    int left = places[ridx].second + lefts[eidx];
                    int right = hb.Bases(eidx) 
                         - places[ridx].second - 88 + rights[eidx];
                    if ( !( left > 20000 && right > 20000 ) ) 
                         bz.push_back( bcs[ids[ridx]] );    }
               UniqueSort(bz);
                         
               // Save.

               if ( bz.nonempty( ) )
               {    int lid2 = tol4[ inv[eid] ].first;
                    #pragma omp critical
                    {    lhitsb[lid].append(bz);
                         lhitsb[lid2].append(bz);   }    }    
               xrid = z - 1;    }    }

     // Count.

     cout << Date( ) << ": counting" << endl;
     for ( int l = 0; l < lines.isize( ); l++ )
     {    UniqueSort( lhitsb[l] );
          lhits[l] = lhitsb[l].size( );    }

     // Index: set lbc[b] = {lines} hitting barcode b.

     cout << Date( ) << ": indexing" << endl;
     vec<vec<int>> lbc(nbc_total);
     for ( int l = 0; l < lines.isize( ); l++ )
     for ( int i = 0; i < lhitsb[l].isize( ); i++ )
          lbc[ lhitsb[l][i] ].push_back(l);

     // Get fractions.

     cout << Date( ) << ": getting fractions" << endl;
     int64_t total_seen = 0;
     vec<int> seen( lines.size( ), 0 );
     for ( int l = 0; l < lines.isize( ); l++ )
     {    total_seen += lhitsb[l].size( );
          for ( int j = 0; j < lhitsb[l].isize( ); j++ )
          {    int b = lhitsb[l][j];
               for ( int i = 0; i < lbc[b].isize( ); i++ )
                    seen[ lbc[b][i] ]++;    }    }
     cout << Date( ) << ": setting fractions" << endl;
     vec<double> frac( lines.size( ) );
     for ( int l = 0; l < lines.isize( ); l++ )
          frac[l] = seen[l] / double(total_seen);

     // Write big files.

     cout << Date( ) << ": writing big files" << endl;
     BinaryWriter::writeFile( odir + "/10X.lhitsb", lhitsb );
     BinaryWriter::writeFile( odir + "/10X.lbc", lbc );
     
     // Write.
     
     cout << Date( ) << ": writing little files" << endl;
     BinaryWriter::writeFile( odir + "/10X.line_hits", lhits );
     BinaryWriter::writeFile( odir + "/10X.bc_frac", frac );    }
