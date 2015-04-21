// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// Program: CorrectByKmers
// 
// Given a fastb file of reads, and the table T of nonunique
// k-mers in them (together with their multiplicities, as produced by 
// FindNonuniqueKmers), attempt to correct sequencing errors in the reads using
// the following experimental scheme.
//
// Scan each read R starting at the left.  Take the k-mer x starting at the
// given position, and binary search T.
//     If x is not in T, look up all 3k of the point mutations of x in T, and if 
// exactly one is present, and has multiplicity >= 5, immediately edit R to 
// incorporate the corresponding correction.
//     Then proceed to the next k-mer in R.
//
// If SHOW_MISTAKES=True, find and mark any error-correction events that in fact
// introduce errors.  The file reads.true.fastb must exist to use this option.
// If SHOW_CHANGES=True, show all changes, including mistakes.

#include "Basevector.h"
#include "kmers/KmerRecord.h"
#include "MainTools.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"

template<int K> inline longlong GetPos( const kmer_with_count<K>& x, 
     const vec< kmer_with_count<K> >& kmers, const vec<longlong>& index, 
     int ishift )
{    unsigned int xtop = *x.Ints( );
     unsigned int xshift = xtop >> ishift;
     longlong pos = lower_bound( kmers.begin( ) + index[xshift], 
          kmers.begin( ) + index[xshift+1], x ) - kmers.begin( );
     return pos;    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA); 
     CommandArgument_String(RUN); 
     CommandArgument_Int(K);
     CommandArgument_Int_OrDefault(K0, 14);
     CommandArgument_Int_OrDefault(MULT, 5);
     CommandArgument_String(READS_IN);
     CommandArgument_String(KMERS);
     CommandArgument_String_OrDefault(READS_OUT, "");
     CommandArgument_Bool_OrDefault(SHOW_MISTAKES, False);
     CommandArgument_Bool_OrDefault(SHOW_CHANGES, False);
     EndCommandArguments; 

     // Check essential input requirements.

     ForceAssertLe( K0, 16 );
     ForceAssertLe( K0, K );
     ForceAssertGe( K, 16 );

     // Check special input requirements, imposed only because code not implemented
     // for all K.

     ForceAssertEq( K, 20 );

     // Load data.

     String run_dir = PRE + "/" + DATA + "/" + RUN;
     vecbasevector reads( run_dir + "/" + READS_IN );
     vecbasevector true_reads;
     if ( SHOW_MISTAKES || SHOW_CHANGES ) 
          true_reads.ReadAll( run_dir + "/reads.true.fastb" );
     vec< kmer_with_count<20> > kmers;
     BinaryReader::readFile( run_dir + "/" + KMERS, &kmers );

     // Set up index to kmers using the first K0 bases.

     longlong four_to_K0 = 1;
     for ( int i = 0; i < K0; i++ )
          four_to_K0 *= 4;
     vec<longlong> index( four_to_K0 + 1, -1 );
     int ishift = 2 * (16 - K0);
     for ( longlong i = kmers.size( ) - 1; i >= 0; i-- )
     {    static basevector b(K0);
          unsigned int x = *kmers[i].Ints( );
          index[ x >> ishift ] = i;    }
     index[four_to_K0] = kmers.size( );
     for ( longlong i = four_to_K0 - 1; i >= 0; i-- )
          if ( index[i] == -1 ) index[i] = index[i+1];

     // Make changes.

     longlong changes = 0;
     double clock = -1.0;
     for ( size_t i = 0; i < reads.size( ); i++ )
     {    basevector& r = reads[i];
          if ( i % 100000 == 0 ) 
          {    if ( i > 0 )
               {    double time_used = WallClockTime( ) - clock;
                    cout << "i = " << i << " of " << reads.size( ) << ", changes = "
                         << changes << ", time used = " << setprecision(3) 
                         << time_used << " seconds" << endl;    }
               clock = WallClockTime( );    }
          for ( int j = 0; j <= r.isize( ) - K; j++ )
          {    static basevector b, m;
               b.SetToSubOf( r, j, K );
               m = b;
               CanonicalizeKmer(m);
               kmer_with_count<20> x( m, 0 );
               longlong pos = GetPos( x, kmers, index, ishift );
               if ( pos == kmers.isize( ) ) continue;
               if ( eq_kmer( x, kmers[pos] ) ) continue;
               int uhit = -1, vhit = -1, counthit = -1;
               for ( int u = 0; u < K; u++ )
               {    for ( int v = 0; v < 4; v++ )
                    {    if ( b[u] == v ) continue;
                         m = b;
                         m.Set( u, v );
                         CanonicalizeKmer(m);
                         kmer_with_count<20> y( m, 0 );
                         longlong mpos = GetPos( y, kmers, index, ishift );
                         if ( mpos == kmers.isize( ) ) continue;
                         if ( !eq_kmer( y, kmers[mpos] ) ) continue;
                         int count = kmers[mpos].Count( );
                         if ( count < MULT || uhit >= 0 ) goto next_pos;
                         uhit = u;
                         vhit = v;    
                         counthit = count;    }    }
               if ( uhit >= 0 ) 
               {    if ( SHOW_MISTAKES || SHOW_CHANGES )
                    {    Bool wrong = ( true_reads[i][j+uhit] != vhit );
                         if ( wrong || SHOW_CHANGES )
                         {    cout << "\n";
                              cout << ( wrong ? "Introducing" : "Correcting" );
                              cout << " error in read " << i 
                                   << " using kmer at bases [" << j << "," << j+K
                                   << ") by changing base " << j+uhit
                                   << " to " << as_base(vhit) << ",\n"
                                   << "yielding a kmer having multiplicity "
                                   << counthit << ".\n";
                              cout << "The current base is " << as_base(r[j+uhit]) 
                                   << " and the true base is " 
                                   << as_base(true_reads[i][j+uhit]) << ".\n";
                              cout << "That read has errors at bases";
                              for ( int u = 0; u < r.isize( ); u++ )
                                   if ( r[u] != true_reads[i][u] ) cout << " " << u;
                              cout << ".\n";    }    }    
                    r.Set( j + uhit, vhit );
                    ++changes;    }
               next_pos: continue;    }    }

     // Save results.

     PRINT(changes);
     if ( READS_OUT != "" ) reads.WriteAll( run_dir + "/" + READS_OUT );    }
