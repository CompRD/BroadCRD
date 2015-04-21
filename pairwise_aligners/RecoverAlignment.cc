// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include "Alignment.h"
#include "system/Assert.h"
#include "Basevector.h"
#include "kmers/KmerShape.h"
#include "CoreTools.h"
#include "pairwise_aligners/MakeAligns.h"
#include "Overlap.h"
#include "Qualvector.h"
#include "ScoreAlignment.h"
#include "ShortVector.h"

// RecoverAlignment: try to find an alignment between EE[id1] and EE[id2] with
// orientation rc; return the answer as ap if successful.

Bool RecoverAlignment( int id1, int id2, const vecbasevector& EE,
     const vecqualvector& Q, Bool rc, alignment_plus& ap,
     int kmer_size, int min_offset, int max_offset )
{    temp_file aligns_file( "/tmp/RecoverAlignment_tmp_XXXXXXX" );
     vecbasevector EE2(2), EE2rc(2);
     vecqualvector Q2(2), Q2rc(2);
     EE2[0] = EE2rc[0] = EE[id1];
     EE2[1] = EE2rc[1] = EE[id2];
     EE2rc[0].ReverseComplement();
     EE2rc[1].ReverseComplement();
     Q2[0] = Q[id1];
     Q2[1] = Q[id2];
     Q2rc[0] = Reverse( Q[id1] );
     Q2rc[1] = Reverse( Q[id2] );
     Ofstream( trash, "/dev/null" );

     Assert( kmer_size == 8 || kmer_size == 12 );
     // if ( min_offset < -100000 && max_offset > 100000 )
     #define CALL_MAKE_ALIGNS(_K)                                              \
     {    if ( _K == 8 )                                                       \
               MakeAligns<2, _K, 50>( 1, 1, EE2,                               \
                    to_compare(FIRST_VS_SECOND, 1), 1000, 1000, aligns_file,   \
                    trash, 200, 2500, 0, 100, 7, 9, False, 3 );                \
          else MakeAligns<2, _K, 50>( 1, 1, EE2,                               \
                    to_compare(FIRST_VS_SECOND, 1), 1000, 1000, aligns_file,   \
                    trash, 200, 255, 0, 100, 3, 4 );    }

     DISPATCH_ON_K( kmer_size, CALL_MAKE_ALIGNS );
     /*
     else
     {    vec< vec< vec<int> > > allowed_offsets(1);
          allowed_offsets[0].resize(1);
          allowed_offsets[0][0].push_back( (min_offset + max_offset)/2 );
          int max_offset_discrep = (max_offset - min_offset)/2 + 1;
          if ( kmer_size == 8 )
               MakeAligns<2, 8, 50>( 1, 1, EE2, EE2rc, 
                    to_compare(FIRST_VS_SECOND, 1), 1000, 1000, aligns_file, 
                    trash, 200, 2500, 0, 100, 7, 9, False, 3, 0, False,
                    allowed_offsets, max_offset_discrep );
          else MakeAligns<2, 12, 50>( 1, 1, EE2, EE2rc, 
                    to_compare(FIRST_VS_SECOND, 1), 1000, 1000, aligns_file, 
                    trash, 200, 255, 0, 100, 3, 4, False, 1, 0, False,
                    allowed_offsets, max_offset_discrep );    }
     */

     int aligns_length2;
     Ifstream( aligns_in2, aligns_file );
     int xi1, xi2, best_i2 = -1;
     Bool rc2 = False;
     float bad, best_bad2 = 0;
     while(1)
     {    BinRead( aligns_in2, aligns_length2 );
          if ( !aligns_in2 ) break;
          avector<alignment> aligns2( aligns_length2 );
          for ( int i2 = 0; i2 < aligns_length2; i2++ )
               aligns2(i2).Read( aligns_in2, xi1, xi2, rc2 );
          if ( rc != rc2 ) continue;
          Assert( xi1 == 0 && xi2 == 1 );
          for ( int i2 = 0; i2 < aligns_length2; i2++ )
          {    alignment& a = aligns2(i2);
               if ( EstimatedOverlap(a, EE2[0], EE2[1]) == 0 ) continue;
               int offset = a.pos1( ) - a.pos2( );
               if ( offset < min_offset || offset > max_offset ) continue;
               if ( !rc )
               {    Regap( a, EE2[0], Q2[0], EE2[1], Q2[1] );
                    bad = ScoreAlignment( a, EE2[0], Q2[0], EE2[1], Q2[1] );    }
               else
               {    Regap( a, EE2[0], Q2[0], EE2rc[1], Q2rc[1] );
                    bad = ScoreAlignment( a, EE2[0], Q2[0], EE2rc[1], Q2rc[1] );    }
               if ( best_i2 == -1 || bad < best_bad2 )
               {    best_i2 = i2;
                    best_bad2 = bad;    }    }
          if ( best_i2 >= 0 )
          {    ap = alignment_plus( id1, id2, EE[id1].size( ), EE[id2].size( ), 
                    rc2, aligns2(best_i2), best_bad2 );
               return True;    }
          return False;    }
     return False;    }
