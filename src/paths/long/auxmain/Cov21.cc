///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Report coverage by 21-mers from Fosmid pool reads of Fosmid 40 at kmers that
// are unique in the Fosmid.

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"

int main( )
{
     RunTime( );
     String dir = "/wga/scr4/jaffe/CompareVars/40";
     vecbasevector bases( dir + "/tmp.fos/frag_reads_orig.fastb" );
     vecbasevector ref;
     FetchReads( ref, 0, dir + "/fos.40.fasta" );
     vec< pair<basevector,int> > kmers_plus;
     const int K = 21;
     for ( int p = 0; p <= ref[0].isize( ) - K; p++ )
     {    basevector b( ref[0], p, K );
          kmers_plus.push( b, p + K/2 );    }
     Sort(kmers_plus);
     vec<basevector> kmers0( kmers_plus.size( ) );
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
          kmers0[i] = kmers_plus[i].first;
     vec<Bool> to_delete( kmers_plus.size( ), False );
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    basevector b = kmers_plus[i].first;
          b.ReverseComplement( );
          Bool have_rc = BinMember( kmers0, b );
          int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          if ( j - i > 1 || have_rc )
          {    for ( int k = i; k < j; k++ )
                    to_delete[k] = True;    }
          i = j - 1;    }
     EraseIf( kmers_plus, to_delete );
     vec<basevector> kmers( kmers_plus.size( ) );
     vec<int> pos( kmers_plus.size( ) );
     for ( int i = 0; i < kmers.isize( ); i++ )
     {    kmers[i] = kmers_plus[i].first;
          pos[i] = kmers_plus[i].second;    }
     Sort(pos);
     vec<int> cov( ref[0].size( ), 0 );
     for ( int i = 0; i < (int) bases.size( ); i++ )
     for ( int j = 0; j <= bases[i].isize( ) - K; j++ )
     {    basevector b( bases[i], j, K );
          int p = BinPosition( kmers, b );
          if ( p >= 0 ) cov[ kmers_plus[p].second ]++;
          b.ReverseComplement( );
          p = BinPosition( kmers, b );
          if ( p >= 0 ) cov[ kmers_plus[p].second ]++;    }
     double total = 0;
     double n = 0;
     const int rl = 250;
     for ( int p = rl; p < ref[0].isize( ) - rl; p++ )
     {    if ( BinMember( pos, p ) )
          {    total += cov[p];
               n++;    }    }
     vec<int> depth(1000, 0);
     const int flank = 10;
     for ( int p = rl; p < ref[0].isize( ) - rl; p++ )
     {    if ( BinMember( pos, p ) )
          {    double d = cov[p]/(total/n);
               cout << p << "  " << cov[p] << "  " << d;
               cout << "  ";
               for ( int j = 0; j < flank; j++ )
                    cout << as_base( ref[0][p-flank+j] );
               cout << "/" << as_base( ref[0][p] ) << "/";
               for ( int j = 0; j < flank; j++ )
                    cout << as_base( ref[0][p+j+1] );
               cout << endl;    
               depth[ int(round(10*d)) ]++;    }    }
     cout << "\n";
     int top = 20;
     for ( int i = 0; i < top; i++ )
     {    if ( depth[i] > 0 )
          {    cout << i/10.0 << "-" << (i+1)/10.0 << "  " << depth[i] 
                    << "  " << 100 * depth[i]/n << "%" << endl;    }    }
     double tail = 0;
     for ( int i = top; i < depth.isize( ); i++ )
          tail += depth[i];
     cout << top/10.0 << "-infinity  " << tail << "  " << 100 * tail/n 
          << "%" << endl;    }
