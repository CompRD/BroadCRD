// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include <fstream>

#include "Nobbit.h"
#include "PackAlign.h"
#include "system/System.h"

ostream& operator<<( ostream& out, const nobbit& n )
{    BinWrite( out, n.pos1 );
     BinWrite( out, n.pos2 );
     int nblocks = n.a.Nblocks( );
     BinWrite( out, nblocks );
     for ( int i = 0; i < n.a.Nblocks( ); i++ )
     {    int gap = n.a.Gaps(i), length = n.a.Lengths(i);
          BinWrite( out, gap );
          BinWrite( out, length );    }
     BinWrite( out, n.RC );
     BinWrite( out, n.length1 );
     BinWrite( out, n.length2 );
     BinWrite( out, n.id1 );
     BinWrite( out, n.id2 );
     BinWrite( out, n.score );
     BinWrite( out, n.Pos1 );
     BinWrite( out, n.Pos2 );
     return out;    }

istream& operator>>( istream& in, nobbit& n )
{    int gap, len, nblocks;
     BinRead( in, n.pos1 );
     n.a.Setpos1(n.pos1);
     BinRead( in, n.pos2 );
     n.a.Setpos2(n.pos2);
     BinRead( in, nblocks );
     n.a.SetNblocks(nblocks);
     for ( int i = 0; i < nblocks; i++ )
     {    BinRead( in, gap );
          n.a.SetGap( i, gap );
          BinRead( in, len );
          n.a.SetLength( i, len );    }
     BinRead( in, n.RC );
     BinRead( in, n.length1 );
     BinRead( in, n.length2 );
     BinRead( in, n.id1 );
     BinRead( in, n.id2 );
     BinRead( in, n.score );
     BinRead( in, n.Pos1 );
     BinRead( in, n.Pos2 );
     return in;    }

void nobbit::Print( ostream& out )
{
  out << id1 << "\t"
      << id2 << "\t"
      << length1 << "\t"
      << length2 << "\t"
      << pos1 << "\t"
      << Pos1 << "\t"
      << pos2 << "\t"
      << Pos2 << "\t"
      << RC << "\t"
      << score << "\n";

}
