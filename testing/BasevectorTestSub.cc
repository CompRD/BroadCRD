#define FAST_ASSERT

#include "system/RunTime.h"
#include "Basevector.h"
#include "TaskTimer.h"

int main( int argc, char **argv )
{
  RunTime();

  for ( unsigned int offset1 = 0; offset1 < 4; ++offset1 )
  {
    for ( unsigned int offset2 = 0; offset2 < 4; ++offset2 )
    {
      for ( unsigned int little = 1; little < sizeof(long) * 4; ++little )
      {
        cout << "Tiny " << little << " sub " << offset1 << offset2 << endl;
        
        unsigned int big_number = (unsigned int)3 * 1000 * 1000 + offset1;
      
        basevector big( big_number + little + little);
        for ( unsigned int i = 0; i < big_number; ++i )
          big.Set( i, (i%7)%4 );
      
        for ( unsigned int i = big_number; i < big_number + little; ++i )
          big.Set( i, (i%5)%4 );
      
        for ( unsigned int i = big_number + little;
              i < big_number + little*2; ++i )
          big.Set( i, (i%7)%4 );
      
        basevector sub;
        sub.SetToSubOf( big, big_number, little );
        
        for ( unsigned int i = 0; i < sub.size(); ++i )
          ForceAssertEq( (int)sub[i], (int)big[big_number+i] );
      }
    }
  }

  for ( unsigned int offset1 = 0; offset1 < 4; ++offset1 )
  {
    for ( unsigned int offset2 = 0; offset2 < 4; ++offset2 )
    {
      cout << "Small sub " << offset1 << offset2 << endl;
      
      unsigned int big_number = (unsigned int)3 * 1000 * 1000 + offset1;
      unsigned int little_number = 1000000 + offset2;
      
      basevector big( big_number + little_number + little_number);
      for ( unsigned int i = 0; i < big_number; ++i )
        big.Set( i, (i%7)%4 );
      
      for ( unsigned int i = big_number; i < big_number + little_number; ++i )
        big.Set( i, (i%5)%4 );
      
      for ( unsigned int i = big_number + little_number;
            i < big_number + little_number*2; ++i )
        big.Set( i, (i%7)%4 );
      
      basevector sub;
      sub.SetToSubOf( big, big_number, little_number );
      
      for ( unsigned int i = 0; i < sub.size(); ++i )
        ForceAssertEq( (int)sub[i], (int)big[big_number+i] );
    }
  }
}
  

  
