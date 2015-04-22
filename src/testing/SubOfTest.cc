#include "Basevector.h"
#include "random/Random.h"
#include "TaskTimer.h"

#define TOTAL 600

int main(int argc, char **argv) {

  basevector a(TOTAL);

  for( int i = 0 ; i < TOTAL ; i++) {
    a.Set(i,randomx()%4);
  }

  TaskTimer * t = new TaskTimer();
  t->Start();

  for ( unsigned int length = 1; length <= TOTAL ; length++ ) {
    for ( unsigned int pos = 0 ; pos < TOTAL ; pos++ ) {
      if ( pos+length > TOTAL ) continue;
      for ( unsigned int extra = 0 ; extra < 17 ; extra++ ) {
	basevector b;
	b.SetToSubOf(a,pos,length, extra);

	if ( b.size() != length ) cerr << "ERROR!!!" << endl;
	basevector::const_iterator b_it = b.Begin();
	basevector::const_iterator b_it_end = b.End();
	basevector::const_iterator a_it = a.Begin(pos);	
	for ( ; b_it != b_it_end; ++b_it, ++a_it  ) {
	  if ( *b_it != *a_it ) {
	    cout << "ERROR!!! pos= " << pos << " len=" << length << endl;
	    cout << a.ToString() << endl;
	    for ( unsigned int skip = 0 ; skip < pos ; skip++ ) cout << " ";
	    cout << b.ToString() << endl;
	  }
	}
      }
    }
  }
  t->Stop();
  cout << "Task time: " << (*t) << endl;
  delete t;
}
