// This code demonstrates the speed of fstream i/o, for purposes of evaluating other compilers.

#include "MainTools.h"
#include "TaskTimer.h"

int main( int argc, char** argv )
{
  using namespace std;

  ios::sync_with_stdio(false);

  BeginCommandArguments;
  CommandArgument_String( FILE );
  EndCommandArguments;

  ifstream in( FILE.c_str() );

  TaskTimer t;

  t.Start();
  char c[8192];
  while ( in )
    in.getline(c, 8192);
  t.Stop();
  cout << t << endl;

  t.Reset();
  in.close();
  in.clear();

  in.open( FILE.c_str() );
  ForceAssert( in.good() );

  t.Start();
  int i;
  while ( in )
    in >> i;
  t.Stop();
  cout << t << endl;
}
    
