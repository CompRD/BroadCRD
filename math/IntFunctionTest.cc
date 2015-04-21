#include "MainTools.h"
#include "Vec.h" 

#include "math/IntFunction.h"

int main(int argc, char *argv[]) 
{
  RunTime();
  BeginCommandArguments;
  CommandArgument_Int_Doc (COV,   "Coverage.");
  
  EndCommandArguments;




  IntFunction<double> f(10, 13, 0);
  f[4] = 2.0;
  f[11] = 3.0;
  f[12] = 4.0;
  f[13] = 1.0;

  f.expand_pos_infinity();
  f.expand_neg_infinity();

  cout << "f(x_min = " << f.x_min() << ") = " << f.f_x_min() << endl;
  cout << "f(x_max = " << f.x_max() << ") = " << f.f_x_max() << endl;


  for (int x = 0; x <= 16; x++)
    cout << "f(" << x << ") = " << f[x] << endl;
  cout << endl;

  IntFunctionPrimitive<double> fp(f);

  
  for (int x = 0; x <= 20; x++)
  {
      int y = x + 2*exp(0);
      cout << "sum_f(" << x << ".." << y << ") = " << fp.f_sum(x, y) << endl;
  }

}
