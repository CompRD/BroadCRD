#include "MainTools.h"
#include "math/Matrix.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv) {
	RunTime();

	cout << "Testing..." << endl;

	/*
	matrix<float> m2x2(2, 2, 0.0);

	m2x2(0,0) = 3.0; m2x2(0,1) = 1.0;
	m2x2(1,0) = 2.0; m2x2(1,1) = 6.0;

	Assert(m2x2.Det() == 16);

	matrix<float> m3x3(3, 3, 0.0);

	m3x3(0,0) = 3.0; m3x3(0,1) = 1.0; m3x3(0,2) = 4.0;
	m3x3(1,0) = 2.0; m3x3(1,1) = 6.0; m3x3(1,2) = 9.0;
	m3x3(2,0) = 1.0; m3x3(2,1) = 2.0; m3x3(2,2) = 1.0;

	Assert(m3x3.Det() == -37);
	*/

	matrix<float> m4x4(4, 4, 0.0);

	m4x4(0,0) =  2.0; m4x4(0,1) =  3.0; m4x4(0,2) =  4.0; m4x4(0,3) =  5.0;
	m4x4(1,0) =  0.0; m4x4(1,1) = -1.0; m4x4(1,2) =  2.0; m4x4(1,3) =  1.0;
	m4x4(2,0) =  0.0; m4x4(2,1) =  0.0; m4x4(2,2) =  2.0; m4x4(2,3) =  4.0;
	m4x4(3,0) =  0.0; m4x4(3,1) =  3.0; m4x4(3,2) = -6.0; m4x4(3,3) =  0.0;

	/*
	Assert(m4x4.Det() == -12);

	m4x4.PrettyPrint(cout);
	cout << m4x4.Cofactor(0,0) << " " << m4x4.Cofactor(0,1) << " " << m4x4.Cofactor(0,2) << " " << m4x4.Cofactor(0,3) << endl;
	cout << m4x4.Cofactor(1,0) << " " << m4x4.Cofactor(1,1) << " " << m4x4.Cofactor(1,2) << " " << m4x4.Cofactor(1,3) << endl;
	cout << m4x4.Cofactor(2,0) << " " << m4x4.Cofactor(2,1) << " " << m4x4.Cofactor(2,2) << " " << m4x4.Cofactor(2,3) << endl;
	cout << m4x4.Cofactor(3,0) << " " << m4x4.Cofactor(3,1) << " " << m4x4.Cofactor(3,2) << " " << m4x4.Cofactor(3,3) << endl;

	m4x4.Invert();
	m4x4.PrettyPrint(cout);
	*/

	matrix<float> m4x4_inv = m4x4, m4x4_identity;
	m4x4_inv.Inverse(m4x4_inv);
	mul(m4x4_inv, m4x4, m4x4_identity);
	m4x4_identity.PrettyPrint(cout);

	cout << "Passed!" << endl;

	return 0;
}
