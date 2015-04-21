/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "solexa/FourBase.h"

int four_base::Call( ) const
{
    int result = 0;
    intens_t max = A();
    if ( C() > max )
    {
        result = 1;
        max = C();
    }
    if ( G() > max )
    {
        result = 2;
        max = G();
    }
    if ( T() > max )
    {
        result = 3;
    }
    return result;
}

float four_base::CallQuality( ) const
{    int call = Call( );
     float call_value = base(call), next_base = 0.0;
     if ( call_value == 0 ) return 1;
     for ( int i = 0; i < 4; i++ )
     {    if ( i == call ) continue;
          next_base = ::Max( next_base, base(i) );    }
     if ( next_base == 0 ) return 1000000000;
     return call_value / next_base;    }

ostream& operator<<( ostream& out, const four_base& b )
{    return out << "[A=" << b.A( ) << ",C=" << b.C( ) << ",G=" << b.G( )
          << ",T=" << b.T( ) << "]";    }

void Call( const VecFourBaseVec& I, vecbasevector& bases )
{    longlong nbases = 0;
     int nseqs = I.size( );
     for ( VecFourBaseVec::size_type i = 0; i < I.size( ); i++ )
          nbases += I[i].size( );
     bases.clear( );
     bases.Reserve( nbases/16 + nseqs, nseqs );
     for ( VecFourBaseVec::size_type i = 0; i < I.size( ); i++ )
     {    static basevector b;
          b.resize( I[i].size( ) );
          for ( int j = 0; j < b.isize( ); j++ )
               b.Set( j, I[i][j].Call( ) );
          bases.push_back(b);    }    }


DeconvolveIntensities::DeconvolveIntensities(String filename) :
  m(filename.empty() ? 0 : 4, filename.empty() ? 0 : 4),
  permut(4)
{
  if (!filename.empty()) {
    // Read a 4x4 matrix from a Solexa Firecrest run.  Here is an example:
    //
    //# Auto-generated frequency response matrix
    //> A
    //> C
    //> G
    //> T
    //0.78 0.11 -0.02 -0.02 
    //0.91 0.95 -0.02 -0.01 
    //-0.06 -0.05 1.08 -0.02 
    //-0.05 -0.05 0.97 1.44 
    //
	/*
    Ifstream(in, filename);
    string line;
    int i, j;
    // Skip 5-line header
    for (i=0; i<5; ++i)
      getline(in, line);
    // Read 16 values into m
    for (i=0; i<4; ++i) {
      for (j=0; j<4; ++j) {
	in >> m(i,j);
      }
    }
	*/

	ifstream in(filename.c_str());

	string line;
	int r = 0;

	while (getline(in, line)) {
		String sline(line);

		if (r < 4 && !sline.Contains("#") && !sline.Contains(">")) {
			istringstream ss(line);

			for (int c = 0; c < 4; c++) {
				ss >> m(r, c);
			}

			r++;
		}
	}

    //if (!in)
      //FatalErr("Failed to read Solexa matrix from " << filename);

    cout << "Read intensity crosstalk matrix:\n";
    m.PrettyPrint(cout);

   // prepare to solve equations mx=b for x--equivalent to inverting
    // m then computing inv(m) * b, but better numerical performance.
    m.lu_decompose(permut);
  }
}

#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
template class SmallVec< four_base, MempoolAllocator<four_base> >;
template class OuterVec<FourBaseVec>;
