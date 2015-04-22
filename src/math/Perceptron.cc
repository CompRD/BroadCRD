/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Perceptron: voting perceptron from 
//      "Large Margin Classification Using the Perceptron Algorithm"
//      Machine Learning, 37(3):277-296, 1999.
//
//      http://www.cs.ucsd.edu/~yfreund/papers/LargeMarginsUsingPerceptron.pdf
//
//      In practice, almost as accurate as SVM, waaaay cheaper.

#include <stdlib.h>
#include <cassert>
#include <fstream>
#include "Basevector.h"
#include "MainTools.h"
#include "Map.h"

vec<char> separators;

void ReadSample(String const& line, int* y, StdMap<int,float>& X)
{
    vec<String> tokens;
    Tokenize(line, separators, tokens);

    *y = tokens[0].Int();
    for (int i = 1; i < tokens.isize(); i++)
    {
        int d = tokens[i].Before(":").Int();
        float x = (float)(tokens[i].After(":").Double());
        X[d] = x;
    }
}

int Predict(StdMap<int,float>& V, StdMap<int,float>& X)
{
    float score = 0.0;

    for (map<int,float>::iterator i = V.begin(); i != V.end(); i++)
    {
        int   d = (*i).first;
        float v = (*i).second;

        score += v * X[d]; 
    }    

    if (score >= 0) { return  1; }
    if (score <  0) { return -1; }
    return 0;
}

StdMap<int,float> Update(StdMap<int,float>& V, int y, StdMap<int,float>& X)
{
    StdMap<int,float> V2;

    for (auto i = V.begin(); i != V.end(); i++)
    {
        int   d = (*i).first;
        float v = (*i).second;
        V2[d] = v + (y*X[d]); 
    }

    for (auto i = X.begin(); i != X.end(); i++)
    {
        int   d = (*i).first;
        float v = (*i).second;
        V2[d] = V2[d] + (y*X[d]); 
    }

    return V2;
}

int Classify(vec<float>& c, vec< StdMap<int,float> >& V, StdMap<int,float>& X)
{
    int s = 0;
    for (int k = 0; k < V.isize(); k++)
    {
        float sum = 0.0;
        for (map<int,float>::iterator i = V[k].begin(); i != V[k].end(); i++)
        {
            int   d = (*i).first;
            float v = (*i).second;

            sum += v * X[d]; 
        }

        if (sum >= 0) { s += (int)( 1 * c[k]); }
        if (sum <  0) { s += (int)(-1 * c[k]); }
    }

    if (s >= 0) { s =  1; }
    if (s <  0) { s = -1; }

    return s;
}

float Classify_real(vec<float>& c, vec< StdMap<int,float> >& V, StdMap<int,float>& X)
{
    float s = 0;
    for (int k = 0; k < V.isize(); k++)
    {
        float sum = 0.0;
        for (map<int,float>::iterator i = V[k].begin(); i != V[k].end(); i++)
        {
            int   d = (*i).first;
            float v = (*i).second;

            sum += v * X[d]; 
        }

        if (sum >= 0) { s +=  1 * c[k]; }
        if (sum <  0) { s += -1 * c[k]; }
    }

    s = s / (float)(V.size());

    return s;
}

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String_OrDefault(TRAINING, "");
  CommandArgument_String_OrDefault(TEST, "");
  CommandArgument_String(MODEL);
  CommandArgument_Int_OrDefault(T, 20);
  CommandArgument_Bool_OrDefault(REAL_OUTPUT, false);
  EndCommandArguments;

  separators.clear( );
  separators.push_back( ' ' );
  separators.push_back( '\t' );
  separators.push_back( '\n' );

    if (TRAINING != "")
	{
            std::ofstream model(MODEL.c_str());
	
	    int k = 0;
	    vec< StdMap<int,float> > V(1);
	    vec<float> c(1);
	
	    for (int t = 0; t < T; t++)
	    {
                std::ifstream training(TRAINING.c_str());
	
	        StdMap<int,float> X;
	        int            y;
	        String line;
	        while ( getline(training,line) )
	        {
	            ReadSample(line,&y,X);
	            int z = Predict(V[k], X);
	            if (z == y) 
	            { 
	                // Correct.
	                c[k] += 1; 
	            }
	            else
	            {
	                // Incorrect.
	                V.push_back(Update(V[k], y, X));
	                c.push_back(1);
	                k += 1;
	            }
	        }
	        ForceAssert(training.eof());
	    }
	
	    for (int k = 0; k < V.isize(); k++)
	    {
	        model << (int)c[k] << ' ';
	        for (map<int,float>::iterator i = V[k].begin(); i != V[k].end(); i++)
	        {
	            int   d = (*i).first;
	            float v = (*i).second;
	
	            model << d << ':' << v << ' ';
	        }
	        model << '\n';
	    }
	    model.close();
	    ForceAssert(model);
	}
    else if (TEST != "")
    {
	    int k;
	    vec< StdMap<int,float> > V;
	    vec<float> c;

	std::ifstream model(MODEL.c_str());
        int weight;
        StdMap<int,float> v;
        String line;
        while( getline(model,line) )
        {
            ReadSample(line, &weight, v);
            c.push_back(weight);
            V.push_back(v);
        }
        ForceAssert(model.eof());

        std::ifstream test(TEST.c_str());
        int y;
        StdMap<int,float> X;

        while( getline(test,line) )
        {
            ReadSample(line, &y, X);
            if (! REAL_OUTPUT) 
            {
                int z = Classify(c, V, X);
                printf("%d\n", z);
            }
            else
            {
                float z = Classify_real(c, V, X);
                printf("%f\n", z);
            }

        }
        ForceAssert(test.eof());
    }
}





