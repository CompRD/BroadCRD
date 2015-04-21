/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "MainTools.h"
#include "FastIfstream.h"
#include "graphics/BuildHistogram.h"
#include "graphics/BarGraph.h"

#ifndef NAN
#error NAN is not defined
#endif

const char *DOC =

   "Generates histogram plot (eps file)";

// See also PlotHistogramB (slightly different functionality).
// sante - 2010.12.07


bool ToInt( const char * ptr_b , const char *ptr_e, int &result ) {
  if ( ptr_e == ptr_b ) {
    cout << "ToInt::Empty field found" << endl;
    return false;
  }
  int sign = 1;
  if ( *ptr_b == '+' || *ptr_b == '-' ) {
    if ( *ptr_b == '-' ) sign = (-1);
    ptr_b++;
    if ( ptr_e == ptr_b ) return false;
  }

  result = 0;
  for ( ; ptr_b < ptr_e ; ++ptr_b ) {
    if ( *ptr_b >= '0' && *ptr_b <= '9' ) result = result*10+static_cast<int>((*ptr_b) - '0');
    else return false;
  }
  result *= sign;
  return true;
}

float ToFloat( const char * ptr_b, const char *ptr_e ) {
  float result = 0.0f;
  float f1 = 10.0f;
  float f2 = 1.0f;

  if ( ptr_e == ptr_b ) {
    cout << "Empty field found" << endl;
    return NAN;
  }

  bool has_point = false;

  int sign = 1;
  if ( *ptr_b == '+' || *ptr_b == '-' ) {
    if ( *ptr_b == '-' ) sign = (-1);
    ptr_b++;
    if ( ptr_e == ptr_b ) return NAN;
  }

  for ( ; ptr_b < ptr_e ; ++ptr_b) {
    if ( *ptr_b >= '0' && *ptr_b <= '9' ) {
      if ( has_point ) f2 = f2*10.0f;
      result = result*f1+static_cast<float>((*ptr_b) - '0')/f2;
    } else {
      if ( *ptr_b == '.' ) {
	has_point = true;
	f1 = 1.0f;
      } else {
	if ( *ptr_b == 'e' || *ptr_b == 'E' ) {
	  int exponent;
	  if ( ToInt(ptr_b+1, ptr_e, exponent) ) {
	    result *= powf(10.0f,static_cast<float>(exponent));
	    return ( result * sign );
	  } else return NAN;
	} else return NAN;
      }
    }
  }
  return result*sign;
}

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandDoc(DOC);
     CommandArgument_String_Doc(DATA, "Space and/or TAB separated data file; if DATA=- reads from stdin");
     CommandArgument_String_Abbr_OrDefault_Doc(PATTERN,P,"","Only lines matching the pattern will be used");
     CommandArgument_Int_OrDefault_Doc(COL,0,"Column in the data file to plot histogram for; zero-based");
     CommandArgument_Int_OrDefault_Doc(NBARS,20,"Number of bars (i.e. bins to break [XMIN, XMAX] interval into)");
     CommandArgument_Bool_OrDefault_Doc(FRACTION,False,"If true, y axis is labeled as fractions of total number of binned values");
     CommandArgument_String_OrDefault_Doc(LABEL,"","Bargraph main label");
     CommandArgument_String_OrDefault_Doc(XLIM,"","Limit plot to [XMIN,XMAX] region (by default, data points "
					  "outside the region are discarded, see COUNT_OUTLIERS)");
     CommandArgument_Double_OrDefault_Doc(BARW,5,"Width of the bars");
     CommandArgument_Int_OrDefault_Doc(FONTSIZE,7,"Font size (labels)");
     CommandArgument_Bool_OrDefault_Doc(COUNT_OUTLIERS,False,"Count data points below XMIN or above XMAX into "
					"the first or the last bar, respectively");
     CommandArgument_String_OrDefault_Doc(OUT,"","Output graphics file name; defaults to DATA.COL.eps");
     CommandArgument_Bool_OrDefault_Doc(LOG,False,"Plot histogram of logarithms of (data value+1)");
     EndCommandArguments;

     istream *data;
     if ( DATA=="-" ) {
       data = &cin ;
     } else {
       if ( !IsRegularFile(DATA) ) FatalErr("Input file " << DATA << " does not exist");
       data = new ifstream(DATA.c_str());
     }
     if ( OUT.empty() ) {
       if ( DATA=="-" ) OUT = "hitogram.eps" ;
       else OUT = DATA + "." + ToString(COL) + ".eps";
     }

     float xmin=0, xmax=0; // will mot be used later unless  explicitly set by XLIM (below)
     if ( ! XLIM.empty() ) {
       if ( XLIM[0] != '[' || XLIM[XLIM.size()-1] != ']' || XLIM.size() < 5 ) {
	 cout << "Illegal XLIM argument format. Expected: [xmin,xmax]" << endl;
	 exit(1);
       }
       XLIM=XLIM.substr(1,XLIM.size()-2); // trim '[' and ']' from the string ends
       int pos = XLIM.Position(","); // find separator (comma)
       if ( pos == (-1) || pos==0 || static_cast<unsigned int>(pos) == XLIM.size()-1 ) {
	 cout << "Illegal XLIM argument format. Expected: [xmin,xmax]" << endl;
	 exit(1);
       }
       xmin = ToFloat( XLIM.c_str(), XLIM.c_str()+pos );
       xmax = ToFloat( XLIM.c_str()+pos+1, XLIM.c_str()+XLIM.size() );
       cout << "Plotting in ["<<xmin<<','<<xmax<<"] region" << endl;       
     }

     vec<float> values;

     String line;
     unsigned int line_no = 0;
     while(1) {
       getline( *data, line );
       if ( data->eof() ) break;

       line_no++;

       if ( ! line.Contains(PATTERN) ) continue;
       const char * ptr_e = line.c_str() - 1;
       const char * ptr_b ;
       int col = -1;
       do {
	 ptr_b = ptr_e + 1; // b points to the first character in column col+1
	 for ( ptr_e = ptr_b ; *ptr_e != '\0' && *ptr_e != ' ' && *ptr_e != '\t' ; ptr_e++);
	 // e now points to the first separator character (including eol) after column col+1
	 col++; // column col is spanned by [ptr_b , ptr_e)
       } while(col < COL && ptr_e != '\0');

       if ( col < COL ) FatalErr( "line " << line_no << " has only " << (col+1) << " columns");

       float val = ToFloat(ptr_b, ptr_e ) ;
       if ( isnan(val) ) {
	 cout << "line " << line_no << ": NOT A NUMBER in column " << COL << endl << endl;
	 exit(1);
       }

       values.push_back( ToFloat(ptr_b,ptr_e) );
     }

     if ( DATA != "-" ) (static_cast<ifstream *>(data))->close();

     cout << "Total number of binned values: " << values.size() << endl;

     PRINT2(Min(values),Max(values));

     if ( LOG ) {
       for ( size_t i = 0 ; i < values.size() ; i++ ) values[i] = log10f(values[i]+1);
     }

     sort(values.begin(),values.end(),greater<float>());
     float_histogram hist( values);
     vec<float> x_bars, y_bars;
 
     hist.SetBarNumber(NBARS);
     if ( ! XLIM.empty() ) hist.SetXLimits(xmin,xmax);
     hist.SetCountOutliers(COUNT_OUTLIERS);
     hist.ProduceHistogram(x_bars,y_bars);
     if ( FRACTION ) {
       for ( unsigned int i = 0 ; i < y_bars.size(); i++ ) y_bars[i]/=values.size();
     }

     ns_psplot::freetext label(LABEL,color(0,0,0),10);
     vec<ns_psplot::freetext> labels;
     labels.push_back(label);

     float bar_sep = 1.5;
     color bar_color = gray;
     float max_bar_height = 150;
     vec<float> bot_ticks;

     cout << "Plotting to " << OUT << "..." ;
     ns_psplot::bargraph graph = ns_psplot::BarPlot(y_bars, x_bars[0], x_bars[x_bars.size()-1],
						    static_cast<float>(BARW), bar_sep, bar_color,  
						    static_cast<short>(FONTSIZE),
						    max_bar_height, labels, bot_ticks, true, 2, true);

     ofstream f_fig (OUT.c_str());
     f_fig << graph << endl;
     cout << "done." << endl;
     
}
       
	 
       

