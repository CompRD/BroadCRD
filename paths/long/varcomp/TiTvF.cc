///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// TiTvF.  Compute Ti/Tv for SNPs that are in the Fosmid truth set within a 
// fixed distance from the nearest indel. Generates a table for different
// indel differences

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "STLExtensions.h"
#include "TokenizeString.h"
#include "util/TextTable.h"
#include "Histogram.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_IntSet_OrDefault_Doc(MIN_DIST_INDEL_SET, "{10,20,40,80,160,320,640,1000}", 
          "exclude variants that are not this close to an indel variant");
     CommandArgument_IntSet_OrDefault_Doc(DIST_INDEL_BIN, "{10,100,1000,10000,100000}",
          "upper-bound of bins of ti/tv log");
     CommandArgument_Bool_OrDefault_Doc(REMOVE_DBSNP, True, 
          "ignore SNPs in dbSNP");
     CommandArgument_String_OrDefault_Doc(FMINUS, "",
          "exclude these Fosmids from the analysis");
     CommandArgument_Bool_OrDefault_Doc(DISPLAY_HISTOGRAM, False,
          "display histogram for distances");
     CommandArgument_Double_OrDefault_Doc(LOG_BIN_SIZE, 0.5,
          "bin size");
     EndCommandArguments;

     vec<int> fminus;
     ParseIntSet( FMINUS, fminus );

     // Load dbSNP.

     vec< pair< pair<String,int>, pair<String,String> > > dbsnp;
     String line, junk, chr, ref, alt, gen;
     int pos;
     if ( IsRegularFile( "dbSNP.cache" ) )
          BinaryReader::readFile( "dbSNP.cache", &dbsnp );
     else
     {    fast_ifstream din( "/humgen/gsa-hpprojects/GATK/bundle/current/b37/"
               "dbsnp_138.b37.excluding_sites_after_129.vcf" );
          vec<String> dlines;
          dlines.reserve(14000000);
          while(1)
          {    getline( din, line );
               if ( din.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               if ( line[0] != 'X' && line[0] != 'Y' 
                    && ( line[0] < '0' || line[0] > '9' ) )
               {     continue;     }
               dlines.push_back(line);    }
          dbsnp.reserve( dlines.size( ) );
          for ( int z = 0; z < dlines.isize( ); z++ )
          {    const String& line = dlines[z];
               istrstream iline( line.c_str( ) );
               iline >> chr >> pos >> junk >> ref >> alt;
               if ( ref.size( ) > 1 || alt.size( ) > 1 ) continue;
               dbsnp.push( make_pair(chr,pos), make_pair(ref,alt) );    }
          ParallelSort(dbsnp);
          BinaryWriter::writeFile( "dbSNP.cache", dbsnp );    }

     // Load Fosmid variants.

     vec< pair< pair<String,int>, pair<String,String> > > fvars;
     fast_pipe_ifstream in( "cat /wga/scr4/jaffe/CompareVars/*/variants.all "
          "| grep Fosmid" );
     String loc;
     vec<String> ids;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          istringstream iline( line.c_str( ) );
          String id;
          iline >> id >> junk >> junk >> loc >> ref >> alt;
          id = id.After( "#" );
          if ( BinMember( fminus, id.Int( ) ) ) continue;
          chr = loc.Before( ":" );
          pos = loc.After( ":" ).Int( );
          fvars.push( make_pair(chr,pos), make_pair(ref,alt) );
          ids.push_back(id);    }
  
     // Keep distances for histogram
     vec<int> distance;

     //distance-decomposition of titv entries
     vec<int> transition_dist,transversion_dist;
     UniqueSort(DIST_INDEL_BIN);
     vec<int> ti_dist_bin(DIST_INDEL_BIN.size()+1,0);
     vec<int> tv_dist_bin(DIST_INDEL_BIN.size()+1,0);
     
     // Count transitions and transversions.
     Sort(MIN_DIST_INDEL_SET);
     int table_size = MIN_DIST_INDEL_SET.size();
     vec<int> transitions(table_size,0);
     vec<int> transversions(table_size,0);
     int transitions_outside = 0;
     int transversions_outside = 0;
     //     int transitions = 0, transversions = 0;
     for ( int i = 0; i < fvars.isize( ); i++ )
     {    if ( REMOVE_DBSNP && BinMember( dbsnp, fvars[i] ) ) continue;
          String R = fvars[i].second.first, A = fvars[i].second.second;
          if ( R.size( ) > 1 || A.size( ) > 1 ) continue;
          int dist_indel = 1000000000;
          for ( int j = i - 1; j >= 0; j-- )
          {    if ( ids[j] != ids[i] ) break;
               if ( fvars[j].second.first.size( ) == 1 
                    && fvars[j].second.second.size( ) == 1 ) continue;
               dist_indel = Min( 
                    dist_indel, fvars[i].first.second - fvars[j].first.second );    }
          for ( int j = i + 1; j < fvars.isize( ); j++ )
          {    if ( ids[j] != ids[i] ) break;
               if ( fvars[j].second.first.size( ) == 1 
                    && fvars[j].second.second.size( ) == 1 ) continue;
               dist_indel = Min( 
                    dist_indel, fvars[j].first.second - fvars[i].first.second );    }

	  distance.push_back(dist_indel);

	  int index = 0;
	  while ( (index < table_size) && (dist_indel > MIN_DIST_INDEL_SET[index]) )
	    index++;

	  size_t bin_idx=0;
	  for(;   bin_idx<DIST_INDEL_BIN.size()
	       && dist_indel >= DIST_INDEL_BIN[bin_idx]
	      ;++bin_idx)
	  {}

          Bool transition = False;
          if ( R == "A" && A == "G" ) transition = True;
          if ( R == "G" && A == "A" ) transition = True;
          if ( R == "C" && A == "T" ) transition = True;
          if ( R == "T" && A == "C" ) transition = True;
          if(transition){
              transition_dist.push_back(dist_indel);
              ++ti_dist_bin[bin_idx];
          }
          else{
              transversion_dist.push_back(dist_indel);
              ++tv_dist_bin[bin_idx];
          }

	  if (index == table_size)  // distance to indel is greater than table maximum
	    if (transition) transitions_outside++;
	    else transversions_outside++;   
	  else
	    if (transition) transitions[index]++;
	    else transversions[index]++;   
     }
     

     TextTable table;
     table << "MinDistToIndel" << Tab << "Transitions" << Tab << "Transversions"
	   << Tab << "Ratio" << EndRow;
     table << DoubleLine;
     
     int transitions_sum = 0;
     int transversions_sum = 0;
     for (int index = 0; index < table_size; index++) {
       transitions_sum += transitions[index];
       transversions_sum += transversions[index];
       table << MIN_DIST_INDEL_SET[index] << Tab << transitions_sum << Tab
	     << transversions_sum << Tab 
	     << double(transitions_sum)/double(transversions_sum) << EndRow;
     }
     // Add > Max indel depth value
     table << ">" << ToString(MIN_DIST_INDEL_SET.back()) << Tab << transitions_outside << Tab
	   << transversions_outside << Tab 
	   << double(transitions_outside)/double(transversions_outside) << EndRow;
     // Add global value
     table << "global" << Tab << transitions_outside + transitions_sum << Tab
	   << transversions_outside + transversions_sum << Tab 
	   << double(transitions_outside + transitions_sum)/double(transversions_outside + transversions_sum) << EndRow;

     table.Print( std::cout, 5, "r");
     
     if (DISPLAY_HISTOGRAM) {
       histogram<int> hist;
       hist.SetAllFromData(distance.begin(), distance.end(), 0, 100, 100);
       cout << hist;
     }
     if (DISPLAY_HISTOGRAM) {
          typedef double xval_t;
          vec<xval_t> ti_log,tv_log;
          ti_log.reserve(transition_dist.size());
          std::transform(transition_dist.begin(),transition_dist.end(),std::back_inserter(ti_log) ,[](int a){return log10(a);});
          tv_log.reserve(transversion_dist.size());
          std::transform(transversion_dist.begin(),transversion_dist.end(),std::back_inserter(tv_log) ,[](int a){return log10(a);});

          histogram<xval_t> hist_ti,hist_tv;
          hist_ti.SetAllFromData(ti_log.begin(), ti_log.end(), 0, LOG_BIN_SIZE, 10./LOG_BIN_SIZE+0.5);
          hist_tv.SetAllFromData(tv_log.begin(), tv_log.end(), 0, LOG_BIN_SIZE, 10./LOG_BIN_SIZE+0.5);
          std::cout << "log(d)\t#ti\t#tv\t#ti/#tv\n";
          std::transform(hist_ti.begin(),hist_ti.end()
                        ,hist_tv.begin(),std::ostream_iterator<String>(std::cout,"\n")
                        ,[](std::pair<xval_t,longlong>a,std::pair<xval_t,longlong>b){
                             ForceAssert(a.first==b.first);
                             return    ToString(a.first) + "\t"
                                     + ToString(a.second) + "\t" + ToString(b.second) + "\t"
                                     + ToString(double(a.second)/double(b.second));
                         }
                        );
     }
     {
         TextTable table;
         table <<DoubleLine;
         table <<"begin" << Tab << "end" << Tab << "Ti" << Tab << "Tv" << Tab << "Ti/Tv" << EndRow;
         table <<DoubleLine;
         for(size_t bb=0;bb<DIST_INDEL_BIN.size();++bb){
             table << (bb==0?0:DIST_INDEL_BIN[bb-1]) << Tab << DIST_INDEL_BIN[bb] << Tab << ti_dist_bin[bb] << Tab << tv_dist_bin[bb] << Tab << double(ti_dist_bin[bb])/double(tv_dist_bin[bb]) << EndRow;
         }
         table << (DIST_INDEL_BIN.size()==0?0:DIST_INDEL_BIN.back()) << Tab << "inf" << Tab
               << ti_dist_bin.back() << Tab << tv_dist_bin.back() << Tab << double(ti_dist_bin.back())/double(tv_dist_bin.back()) << EndRow;
         table.Print(std::cout);

     }
}
