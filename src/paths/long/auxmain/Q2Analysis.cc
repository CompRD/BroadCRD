///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Performs analysis on variants found by DISCOVAR and present in the Fosmid,
// but not found by GATK-250
//
// The input file should be prepared from the VCF file using:
//   echo "#fosmid loc-local loc-global ref alt" > $filtered
//   grep -v -e "GATK-250" my.vcf  | grep -e "Fosmid" | grep -e "DISCOVAR" | awk '{print $1, $3, $4, $5, $6}' >> filtered_for_q2
//
// Description of results:
//
// Events: fosmid variants called by Discovar but not GATK.
//
// alignable - two or more reads in each direction span a window - extended 10 bases on each side of variant.
// non-Q2 alignable - as above, but only if there are no Q2 bases in the window.
// de nevo - less than two reads in each direction span the window


// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable


#include "MainTools.h"
#include "Basevector.h"
#include "FastIfstream.h"
#include <iostream>
#include "String.h"
#include "lookup/LookAlign.h"


void tokenize( String buf, vec<String>& fields, char separator = ' ' ) {
  fields.clear();

  String::iterator last = buf.begin();
  for ( String::iterator it=buf.begin(); it!=buf.end(); ++it)
    if( *it == separator) {
      fields.push_back(String(last, it));
      last = it + 1;
    }
  if (last < buf.end())
    fields.push_back(String(last, buf.end()));
}


class Analyser {

private:
  int m_last_fosmid_id, m_fosmid_id, m_fosmid_pos;
  String m_chr_name;
  int m_chr_pos1;   // one based co-ordinate
  String m_ref, m_alt;

  int m_ref_start0, m_ref_stop0;   // zero based co-ordinates

  String m_bam_filename;
  String m_ref_filename;

  int m_window = 10;

  int m_log_level = 0;

  vecbasevector ref;

  int m_event_count = 0;
  int alignment_callable = 0;  // callable by alignment 
  int nonq2_alignment_callable = 0; // callable by aligment without using q2 bases
  int denovo_callable = 0; // not callable by alignment

  uint64_t m_total_bases;
  uint64_t m_q2_bases;

  bool parse( String line ) {
    vec<String> fields;
    tokenize(line, fields, ' ');
    if (fields.empty() || fields[0] == "#fosmid" || fields[0][0] != '#')
      return false;
    else {
      int index = 0;
      m_fosmid_id = fields[index++].SafeAfter("#").Int();
      m_fosmid_pos = fields[index++].Int();
      m_chr_name = fields[index].SafeBefore(":");
      m_chr_pos1 = fields[index++].SafeAfter(":").Int();
      m_ref = fields[index++];
      m_alt = fields[index++];

      return true;
    }
  }

  void prepare () {
      m_ref_start0 = m_chr_pos1 - 1 - m_window;
      m_ref_stop0 = m_chr_pos1 - 1 + m_ref.size() + m_window;
  }

  void import_reference() {
    ref.ReadAll( m_ref_filename );
  }

  void import_fosmid_data () {
    cout << "Importing data for Fosmid: " << m_fosmid_id << endl;
  }

  void extract_reads_from_bam () {
    String range = m_chr_name + ":" + ToString(m_ref_start0 + 1) + "-" + ToString(m_ref_stop0 + 1);
    String out_head = "testy";
    String samtools_cmd = "samtools view -h " + m_bam_filename + " " + range;
    String SAM2CRDDump_cmd = "SAM2CRDDump OUT_HEAD=" + out_head +
      " USE_OQ=True NH=True LOG_TO_CERR=False REQUIRE_LIBINFO=False WRITE_ALIGNS=True" +
      " WRITE_NAMES=False WRITE_PAIRS=False MAPPED_PAIRS_ONLY=False";

    cout << "Extracting reads from region " << range << " from bam." << endl;

    int status = System( samtools_cmd + " | " + SAM2CRDDump_cmd );
    if ( status != 0 ) 
      FatalErr("Unable to extract reads from the bam file");

    uint64_t read_count = MastervecFileObjectCount( out_head + ".fastb" );
    cout << "Found " << read_count << " reads." << endl;
  }

  void compute_stats() {
    vec<look_align> aligns; 
    String out_head = "testy";
    
    m_event_count++;

    int nonq2_read_count_fw = 0;
    int read_count_fw = 0;
    int nonq2_read_count_rc = 0;
    int read_count_rc = 0;

    int chr_id = ( m_chr_name == "X" ? 22 : (m_chr_name == "Y" ? 23 : m_chr_name.Int() - 1 ));

    // Import reads, quals and aligns
    vecbasevector bases( out_head + ".fastb" );
    vecqualvector quals( out_head + ".qualb" );
    LoadLookAligns( out_head + ".qltout", aligns );

    // Examine each alignment
    for ( size_t index = 0; index < aligns.size(); ++index) {
      const look_align& lalign = aligns[index];
      int soq = lalign.StartOnQuery();
      int eoq = lalign.EndOnQuery();
      int sot = lalign.StartOnTarget();
      int eot = lalign.EndOnTarget();
      bool rc = lalign.Rc1();
      int query_id = lalign.QueryId();

      if (m_log_level > 1) {
	cout << query_id << " ";
	cout << soq << " ";
	cout << eoq << " ";
	cout << (eoq - soq) << " ";
	cout << sot << " ";
	cout << eot << " ";
	cout << (eot - sot) << " ";
	cout << ToStringBool(rc) << " ";
	cout << m_ref << " ";
	cout << m_alt << " ";
	cout << endl;
      }

      if (sot > m_ref_start0 || eot < m_ref_stop0) {
	if (m_log_level > 1) cout << "Skipping: window not fully covered" << endl;
	continue;
      }
      

      align a = lalign.a;

      // Find read start and stop bases corresponding to the reference window
      int ref_start = m_ref_start0;
      int ref_stop = m_ref_stop0;
      int query_start = -1;
      int query_stop = -1;
      while (query_start == -1 && ref_start < ref_stop)
	query_start = a.PosOn1(ref_start++);
      while (query_stop == -1 && ref_stop > ref_start)
	query_stop = a.PosOn1(ref_stop--);
      if (ref_stop == ref_start) {
	if (m_log_level > 1) cout << "Skipping: No aligned bases in the window." << endl;
	continue;
      }

      // grab the read
      basevector read(bases[query_id]);
      qualvector qual(quals[query_id]);
      if (rc) {
	read.ReverseComplement();
	std::reverse(qual.begin(), qual.end());
      }

      // extract the region around the window
      basevector window_query(read.begin() + query_start,read.begin() + query_stop);
      basevector window_target(ref[chr_id].begin() + m_ref_start0,ref[chr_id].begin() + m_ref_stop0);
      qualvector window_quals(qual.begin() + query_start,qual.begin() + query_stop);

      m_total_bases+=window_quals.size();
      m_q2_bases+=std::count_if(window_quals.begin(),window_quals.end(),[](typename qualvector::value_type in){return in==2;});

      if (m_log_level > 2) {
	cout << window_query.ToString() << endl;
	cout << window_target.ToString() << endl;
	Print(cout, window_quals, "name");
	//      lalign.PrintVisual( cout, bases[index], quals[index], ref[chr_id], 0);
	cout << "=============================================================" << endl;
      }

      if (rc)
	read_count_rc++;
      else
	read_count_fw++;
      if (std::find(window_quals.begin(), window_quals.end(), 2) == window_quals.end())
	if (rc)
	  nonq2_read_count_rc++;
	else
	  nonq2_read_count_fw++;

    }
    if (m_log_level > 0) {
      cout << read_count_fw << " " << read_count_rc << endl;
      cout << nonq2_read_count_fw << " " << nonq2_read_count_rc << endl;
    }
    
    if (read_count_fw >= 2 && read_count_rc >=2) {
      alignment_callable++;
      if (nonq2_read_count_fw >= 2 && nonq2_read_count_rc >=2)
	nonq2_alignment_callable++;
    } else
      denovo_callable++;
  }
      
  
public:
  
  Analyser(const String& bam_name, const String& ref_name, int log_level) :
    m_bam_filename(bam_name), m_ref_filename(ref_name) , m_log_level(log_level) ,
    m_total_bases(0), m_q2_bases(0){
    import_reference();
  };
  
  void analyse (String line, int target_fosmid) {
    if (parse (line) ) {
      if (target_fosmid != 0 && m_fosmid_id != target_fosmid)
	return;
      if (m_last_fosmid_id != m_fosmid_id) 
	import_fosmid_data();
      prepare();
      extract_reads_from_bam();
      compute_stats();
      m_last_fosmid_id = m_fosmid_id;
    }
  }
  
  void output() {
    cout << "Event count: " << m_event_count << endl;
    cout << "Non-Q2 alignment callable: " << nonq2_alignment_callable << endl;
    cout << "Alignment callable: " << alignment_callable << endl;
    cout << "De novo callable: " << denovo_callable << endl;
    cout << "Fraction of q2 base calls: " << double(m_q2_bases)/double(m_total_bases) << endl;
  }
  
};




int main(int argc, char *argv[])
{
  RunTime( );
  // essential arguments
  BeginCommandArguments;
  // logging options
  CommandArgument_String(INPUT);
  CommandArgument_Int_OrDefault(FOSMID,0);
  CommandArgument_Int_OrDefault(LOG_LEVEL,0);
  EndCommandArguments;
  // Check command line arguments
  
  Analyser analyser("/wga/scr4/NA12878_calls/H01UJADXX/realign-bwa-mem/reverted.12.aligned.wholegenome.sorted.bam",
		    "/wga/dev/references/Homo_sapiens/genome.fastb",
		    LOG_LEVEL);

  String line;
  ifstream myfile (INPUT);
  if (myfile.is_open()) {
    while ( getline (myfile,line) )   {
      analyser.analyse(line, FOSMID);
    }
    myfile.close();
  } else 
    cout << "Unable to open file"; 

  analyser.output();

  cout << Date() << ": Done!" << endl;    
}
