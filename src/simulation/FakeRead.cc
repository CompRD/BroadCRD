///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "simulation/FakeRead.h"
#include <fstream>

void CFakeRead::Set(const basevector & bases,
		    const qualvector & quals,
		    int from, int to,
		    bool r,
		    int size,
		    const String & lib,
		    const String & plate,
		    const String & name,
		    const String & tmp,
		    const String & well,
		    bool bPaired)
{
  
  m_bPaired = bPaired;

  if (r)
    m_fr = "R";
  else
    m_fr = "F";
  
  m_lib = lib;
  m_plate = plate;
  m_name = name;
  m_template = tmp;
  
  if (from < 0)
    from = 0;
  if (to > (int)bases.size())
    to = (int)bases.size();
  
  m_bases.SetToSubOf(bases, from, to - from);
  m_quals.SetToSubOf(quals, from, to - from);
  
  //for (int i=0; i<(int)m_bases.size(); i++) {
    //m_quals[i] = 30;
    //int n = (big_random() % 1000000);
    //double w = (double)n / 1000000;
    //if (w < 0.004) {
      //cout << "Making read error..." << endl;
      //char base = m_bases[i];
      //base += (1 + big_random() % 3);
      //base = (base % 4);
      //m_bases.Set(i, base);
      //m_quals[i] = 20;
    // }
  //}

  if (r) {
    m_bases.ReverseComplement();
    ReverseThis(m_quals);
  }
  
  m_size = size;
  m_well = well;
}

void CFakeRead::Print( ostream& xmls, ostream& bases, ostream& quals ) const {
  xmls << "  <trace>\n";
  xmls << "        <trace_name>" << m_name << "</trace_name>\n";
  xmls << "        <plate_id>" << m_plate << "</plate_id>\n";

  static int g_well;
  if ( m_well.empty() ) {
    xmls << "        <well_id>A" << g_well++ << "</well_id>\n";
  } else {
    xmls << "        <well_id>" << m_well << "</well_id>\n";
  }

  xmls << "        <template_id>" << m_template << "</template_id>\n";
  xmls << "        <trace_end>" << m_fr << "</trace_end>\n";
  xmls << "        <clip_vector_left>1</clip_vector_left>\n";
  xmls << "        <clip_vector_right>" << m_bases.size() << "</clip_vector_right>\n";
  //xmls << "        <clip_quality_right>" << m_bases.size() << "</clip_quality_right>\n";
  xmls << "        <library_id>" << m_lib << "</library_id>\n";

  if ( m_bPaired ) {
    xmls << "        <insert_size>" << m_size << "</insert_size>\n";
    xmls << "        <insert_stdev>" << m_size/15 << "</insert_stdev>\n";
  } else {
    xmls << "        <type>unpaired_production</type>\n";
  }

  xmls << "  </trace>\n";

  m_bases.Print(bases,m_name);
  ::Print(quals,m_quals,m_name);
}

void PrintFakeReads( vec<CFakeRead> const& fakeReads, String const& xmlFile,
                            String const& baseFile, String const& qualFile )
{
    std::ofstream xmls(xmlFile.c_str());
    std::ofstream bases(baseFile.c_str());
    std::ofstream quals(qualFile.c_str());

    xmls << "<?xml version=\"1.0\"?>\n";
    xmls << "<trace_volume>\n";
    xmls << "<volume_name>optimap-redoids</volume_name>\n";
    xmls << "<volume_date>2005-06-05</volume_date>\n";
    xmls << "<volume_version>optimap.1</volume_version>\n";

    typedef vec<CFakeRead>::const_iterator Itr;
    for ( Itr itr(fakeReads.begin()), end(fakeReads.end()); itr != end; ++end )
        itr->Print(xmls,bases,quals);

    xmls << "</trace_volume>\n";
    xmls.close();
    ForceAssert(xmls);

    bases.close();
    ForceAssert(bases);
    quals.close();
    ForceAssert(quals);
}
