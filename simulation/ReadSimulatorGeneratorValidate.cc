///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Sep 2012 
//

#include "MainTools.h"
#include "simulation/Framework.h"
#include "simulation/ReferenceTools.h"
#include "simulation/ReadSimulatorSimpleGenerator.h"
#include "Qualvector.h"
#include "simulation/ReadSimulatorSimpleCore.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "simulation/ReferenceIterator.h"
#include <memory>

namespace {

struct Log10Hist {
  vector<unsigned long> _hist;
  explicit Log10Hist(int size) : _hist(size,0) {};
  void addTo( unsigned int val ) {
    size_t bin = static_cast<size_t>( ( val == 0 ) ? 0 : log10(val) );
    if ( bin > _hist.size() - 1 )
      bin = _hist.size() - 1;
    _hist[bin]++;
  }
  void Print( std::ostream& str ) {
    unsigned bin = 1;
    for ( size_t i = 0; i < _hist.size(); ++i ) {
      cout << bin << "-" << bin*10 << "\t\t";
      bin*=10;
    }
    cout << endl;
    for ( size_t i = 0; i < _hist.size(); ++i ) {
      cout << _hist[i] << "\t\t";
    }
    cout << endl;
  }
};

void invert_locs(const String& LOCS, const String& PERFECT_GENOME, const BaseVecVec& reads_bvv )
{
  // read the locs
  RefLocusVec refv(0);
  BinaryReader refv_reader(LOCS);
  readBinary(&refv, refv_reader);

  // read the genome
  BaseVecVec genome_bvv( PERFECT_GENOME );

  cout << Tag() << "sampling the genome and comparing to " << reads_bvv.size() << " reads." << endl;

  size_t bad = 0;

  // compare each read base by base to the genome
  for ( size_t j = 0; j < reads_bvv.size(); ++j ) {

    size_t read_size = reads_bvv[j].size();

    ReferenceIterator riter( genome_bvv[ refv[j].mRefID], refv[j].mOffset + (refv[j].mRC ? read_size - 1 : 0), refv[j].mRC );

    bool bad = false;
    for ( size_t i = 0; i < read_size; ++i, ++riter ) {
      if ( reads_bvv[j][i] != *riter )
	bad++;
    }
  }

  cout << Tag() << "there should be NO bad reads: " << bad << " of " << reads_bvv.size() << endl;

}

void compare_locs_to_qltout(const String& LOCS, const BaseVecVec& reads_bvv, const String& QLTOUT, size_t fuzz = 0)
{
  // reads the locs
  RefLocusVec refv(0);
  BinaryReader refv_reader(LOCS);
  readBinary(&refv, refv_reader);

  // locs and reads should always be the same size
  ForceAssertEq( refv.size(), reads_bvv.size() );

  // read the qltout
  vec<vec<align_id_t> > aligns_index;
  vec<look_align> aligns;

  const String qltout = QLTOUT;
  ForceAssert(CountLookAligns(qltout));
  LoadLookAligns(qltout, aligns, aligns_index, reads_bvv.size());
  cout << "number of look aligns loaded=" << CountLookAligns(qltout) << endl;

  Log10Hist hist(4);

  // compare the qltout to the locs
  cout << Tag() << "validating read locs " << LOCS << " versus qltout lookaligns " << QLTOUT << endl;
  size_t nreflocs = refv.size();
  size_t good = 0U;
  size_t multiples = 0U;
  for (size_t i = 0; i < nreflocs; ++i) {

    if (i % 1000 == 0) {
      cout << ".";
      flush(cout);
    }

    const RefLocus& loc = refv[i];
    //    cout << refv[i].Str() << endl;
    size_t naligns = aligns_index[i].size();
    bool found = false;

    for (size_t j = 0;j < naligns;++j) {
      size_t ialign = aligns_index[i][j];
      const look_align& la = aligns[ialign];

      // the following is simplistic for the breadth of stuff possible in look_aligns and aligns, but should be good enough for perfect reads
      size_t start_on_target = static_cast<size_t>(la.StartOnTarget() );
      size_t absdiff = (start_on_target > loc.mOffset) ? start_on_target - loc.mOffset : loc.mOffset - start_on_target;
      if (loc.mRefID != static_cast<size_t>(la.target_id) || absdiff > fuzz|| loc.mRC != la.rc1) {

	cout << loc << endl;
	cout << j << ": " << la.query_id << "," << la.target_id << "," << la.a.StartOnTarget() << "," << (la.rc1 ? "rc" : "fw") << endl;

	if ( loc.mRefID == static_cast<size_t>(la.target_id) && loc.mRC == la.rc1 ) {
	  hist.addTo(absdiff);
	}
      } else {
	hist.addTo(absdiff);
	found = true;
      }
    }

    if (naligns > 1) {
      multiples++;
      cout << "multiple alignments found";
      if (found)
	cout << ", but one is correct.";
      else
	cout << ", and none are correct.";

      cout << endl;
    }

    if ( found ) ++good;


  } // end for reflocs...

  cout << endl << Tag() << "done!" << endl;

  cout << Tag() << static_cast<float>( good ) / nreflocs * 100 << "% are good within " << fuzz << " bases." << std::endl;
  cout << Tag() << static_cast<float>( multiples ) / nreflocs * 100 << "% had multiple alignments" << std::endl;

  hist.Print(cout);
}



void invert_errors( const BaseVecVec& errs_bvv, const ReadErrorVecVec& err_locs, const BaseVecVec& reads_bvv, bool verbose = false )
{
  ForceAssertEq(errs_bvv.size(), reads_bvv.size());
  ForceAssertEq(errs_bvv.size(), err_locs.size());

  for ( size_t j = 0; j < reads_bvv.size(); ++j ) {
    const BaseVec& errs = errs_bvv[j];
    const BaseVec& ref = reads_bvv[j];
    const ReadErrorVec& errloc = err_locs[j];

    // make sure that the error locations are sorted
#if 1
    int last = -1;
    for ( size_t i = 0; i < errloc.size(); ++i ) {
      ForceAssertGe(static_cast<int>(errloc[i].getLocation()), last);
      last = errloc[i].getLocation();
    }
#endif

    if ( verbose ) {
      cout << "original:" << endl;
      for ( size_t i = 0; i < std::min(80U,ref.size()); ++i )
	if ( i % 10 == 0 )
	  cout << i / 10;
	else
	  cout << " ";
      cout << endl;
      ref.Print(cout);
      cout << endl << "modified: " << endl;
      for ( size_t i = 0; i < std::min(80U,ref.size()); ++i )
	if ( i % 10 == 0 )
	  cout << i / 10;
	else
	  cout << " ";
      cout << endl;
      errs.Print(cout);
      cout << "locs: " << endl;
      for ( size_t i = 0; i < errloc.size(); ++i )
	cout << errloc[i] << endl;
    }
    BaseVec corrected(errs.size());

    size_t iref=0;
    size_t ierr=0;
    size_t iloc = 0;
    bool bad = false;
    while( ierr < errs.size() && iref < ref.size() ) {	// can't go past ref.size() without access to the actual reference genome...
      if ( errloc[iloc].getLocation() == iref ) {
	switch ( errloc[iloc].getType() ) {
	case ReadError::SUBSTITUTION:
	  if ( errs[ierr] != errloc[iloc].getReadBase() ) {
	    cout << "bad subs" << endl;
	    bad = true;
	  }
	  iref++;
	  ierr++;
	  iloc++;
	  break;
	case ReadError::INSERTION:
	  if ( errs[ierr] != errloc[iloc].getReadBase() ) {
	    cout << "bad ins" << endl;
	    bad=true;
	  }
	  ierr++;
	  iloc++;
	  break;
	case ReadError::DELETION:
	  if ( ref[iref] != errloc[iloc].getRefBase() ) {
	    cout << "bad del" << endl;
	    bad = true;
	  }
	  iref++;
	  iloc++;
	  break;
	}


      } else {
	if ( ref[iref] != errs[ierr] ) {
	  cout << "bad base" << endl;
	  bad =true;
	}
	iref++;
	ierr++;
      }

    }

    if ( bad )
      cout << "error reconstructing read0=" << j << endl;
  }
}

}

int main(int argc, char* argv[])
{
  RunTime();
  BeginCommandArguments;
  CommandArgument_String_Doc(PERFECT_READS, "Perfect reads input filename. If specified, use to perturb.");
  CommandArgument_String_Doc(LOCS, "Locations");
  CommandArgument_String_OrDefault_Doc(PERFECT_QLTOUT, "", ".qltout");
  CommandArgument_String_OrDefault_Doc(PERFECT_GENOME, "", "genome that created the reads");
  CommandArgument_String_Doc(ERROR_LOCS, "Error locations");
  CommandArgument_String_OrDefault_Doc(ERROR_QLTOUT, "", "error.qltout");
  CommandArgument_String_Doc(ERROR_READS, "Error reads input filename.");
  CommandArgument_UnsignedInt_OrDefault_Doc(FLUFF,10U,"how many bases slip to allow for the error alignments.")
  EndCommandArguments;


  // read the reads
  BaseVecVec reads_bvv;
  reads_bvv.ReadAll(PERFECT_READS);

  // compare the locs to alignments from the perfect reads
  if ( PERFECT_QLTOUT != "" )
    compare_locs_to_qltout(LOCS, reads_bvv, PERFECT_QLTOUT);

  // invert the read locs
  if ( PERFECT_GENOME != "" )
    invert_locs(LOCS, PERFECT_GENOME, reads_bvv);

  // compare the locs to alignments from the error reads
  if ( ERROR_QLTOUT != "" )
    compare_locs_to_qltout(LOCS, reads_bvv, ERROR_QLTOUT, 10U);

  // compare the recorded errors to the errors recovered during alignment
//  compare_error_locs_to_error_qltout( ERROR_LOCS, reads_bvv, ERROR_QLTOUT );

  // invert the errors and compare the reads -- should be exact
  BaseVecVec errs_bvv;
  errs_bvv.ReadAll(ERROR_READS);

  ReadErrorVecVec err_locs;
  err_locs.ReadAll(ERROR_LOCS);

  cout << "inverting the errors and comparing error reads to perfect reads" << endl;

  invert_errors(errs_bvv, err_locs, reads_bvv );


  return 0;

}
