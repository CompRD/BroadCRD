///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* AssemblyBadCoverage: find positions in the reference that are not covered by any contig in the assembly */

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "paths/ReadLoc.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "CoverageAnalyzer.h"
#include "FetchReads.h"
#include <omp.h>
#include <stdexcept>
#include <sstream>

class AssemblyBadCoverage {
public:
  AssemblyBadCoverage( const String& qlt_file, const String& ref_file )  {
    LoadAligns( qlt_file );	// initializes aligns_
    LoadReftigs( ref_file );	// initializes refnames_, reftigs_
  }

  void FindUncovered() {
    // convert lookaligns to sequence intervals for CoverageAnalyzer
    cout << Date() << ": converting " << aligns_.size() << " lookaligns to sequence intervals..." << endl;
    vector<seq_interval> seqv;
    seqv.resize(aligns_.size());

    for ( size_t i = 0; i < aligns_.size(); ++i ) {
      seqv[i].Set(i,aligns_[i].TargetId(),aligns_[i].pos2(),aligns_[i].Pos2());
    }


    // find uncovered portions of the genome
    cout << Date() << ": Coverage analysis..." << endl;
    CoverageAnalyzer analysis( &cerr );
    analysis.CreateCoverages(seqv);
    analysis.GetCoveragesAtMost(0, uncoveredv_);
    DPRINT(uncoveredv_.size());
  }

  void WriteUncoveredToRangeFile( const String& out_name ) {
    WriteSeqIntervalToRangeFile(uncoveredv_, out_name);
  }


  void WriteSeqIntervalToRangeFile(const vec<seq_interval>& uncoveredv, const String& out_name ) {
    // write output in a format suitable for BadCoverage TEST_INTERVALS=...
    ofstream outfile( out_name.c_str() );

    if ( !outfile.is_open() ) {
      ostringstream serr;
      serr << "error opening file " << out_name << " for writing.";
      throw runtime_error( serr.str() );
    }
    cout << Date() << ": writing output to " << out_name << "...";
    for ( size_t i = 0; i < uncoveredv.size(); ++i ) {
      const seq_interval& s = uncoveredv[i];
      outfile << refnames_[s.SeqId()];
      outfile << ":" << s.Begin() << "-" << s.End() << endl;
    }
    cout << "done!" << endl;

    outfile.flush();
    if ( ! outfile.good() )
      throw runtime_error( "WARNING: an I/O error has occurred writing the output." );

    outfile.close();
  }

  
  void FindErrors() {
    cout << endl << Date() << ": caculating high error portions of assemblies" << endl;

    if ( querytigs_.size() == 0 )
      throw logic_error("FindErrors() called before query contigs were loaded...");

    // fat, wasteful, but does the job for now
    tigerr_.resize( reftigs_.size() );
    for ( size_t rtig = 0; rtig < reftigs_.size(); rtig++ )
      tigerr_[rtig].resize( reftigs_[rtig].size(), 0 );

    // for each lookalign, find the reference and query sequences and catalog errors per-reftig-base
    int naligns_done = 0;
#pragma omp parallel for
    for ( size_t i = 0; i < aligns_.size(); ++i ) {
#pragma omp critical
      {
	naligns_done++;
	if ( aligns_.size( ) < 100 || naligns_done % ( aligns_.size()/100 ) == 0 )
	  Dot(cout, 100.0 * naligns_done/(double)aligns_.size() );
      }

      look_align& la = aligns_[i];


      const BaseVec& ref = reftigs_[la.TargetId()];	// alias the reference
      BaseVec query = querytigs_[la.QueryId()];		// copy the query as we may RC() it

      //      la.PrintVisual( cout, query, ref	);

      if ( la.IsQueryRC())
	query.ReverseComplement();

      //      query.PrintBases(cout, la.pos1(), min(50,la.GaplessLength()));
      //      ref.PrintBases(cout, la.pos2(), min(50,la.GaplessLength()));

#if NEIL_DEBUG
      bool intrusive = ( la.TargetId() == 11 && la.StartOnTarget() > 2040000 && la.EndOnTarget() < 2045000 );
#endif
      for ( int i = la.StartOnQuery(); i < la.EndOnQuery(); ++i ) {
	const int& ri = la.a.PosOn2(i);
#if NEIL_DEBUG
	if ( intrusive && ri >= 0) {
	  cout << refnames[la.TargetId()] << ": " << ri << " " << Base::val2Char(ref.at(ri)) << Base::val2Char(query.at(i)) << endl;
	}
#endif
	if ( ri >=  0 && ref.at(ri) != query.at(i) ) {
#pragma omp critical
	  {
	    tigerr_[la.TargetId()][ri]=1;		// was ++ to account for errors on overlapping contigs, but that causes accounting woe
	  }
	}
      }
    } // for lookaligns
    cout << endl;

  }


  void WriteErrorsToIGVFile( int window_radius, double error_thresh, const String& bad_file ) {
    DPRINT2( window_radius, error_thresh );

    ofstream badfile( bad_file.c_str() );
    if ( !badfile.is_open() ) {
      ostringstream serr;
      serr << "error opening file " << bad_file << " for writing assembly errors.";
      throw runtime_error( serr.str() );
    }


    // slide a window across the reftigs and emit ranges with suprathreshold errors
    cout << Date() << " assembly contigs ";
    unsigned short errThreshCount = static_cast<unsigned short>(error_thresh * ( 2*window_radius+1 ));

    DPRINT(errThreshCount);
    ostringstream err_name;
    err_name << "errors-rad" << window_radius << "-thr" << error_thresh;
    badfile << "chr\tstart\tend\tdummy\t" << err_name.str() << endl;
    for ( unsigned int ti = 0; ti < tigerr_.size(); ++ti ) {
      Dot(cout, ti/tigerr_.size());
      int start = -1*window_radius;
      int end = reftigs_[ti].size();
      int count = 0;
      int total_count = 0;
      int err_range = -1;
      int max_err = 0;

      for ( int i = start; i < end; ++i  ) {
	if ( i + window_radius < end ) {
	  count += tigerr_[ti][i+window_radius];
	  total_count++;
	}
	if ( i - window_radius >= 0 ) {
	  count -= tigerr_[ti][i-window_radius];
	  total_count--;
	}


	//	if ( static_cast<double>(count)/total_count > ERROR_THRESH && err_range < 0 )	// starting a range of errors
	if ( count > errThreshCount && err_range < 0 )
	  err_range = i;
	if ( count > errThreshCount && count > max_err )
	  max_err = count;
	if ( count <= errThreshCount && err_range >= 0) {	// ending a range of errors
	  badfile << refnames_[ti] << "\t" << err_range << "\t" << i << "\t" << "errors" << "\t" << static_cast<float>(max_err) / total_count << endl;
//	  badfile << "\t# " << refnames_[ti] << ":" << err_range << "-" << i << endl;
	  err_range = -1;		// mark us as being in an error-free region
	  max_err = 0;
	}
      }
    }

#if NEIL_DEBUG
    for ( unsigned int bi = 2044725; bi < 2045026; bi++ ) {
      cout << bi << ":" << tigerr_[11][bi] << endl;
    }
#endif


    badfile.close();
    if ( !badfile.good() ) {
      ostringstream serr;
      serr << "WARNING: an I/O error occurred while writing the badfile: " << bad_file << endl;
      throw runtime_error( serr.str() );
    }

  }


  void WriteGCToIGVFile( int window_radius, const String& gc_file ) {
      DPRINT(window_radius);

      ofstream gcfile( gc_file.c_str() );
      if ( !gcfile.is_open() ) {
        ostringstream serr;
        serr << "error opening file " << gc_file << " for writing assembly errors.";
        throw runtime_error( serr.str() );
      }


      // slide a window across the reftigs and emit single-base ranges with moving average GC content
      gcfile << "#track color=000,200,000" << endl;
      gcfile << "chr\tstart\tend\tdummy\tGCwin" << window_radius << endl;
      for ( unsigned int ti = 0; ti < reftigs_.size(); ++ti ) {
        Dot(cout, ti/reftigs_.size());
        int start = -1*window_radius;
        int end = reftigs_[ti].size();
        int count = 0;
        int total_count = 0;
        int q; 	// query point
        for ( int i = start; i < end; ++i  ) {
          q = i+window_radius;
          if ( q < end ) {
            if ( reftigs_[ti][q] == Base::G.val() || reftigs_[ti][q] == Base::C.val() )
              count++;
            total_count++;
          }
          q = i-window_radius;
          if ( q >= 0 ) {
            if ( reftigs_[ti][q] == Base::G.val() || reftigs_[ti][q] == Base::C.val() )
              count--;
            total_count--;
          }

          if ( i >= 0 && i < end ) {
            gcfile << refnames_[ti] << "\t" << i << "\t" << i << "\t" << "GC" << "\t" << static_cast<double>(count)/total_count << endl;
//            gcfile << "\t# " << refnames_[ti] << ":" << i << "-" << i << endl;
          }
        }
      }

  #if NEIL_DEBUG
      for ( unsigned int bi = 2044725; bi < 2045026; bi++ ) {
        cout << bi << ":" << tigerr_[11][bi] << endl;
      }
  #endif


      gcfile.close();
      if ( !gcfile.good() ) {
        ostringstream serr;
        serr << "WARNING: an I/O error occurred while writing the badfile: " << gc_file << endl;
        throw runtime_error( serr.str() );
      }

    }


  void LoadQueries( const String& query_file ) {
    // read assembly 'tigs in
    cout << endl << Date() << ": loading assembly contigs" << endl;
    ForceAssert( IsRegularFile(query_file) );
    FetchReads(querytigs_, querynames_, query_file);
    DPRINT(querytigs_.size());
  }

private:
  void LoadAligns( const String& qlt_file ) {
    // load alignments...
    cout << Date() << ": loading alignments from " << qlt_file << endl;
    ForceAssert( IsRegularFile(qlt_file));
    LoadLookAligns(qlt_file,aligns_);
    cout << endl;
  }

  void LoadReftigs( const String& ref_file ) {
    ForceAssert( IsRegularFile(ref_file) );
    FetchReads(reftigs_, refnames_, ref_file);

    DPRINT(reftigs_.size());
    cout << endl;
  }



  AssemblyBadCoverage() {};

  vec<look_align> aligns_;
  vecbasevector reftigs_;
  vecString refnames_;
  vec<seq_interval> uncoveredv_;
  vec<vec<short> > tigerr_;
  vecbasevector querytigs_;
  vecString querynames_;
};

int main(int argc, char *argv[]){
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String_Doc(QLT_FILE, "parsable output of QueryLookupTable relating assembly contigs to the reference");
  CommandArgument_String_Doc(REF_FILE, "reference genome file in FASTA format");
  CommandArgument_String_Doc(OUT_FILE, "file to which uncovered ranges are written.");
  CommandArgument_String_OrDefault_Doc(BAD_FILE, "", "file to which high-error assembly ranges are written");
  CommandArgument_Int_OrDefault_Doc( NUM_THREADS, 16, "number of threads to be used" );
  CommandArgument_String_OrDefault_Doc(QUERY_FILE, "", "assembly contigs FASTA that QLT_FILE relates to the reference, if BAD_FILE is specified");
  CommandArgument_Int_OrDefault_Doc(WINDOW_RADIUS, 50, "radius of a sliding window in which to compute errors, if BAD_FILE is specified");
  CommandArgument_Double_OrDefault_Doc(ERROR_THRESH, 0.09, "fraction of errors within a sliding window for it to be considered part of a range of 'bad' windows");
  CommandArgument_String_OrDefault_Doc(GC_FILE, "", "file to which GC-content track is written, if specified");
  EndCommandArguments;
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads(NUM_THREADS);


  AssemblyBadCoverage abc(QLT_FILE, REF_FILE);
  abc.FindUncovered();
  abc.WriteUncoveredToRangeFile(OUT_FILE);

  // if requested, find high-error portions of assemblies
  if ( BAD_FILE != "" ) {
    abc.LoadQueries(QUERY_FILE);
    abc.FindErrors();
    abc.WriteErrorsToIGVFile(WINDOW_RADIUS, ERROR_THRESH, BAD_FILE);
    if ( GC_FILE != "") abc.WriteGCToIGVFile(WINDOW_RADIUS, GC_FILE);
  }

  cout << endl << Date() << ": Done with AssemblyBadCoverage!" << endl;
  return 0;
}    
    
