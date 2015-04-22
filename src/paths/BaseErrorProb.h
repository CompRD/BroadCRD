///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BASEERRORPROB
#define BASEERRORPROB

#include "MainTools.h"
#include "system/System.h"
#include "math/Combinatorics.h"
#include "Qualvector.h"

/**
   \file

   Code related to error profiles in reads: for each read position,
   what is the probability of error at that position?

   \ingroup grp_quals
*/

class BaseErrorProb;   // Forward Declaration

/// Table of Base Error Probabilities: Represents a list of possible
/// error position combinations together with the probability of each
/// combination, sorted from most probable to least probable.
/// Produced by BaseErrorProbProfile::getErrorTable().
/// \ingroup grp_quals
typedef vec<BaseErrorProb> BaseErrorTable;


/*
    Class: BaseErrorProb

    Holds the probability of a particular combination of base errors
    occuring.  The "combination" here refers only to error positions
    in the read, not to the read content; so, an object of this class
    might represent the probability that a read has errors at
    positions 10 and 15.
    
*/
class BaseErrorProb {

public:

  /// Positions of bases in the read: this class represents the
  /// probability that all positions listed here have errors.
  vec<int> base_pos;
  /// The length of the read, to which the base positions in \p
  /// base_pos refer.
  unsigned int read_length;
  /// The probability that all positions in \p base_pos have errors.
  double error_prob;

  BaseErrorProb() {
    read_length = 0;
    error_prob = 0.0;
  }

  /// Constructor using vec<bool> base error map
  BaseErrorProb(const vec<bool>& base_map, double error_probability);

  /// Constructor using string of 0's (no error) and 1's (error) as error map
  BaseErrorProb(const string& base_map, double error_probability);

  /// Constructor using a vector of error base positions and read length
  BaseErrorProb(const vec<int>& base_positions, int length, double error_probability)
    : base_pos(base_positions) {
    error_prob = error_probability;
    read_length = length;
  }

  int GetReadLength() const { return read_length; }

  /// Returns error map as string of 0's and 1's, representing the
  /// combination of read positions: this object represents the
  /// probability that all the positions indicated by 1's are
  /// simultaneously wrong in a read.
  string GetBaseMap() const;

  friend bool operator> (BaseErrorProb a, BaseErrorProb b) {
    return a.error_prob > b.error_prob;
  }

  friend bool operator< (BaseErrorProb a, BaseErrorProb b) {
    return a.error_prob < b.error_prob;
  }

  friend ostream& operator<<( ostream& o, const BaseErrorProb& p );
  friend istream& operator>>( istream& i,  BaseErrorProb& p );
  friend istream& operator>>( istream& i,  BaseErrorTable& t );
  
};


/*
  Class: BaseErrorProbProbile
  
  Class to contain the per base error probability for a read and to
  generate BaseErrorTable - a table of possible base error combinations
  sorted by decreasing probability.
*/
class BaseErrorProbProfile {

private:

  /// For each position in a read, the probability of error at that position. 
  vec<double> probs;

  /// Length of reads, for which this object represents the probability of
  /// errors at various positions.  Should equal probs.size() .
  unsigned int table_size;

public:

  BaseErrorProbProfile() {
    table_size = 0;
  }

  /// Constructor using a vec of base error probability, ordered by base position
  BaseErrorProbProfile(vec<double> probabilities) 
    : probs(probabilities) {
    table_size = probs.size();
  }

  /// Constructor to create a flat error profile for reads of up to num_bases in length
  BaseErrorProbProfile(int num_bases) 
    : probs(num_bases, 0.01) {
    table_size = probs.size();
  }

  /// Constructor to create an error profile from the quality scores of an individual read.
  /// In this case the probability of error at each read position is based not only on
  /// the position, but on what actually happened during that particular read.
  /// \sa base_prob
  BaseErrorProbProfile(const qualvector& quals);

  /// Constructor that imports a list of error probabilites (given in %) from a file
  BaseErrorProbProfile(string filename); 

  /// Returns a BaseErrorTable for reads of length \p num_bases, sorted by probability
  /// of the error (from most probable to least probable).  Each entry in the returned
  /// BaseErrorTable represents a particular combination of errors, that is, errors
  /// occurring simultaneously at a particular set of positions in the read.
  /// \param num_bases length of reads considered.  if zero, \p table_size (the read length
  ///        given when this object was constructed) is used.
  /// \param min_errors we consider combinations of at least this many errors
  /// \param max_errors we consider combinations of at most this many errors
  /// \param max_entries if non-zero, limits the number of entries in the returned
  /// \param peek_ahead if true then calculate probability of most likely max_errors+1
  ///        combination and append to table and trim table to this probability
  BaseErrorTable getErrorTable(unsigned int num_bases, unsigned int min_errors, 
			       unsigned int max_errors, const unsigned int max_entries,
			       const bool peek_ahead = false) const;


  /// \see getErrorTable(unsigned int,unsigned int,unsigned int,const unsigned int,const bool)const
  BaseErrorTable getErrorTable(const unsigned int num_bases,
			       const unsigned int max_errors = 2,
			       const unsigned int max_entries = 0,
			       const bool peek_ahead = false) const {
    return getErrorTable(num_bases, 1, max_errors, max_entries, peek_ahead);
  }

  

  unsigned int getProbProfileSize() const {
    return table_size;
  }


private:

  double prob(const int num_bases, const int num_errors,
	      BaseErrorTable& error_table,
	      const vec<double>& per_base_err_rate ) const;
  
};


#endif


