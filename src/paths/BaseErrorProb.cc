///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/BaseErrorProb.h"
#include "VecUtilities.h"
#include "MainTools.h"
#include "system/System.h"

/**
   \file
   \ingroup grp_quals 
   \copydoc BaseErrorProb.h
*/

//--------------------------------------------------
// BaseErrorTable Class:

BaseErrorProb::BaseErrorProb(const vec<bool>& base_map, double error_probability) {
  error_prob = error_probability;
  read_length = base_map.size();
  for (unsigned int i = 0; i < base_map.size(); ++i) {
    if (base_map[i])
      base_pos.push_back(i);
  }
}

BaseErrorProb::BaseErrorProb(const string& base_map, double error_probability) {
  error_prob = error_probability;
  read_length = base_map.size();
  for (unsigned int i = 0; i < base_map.size(); ++i) {
    if (base_map[i] == '1')
      base_pos.push_back(i);
  }
}

 
string BaseErrorProb::GetBaseMap() const {
  string bit_vec_str(read_length,'0');
  for (unsigned int i = 0; i < base_pos.size(); ++i) {
    bit_vec_str[base_pos[i]] = '1';
  }
  return bit_vec_str;
}


ostream& operator<<( ostream& o, const BaseErrorProb& p )
{
  return o << p.GetBaseMap() << " " << p.error_prob << endl;    
}

istream& operator>>( istream& i,  BaseErrorProb& p )
{
  string base_map;
  double error_prob;
  i >> base_map >> error_prob;
  p = BaseErrorProb(base_map, error_prob);
  return i;
}

istream& operator>>( istream& i,  BaseErrorTable& t )
{
  BaseErrorProb p;
  int length;
  i >> length;
  for ( int j = 0; j < length; j++ ) {
    i >> p;
    t.push_back(p);
  }
  return i;
}


//--------------------------------------------------
// BaseErrorProbProfile Class:

BaseErrorProbProfile::BaseErrorProbProfile(string filename) {
  ifstream in;
  OpenIfstream( in, filename );
  
  double p;
  while( in >> p  )
    probs.push_back( p/100.0 );
  
  table_size = probs.size();
}

BaseErrorProbProfile::BaseErrorProbProfile(const qualvector& quals) {
  table_size = quals.size();
  probs.reserve(table_size);
  double p;
  for (unsigned int i = 0; i < table_size; ++i) {
    p = pow(10, double(quals[i])/(-10.0));
    probs.push_back( p );
  }
}

BaseErrorTable BaseErrorProbProfile::getErrorTable(unsigned int num_bases, 
						   unsigned int min_errors, 
						   unsigned int max_errors, 
						   const unsigned int max_entries,
						   const bool peek_ahead) const {
  if (num_bases == 0) {
    // Use entire probability profile
    num_bases = table_size;
  }
  
  if (table_size < num_bases) {
    cout << "Requested base error table size " << num_bases << "\n"
	 << "is larger than probability profile size " << table_size << "\n"
	 << "Aborting.\n";
    exit(1);
  }

  if (max_errors > num_bases) {
    cout << "Warning: max_errors reset to requested number of bases: " 
	 << num_bases << "\n\n";
    max_errors = num_bases;
  }

  if (min_errors > max_errors) {
    cout << "Warning: min_errors greater than max_errors. Reset to: 0\n\n ";
    min_errors = 0;
  }


  BaseErrorTable error_table;
  double pleft=1.0;
  for(unsigned int errs = min_errors; errs <= max_errors; errs++) {
    pleft -= prob(num_bases, errs, error_table, probs);

    if (max_entries && error_table.size() >= max_entries) {
      // Met or exceeded required number of entries. Stop.
      break;
    }
  }

  // Sort the table in reverse order, most probable first
  sort(error_table.begin(), error_table.end(), greater<BaseErrorProb>());

  // find most probable max_error+1 entry and use as last entry in table
  if (peek_ahead && max_entries == 0) {
    vec<double> sorted_probs(probs);
    sorted_probs.resize(num_bases);
    vec<int> index(num_bases, vec<int>::IDENTITY); 
    SortSync(sorted_probs, index);

    double this_prob = 1.0;
    for(unsigned int i = 0; i < num_bases - (max_errors + 1); i++)
      this_prob *= 1.0-sorted_probs[i] ;
    for(unsigned int i = num_bases - (max_errors + 1); i < num_bases; i++) 
      this_prob *= sorted_probs[i];

    int trim_to;
    for (trim_to = error_table.isize() ; trim_to > 0; trim_to--)
      if (error_table[trim_to - 1].error_prob > this_prob)
	break;
    error_table.resize(trim_to);
    
    vec<bool> err;
    err.resize(num_bases, 0);
    for(unsigned int i = num_bases - (max_errors + 1); i < num_bases; i++) 
      err[index[i]] = 1;
    error_table.push_back(BaseErrorProb(err, this_prob));
  }

  if (max_entries && error_table.size() > max_entries) {
    // Exceeded required number of entries. Trunc table.
    error_table.erase(error_table.begin() + max_entries, error_table.end());
  }

  // Return the table
  return error_table;
}

/**
   For each possible combination of \p num_error error positions in reads of length
   \p num_bases, compute the probability of that combination based in the probabilities
   of errors at individual read positions as given in \p per_base_err_rate.

   \param[in] num_bases length of reads considered here
   \param[in] num_errors we consider combinations of exactly this many errors,
      that is, errors occurring at this many positions in the read simultaneously
   \param[in] per_base_err_rate for each read position, the individual probability
      of error at that position
   \param[out] error_table to this table we add an entry for each combination of \p num_errors
      errors, representing the combination and the probability of its occurrence.

   \return the total probability that \e some combination of \p num_errors errors will
      occur in a read of length \p num_bases, given the probabilities of errors at individual
      positions in \p per_base_err_rate
*/
double BaseErrorProbProfile::prob(const int num_bases, const int num_errors,
			       BaseErrorTable& error_table,
			       const vec<double>& per_base_err_rate ) const {
    
  double tot_prob = 0.0, this_prob;
  int n = num_bases;
    
  // Step through all n-choose-num_errors ways of having the errors.
    
  vec<bool> err;
  FirstCombination(err, n, num_errors);
  do {
    this_prob = 1.0;
    for(int i=0; i<n; i++)
      this_prob *= (err[i] ? per_base_err_rate[i] : 1.0-per_base_err_rate[i] );
    tot_prob += this_prob;
    error_table.push_back(BaseErrorProb(err, this_prob));
  } while( NextCombination(err) );
    
  return tot_prob;
}
