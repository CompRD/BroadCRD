// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#ifndef TABLE_H
#define TABLE_H



/*
 * Cambridge, July 9, 2001
 *
 * Classes
 *  table (template)
 *
 * Template class for tables. This is mainly used in the Report generation,
 * and the idea is that an implemenation of a table is an object with
 * n rows and m columns (the 'entries',) and with optional title, and
 * labels for the rows and columns.
 *
 * Class T must have ==, <<, and >> operators.
 */
#include <algorithm>
#include <strstream>

#include "String.h"
#include "Vec.h"



/*
 * Template class table
 *
 * A table of classes class. Each table may have a title (one big row at
 * the top,) and a special first row/column, of "labels." Other than that,
 * it is a simple list of rows.
 *
 * Remark:
 *  column labels can be broken in multiple rows, by using the character
 *  "/" to separate entries (eg a/b/c will split the label in three rows,
 *  a, be, and c.)
 */
template <class T>
class table
{
public:
  /*
   * table
   * constructor
   * public
   *
   * Legenda:
   *  n_rows: number of rows;
   *  n_cols: number of columns.
   */
  table(int n_rows, int n_cols) :
    n_rows_ ( n_rows ),
    n_cols_ ( n_cols )
    {
      // Must have at least one row and one column.
      Assert ( n_rows > 0 );
      Assert ( n_cols > 0 );
      
      // Deafults.
      title_.resize(0);
      row_labels_.clear();
      col_labels_.clear();
      
      // Resize.
      entries_.resize(n_rows);
      for (int jj=0; jj<n_rows; ++jj)
	entries_[jj].resize(n_cols);
    };
  
  void SetTitle(String title) { title_ = title; };
  
  void SetLabels(vec<String>& row_labels, vec<String>& col_labels) {
    row_labels_ = row_labels;
    col_labels_ = col_labels;
  };
  
  vec<T>& GetRowLabels() { return row_labels_; }

  vec<T>& GetColLabels() { return col_labels_; }

  void SetEntries(vec< vec<T> >& entries) { entries_ = entries; };
  
  vec< vec<T> >& GetEntries() { return entries_; }

  // If all the entries in a column (row) are the same, take off that
  //  column (row.) This might be useful to get rid of all 0's columns
  //  in report tables.
  void Compactify() {
    // Check columns (if there are two rows or more.)
    if ( n_rows_ > 1 ) {
      int jj = 0;
      while ( jj < n_cols_ ) {
	bool kill_column = true;

	for (int ii=1; ii<n_rows_; ++ii) {
	  if ( entries_[ii][jj] != entries_[ii-1][jj] ) {
	    kill_column = false;
	    break;
	  }	  
	}
      
	if ( kill_column ) {
	  if ( col_labels_.size() > 0 )
	    col_labels_.erase( col_labels_.begin() + jj );
	
	  for (int ii=0; ii<n_rows_; ++ii)
	    entries_[ii].erase( entries_[ii].begin() + jj );

	  --n_cols_;
	}
	else
	  ++jj;
      }
    }

    // Check rows (if there are two columns or more.)
    if ( n_cols_ > 1 ) {
      int ii = 0;
      while ( ii < n_rows_ ) {
	bool kill_row = true;

	for (int jj=1; jj<n_cols_; ++jj) {
	  if ( entries_[ii][jj] != entries_[ii][jj-1] ) {
	    kill_row = false;
	    break;
	  }	  
	}
      
	if ( kill_row ) {
	  if ( col_labels_.size() > 0 )
	    row_labels_.erase( row_labels_.begin() + ii );

          entries_.erase( entries_.begin() + ii );

	  --n_rows_;
	}
	else
	  ++ii;
      }
    }

  }

  // Sort rows by row labels (if row labels exist.)
  void SortTable() {
    // If there are no entries.
    if ( 0 == n_rows_ )
      return;

    // If there are no row labels.
    if ( n_rows_ != (int)row_labels_.size() )
      return;
    
    // Map and sort row labels (original position of labels.) Warning:
    //  the last row may be a 'total' row!
    map <String, int> label_to_id;
    for (int ii=0; ii<n_rows_; ii++)
      label_to_id[row_labels_[ii]] = ii;

    if ( "total" == row_labels_[n_rows_-1] )
      sort( row_labels_.begin(), row_labels_.end() - 1);
    else
      sort( row_labels_.begin(), row_labels_.end() );

    // Create a new (sorted) entries table.
    vec< vec< T > > temp_entries;
    temp_entries.reserve( entries_.size() );
    
    map<String, int>::iterator iter;
    for (int ii=0; ii<n_rows_; ii++) {
      iter = label_to_id.find( row_labels_[ii] );
      if ( iter != label_to_id.end() )
	temp_entries.push_back( entries_[ iter->second ] );      
    }
      
    // Replace entries_.
    entries_ = temp_entries;

  }

  // Remove a row (takes care of row labels.)
  void RemoveRow(int del_row) {
    if ( del_row < 0 || del_row > n_rows_-1 )
      return;
    
    entries_.erase( entries_.begin() + del_row );
    if ( row_labels_.size() > 0 )
      row_labels_.erase( row_labels_.begin() + del_row );

    --n_rows_;

  }

  // Remove a col (takes care of col labels.)
  void RemoveCol(int del_col) {
    if ( del_col < 0 || del_col > n_cols_-1 )
      return;

    for (int ii=0; ii<n_rows_; ii++)
      entries_[ii].erase( entries_[ii].begin() + del_col );

    if ( col_labels_.size() > 0 )
      col_labels_.erase( col_labels_.begin() + del_col );

    --n_cols_;

  }

  /*
   * table
   * GetPrintableTable
   * public
   *
   * Return:
   *  the full table, as a vector of vector of strings. The first
   *  row may contain only one entry, the title. All other rows
   *  have the same size. In case of errors, it returns an empty
   *  vector.
   *
   * Remark:
   *  the column labels can be broken in multiple rows (the same
   *  number of rows for all the labels!) The separator is the
   *  special character "/".
   */
  vec< vec<String > > GetPrintableTable() {
    // The result, given as a vector of rows.
    vec< vec<String> > result;

    // Validation: sizes.
    if ( n_rows_ != (int)row_labels_.size() )
      row_labels_.clear();
    if ( n_cols_ != (int)col_labels_.size() )
      col_labels_.clear();

    bool bad_data = false;
    if ( n_rows_ > 0 && n_rows_ != (int)entries_.size() )
      bad_data = true;
    if ( n_cols_ > 0 && !bad_data ) {
      for (int ii=0; ii<n_rows_; ++ii) {
	if ( n_cols_ != (int)entries_[ii].size() ) {
	  bad_data = true;
	  break;
	}
      }
    }

     /*
     if (bad_data)
     {    cout << "bad data:\n";
          for (int ii=0; ii<n_rows_; ++ii)
          {    for ( int jj = 0;  jj < (int) entries_[ii].size( ); jj++ )
               {    cout << entries_[ii][jj];
                    if ( jj < (int) entries_[ii].size( ) - 1 )
                         cout << ", ";    }
               cout << "\n";    }
          cout << "\n";    }
     */
        
    if ( bad_data )
      return result;

    // Reserve, not resize, result. 
    result.reserve( n_rows_ + 2 );

    // a_row: generic row.
    vec<String> a_row;

    // Title.
    if ( title_.size() > 0 ) {
      a_row.clear();
      a_row.push_back( title_ );
      result.push_back( a_row );
    }

    // Column labels.
    if ( col_labels_.size() > 0 ) {
      // Copy col_labels_ onto a local version.
      vec<String> temp_labels( col_labels_.size() );
      for (int ii=0; ii<(int)col_labels_.size(); ii++)
	temp_labels[ii] = col_labels_[ii];

      // The character "/" breaks a labels in separate rows. Do not
      //  worry about error checking, because Before() takes care of it.
      while ( temp_labels[0].Contains("/") ) {
	a_row.clear();
	a_row.reserve( 1 + (int)temp_labels.size() );
	
	if ( row_labels_.size() > 0 )
	  a_row.push_back("");
	
	for (int ii=0; ii<(int)temp_labels.size(); ++ii)
	  a_row.push_back( temp_labels[ii].Before("/") );
	
	result.push_back( a_row );
	
	for (int ii=0; ii<(int)temp_labels.size(); ++ii)
	  temp_labels[ii] = temp_labels[ii].After( "/" );
      }

      // There are no more "/"'s (last entry.)
      a_row.clear();
      a_row.reserve( 1 + (int)temp_labels.size() );
      
      if ( row_labels_.size() > 0 )
	a_row.push_back("");
      
      for (int ii=0; ii<(int)temp_labels.size(); ++ii)
	a_row.push_back( temp_labels[ii] );
	
      result.push_back( a_row );
    }

    // Entries and row labels.
    for (int ii=0; ii<n_rows_; ++ii) {
      a_row.clear();
      a_row.reserve( 1 + (int)col_labels_.size() );

      if ( row_labels_.size() > 0 )
	a_row.push_back( row_labels_[ii] );

      for (int jj=0; jj<n_cols_; ++jj) {
	String sz_entry;
	strstream str;
	str << entries_[ii][jj] << endl;
	getline( str, sz_entry );
	a_row.push_back( sz_entry );
      }
      
      result.push_back( a_row );
    }

    // All rigth, can return.
    return result;
  };


private:
  /*
   * table private variables:
   *
   * title_: title string. It may be empty (no title;)
   * row_labels_: row labels (it may be empty;)
   * col_labels_: column labels (it may be empty;)
   * entries_: data in the table;
   * n_rows_: number of rows;
   * n_cols_: number of columns.
   */
  String title_;
  vec<String> row_labels_;
  vec<String> col_labels_;
  vec< vec< T > > entries_;
  int n_rows_;
  int n_cols_;
};



#endif
