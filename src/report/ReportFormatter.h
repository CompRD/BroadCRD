// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef REPORT_FORMATTER_H
#define REPORT_FORMATTER_H

// class report_formatter
//
// This is an abstract superclass defining an interface through which you 
// can print elements of a report in a specific format such as LaTeX or HTML.

#include <iostream>

#include "String.h"
#include "Vec.h"

class report_formatter
{
public:
  report_formatter( const String &report_dir,
		    const String &report_name,
		    const String &report_pretty_name ) :
    report_dir_ ( report_dir ),
    report_name_ ( report_name ),
    report_pretty_name_ ( report_pretty_name ),
    ostream_ptr_( 0 ) { };

  virtual
  ~report_formatter() { };

  virtual void
  PrintDocumentBeginning(  const String& title ) = 0;

  // Headers.
  virtual void
  PrintSectionHeader(  const String& section_title ) = 0;

  virtual void
  PrintSubSectionHeader(  const String& sub_section_title ) = 0;

  virtual void
  PrintSubSubSectionHeader(  const String& sub_sub_section_title ) = 0;

  // Layout.
  virtual void
  PageBreak(  ) = 0;

  virtual void
  ParagraphBreak(  const bool indented = true ) = 0;

  // Font.
  virtual void
  UseTinyFont(  ) = 0;

  virtual void
  UseSmallFont(  ) = 0;

  virtual void
  UseNormalFont(  ) = 0;

  virtual void
  UseLargeFont(  ) = 0;

  virtual void
  StartPreformattedText() = 0;

  virtual void
  EndPreformattedText() = 0;

  // Text.
  virtual void
  PrintText(  const String& text ) = 0;

  virtual void
  PrintPreformattedText(  const String& preformatted_text ) = 0;

  /*
   * PrintPlusOrMinus
   *
   * Prin the plus or minus symbol..
   */
  virtual void
  PrintPlusOrMinus( ) = 0;

  /*
   * PrintList
   *
   * Appends a list to the document in the subclass-specific format.
   * The input is given as a vec<vec<String>>, because each list
   * can be on its own way a list. In other words, PrintList can
   * print a list of lists (only two levels, though: maybe I will
   * make it more general sometime in the future.)
   *
   * Remark: the second level list is interpreted as follows: the
   * first entry is a line of explanation, and then the indented
   * list starts.
   *
   * Legenda:
   *  list: list of lists to be appended.
   */
  virtual void
  PrintList(  const vec< vec< String > >& list ) = 0;

  /*
   * PrintLegenda
   *
   * Appends a legenda to the document in the proper format.
   * A legenda is akin to a list, only it is more compact, and
   * it only has one 'level' (instead of two.)
   *
   * Legenda:
   *  title: title of the legenda;
   *  legenda: the actual legenda.
   */
  virtual void
  PrintLegenda( const String &title, const vec< String > &list ) = 0;

  /*
   * PrintTable
   *
   * Appends a table to the document in the subclass-specific format.
   * The input is given as a vec<vec<String>>, one row at a time.
   *
   * Remark: the first vector (table[0]) must contain only one entry,
   * which will become the title of the table. All the other vectors
   * must have the same number of entries.
   *
   * Legenda:
   *  table: table to be appended;
   *  justification: a vec of 'r', 'c', or 'l' (columns justification;)
   *  rows_in_col_labels: how many rows in thge column labels.
   */
  virtual void
  PrintTable(  const vec< vec< String > >& table,
	       const vec< String >& justification,
	       const int rows_in_col_labels = 1 ) = 0;
  
  /*
   * Print a figure.
   */
  virtual void
  PrintFigureReference(  const String& figure_location ) = 0;

  /*
   * Print two figures, side by side.
   */
  virtual void
  PrintSideBySideFiguresReference( const String& figure_location_1,
				   const String& figure_location_2 ) = 0;

  // Prints a list of numbers with arrows between them pointing right.
  virtual void
  PrintSequenceWithArrows(  const vec<int>& sequence ) = 0;

  // Prints a full known contig alignment.
  virtual void
  PrintKnownContigAlignment( const String& title,
			     const vec<String>& alignments,
			     const vec< vec<int> >& link_lists,
			     const vec< vec<int> >& control_lists ) = 0;
  
  virtual void 
  PrintDocumentEnd(  ) = 0;
  
  // Get dir and report (relative, i.e. not full) names.
  virtual void
  GetFileNames(String& dir_name, String& report_name) = 0;
  
  /*
   * GetPrettyName
   *
   * Legenda:
   *  the short name for the report, as it will appear in the report
   *  itself (eg: "ArachneReport.pdf".)
   */
  virtual void
  GetPrettyName(String& report_pretty_name) = 0;

  /*
   * PrintDocumentLayout
   *
   * Final touches (system calls, etc.)
   *
   * Legenda:
   *  pre_dir: full pre directory;
   *  data_dir: data directory (as a subdirectory of pre_dir;)
   *  run_dir: run directory (as a subdirectory of data_dir.)
   */
  virtual void
  PrintDocumentLayout( String& pre_dir, String& data_dir,
		       String& run_dir, String& sub_dir ) = 0;

  void SetReportStreamPtr( ostream* report_stream_ptr )
    { ostream_ptr_ = report_stream_ptr; }      

protected:
  ostream& ReportStream()
    {
      Assert( ostream_ptr_ != 0 );
      return *ostream_ptr_;
    }

protected:
  String report_dir_;
  String report_name_;
  String report_pretty_name_;
  
private:
  ostream* ostream_ptr_;
};

#endif
