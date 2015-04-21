// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef LATEX_REPORT_FORMATTER_H
#define LATEX_REPORT_FORMATTER_H

#include "report/ReportFormatter.h"

// class latex_report_formatter
//
// This is a class through which you can print elements of a report in 
// LaTeX format.

class latex_report_formatter : public report_formatter
{
 public:
  latex_report_formatter( const String &report_dir,
			  const String &report_name,
			  const String &report_pretty_name ) :
    report_formatter( report_dir, report_name, report_pretty_name ) { } 
  
  // The default dtor is sufficient.
  //  ~latexreport_formatter();
  
  void
  PrintDocumentBeginning( const String& title );

  void
  PrintSectionHeader( const String& section_title );

  void
  PrintSubSectionHeader( const String& sub_section_title );

  void
  PrintSubSubSectionHeader( const String& sub_sub_section_title );

  void
  PageBreak();

  void
  PageBreakIfPicasBelow( const int picas );

  void
  ParagraphBreak( const bool indented = true );

  void
  UseTinyFont();

  void
  UseSmallFont();

  void
  UseNormalFont();

  void
  UseLargeFont();

  void
  StartPreformattedText();

  void
  EndPreformattedText();

  void
  PrintText( const String& text );

  void
  PrintPreformattedText( const String& preformatted_text );

  void
  PrintPlusOrMinus();

  void
  PrintList( const vec< vec< String > >& list );

  void
  PrintLegenda( const String &title, const vec< String > &legenda );

  void
  PrintTable( const vec< vec< String > >& table,
	      const vec<String>& justifications,
	      const int rows_in_col_labels = 1 );

  void
  PrintFigureReference( const String& figure_location );
  
  void
  PrintSideBySideFiguresReference( const String& figure_location_1,
				   const String& figure_location_2 );
  
  void
  PrintSequenceWithArrows( const vec<int>& sequence );

  void
  PrintKnownContigAlignment( const String& title,
			     const vec<String>& alignments,
			     const vec< vec<int> >& link_lists,
			     const vec< vec<int> >& control_lists );
  
  void 
  PrintDocumentEnd();

  void
  GetFileNames(String& dir_name, String& report_name);
  
  void
  GetPrettyName(String& report_pretty_name);
  
  void
  PrintDocumentLayout( String& pre_dir,
		       String& data_dir,
		       String& run_dir,
		       String& sub_dir );
  
private:

  const String 
  LatexSafeString( const String& a_string);

};

#endif
