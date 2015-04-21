// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include "report/LatexReportFormatter.h"
#include "String.h"
#include "system/System.h"

void
latex_report_formatter::PrintDocumentBeginning( const String& title )
{
  ostream& report = ReportStream();
  
  report << "\\documentclass{article}\n"
	 << "\\usepackage{epsfig}\n"
	 << "\\usepackage{latexsym}\n"
	 << "\\usepackage{supertabular}\n"
	 << "\n"
	 << "\\batchmode\n"
	 << "\n"
	 << "\\textwidth 7.0in\n"
	 << "\\textheight 9in\n"
	 << "\\hoffset=-0.5in\n"
	 << "\\voffset=-0.5in\n"
	 << "\\oddsidemargin=0pt\n"
    //	 << "\\parskip=4pt plus 2pt\n"
	 << "\\author{}\n"
	 << "\\date{}\n"
	 << "\n"
	 << "\\title{Report for Arachne assembly of {\\tt "
	 << this->LatexSafeString(title) << "} }\n"
	 << "\n"
	 << "\\begin{document}\n"
	 << "\n"
	 << "\\maketitle\n"
	 << endl;
}

void
latex_report_formatter::PrintSectionHeader( const String& section_title )
{
  ostream& report = ReportStream();
  
  report << "\\section{" << this->LatexSafeString(section_title) << "}\n"
	 << endl;
}

void
latex_report_formatter::PrintSubSectionHeader( const String& sub_section_title )
{
  ostream& report = ReportStream();
  
  report << "\\subsection{" << this->LatexSafeString(sub_section_title) << "}\n"
	 << endl;
}

void
latex_report_formatter::PrintSubSubSectionHeader( const String& sub_sub_section_title )
{
  ostream& report = ReportStream();
  
  report << "\\subsubsection{" << this->LatexSafeString(sub_sub_section_title) << "}\n"
	 << endl;
}

void
latex_report_formatter::PageBreak( )
{
  ostream& report = ReportStream();
  
  report << "\\newpage\n"<< endl;
}

void
latex_report_formatter::PageBreakIfPicasBelow( const int picas )
{
  ostream& report = ReportStream();
  
  report << "{\\newdimen\\pageleft \\pageleft=-\\pagetotal"
	 << " \\advance\\pageleft by \\pagegoal\n \\ifdim\\pageleft < "
	 << picas << "pc \\pagebreak \\fi}" << endl;
}


void
latex_report_formatter::ParagraphBreak( const bool indented )
{
  ostream& report = ReportStream();
  
  if ( indented )
    report << "\\par" << endl;
  else
    report << "\\par\\noindent" << endl;
}

void
latex_report_formatter::UseTinyFont( )
{
  ostream& report = ReportStream();

  report << "\\scriptsize" << endl;
}
  
void
latex_report_formatter::UseSmallFont( )
{
  ostream& report = ReportStream();

  report << "\\small" << endl;
}
  
void
latex_report_formatter::UseNormalFont( )
{
  ostream& report = ReportStream();

  report << "\\normalsize" << endl;
}

void
latex_report_formatter::UseLargeFont( )
{
  ostream& report = ReportStream();

  report << "\\large" << endl;
}

void
latex_report_formatter::StartPreformattedText()
{
  ostream& report = ReportStream();

  report << "\n"
	 << "\\begin{verbatim}"
	 << endl;
}

void
latex_report_formatter::EndPreformattedText()
{
  ostream& report = ReportStream();

  report << "\\end{verbatim}"
	 << endl;
}

void
latex_report_formatter::PrintText( const String& text )
{
  ostream& report = ReportStream();
  
  report << this->LatexSafeString(text);
}

void
latex_report_formatter::PrintPreformattedText( const String& preformatted_text )
{
  ostream& report = ReportStream();
  
  report << "\n"
	 << "\\begin{verbatim}\n"
	 << preformatted_text
	 << "\\end{verbatim}"
	 << endl;
}

void 
latex_report_formatter::PrintPlusOrMinus()
{
  ostream& report = ReportStream();
  
  report << "$\\;\\pm\\;$";
}

void 
latex_report_formatter::PrintList( const vec< vec<String> >& list_data)
{
  ostream& report = ReportStream();
  
  report << "\\begin{itemize}\n";
  
  for (int ii=0; ii<(int)list_data.size(); ++ii)
  {
    report << "\\item\n"
	    << this->LatexSafeString(list_data[ii][0]) << "\n";
    
    if ( (int)list_data[ii].size() > 1 )
    {
      report << "\\begin{itemize}\n";
      for (int jj=1; jj<(int)list_data[ii].size(); ++jj)
	report << "\\item\n"
		<< this->LatexSafeString(list_data[ii][jj]) << "\n";
      report << "\\end{itemize}\n";
    }
  }
  
  report << "\\end{itemize}\n"
	 << endl;
}

void 
latex_report_formatter::PrintLegenda( const String & title,
				      const vec<String> &legenda)
{
  ostream& report = ReportStream();
  
  report << "\\begin{verse}\n";    

  report << "\n{\\em " << this->LatexSafeString( title ) << "}\\\\\n";

  for (int ii=0; ii<(int)legenda.size(); ++ii)
    report << "\\--- " << this->LatexSafeString( legenda[ii] ) << "\\\\\n";

  report << "\\end{verse}\n"
	 << endl;
}

// TODO - documentation! justifications here is wrong-wrong-wrong. Must come
//  back and adjust it.
void 
latex_report_formatter::PrintTable( const vec< vec<String> >& table,
				    const vec<String>& justifications,
				    const int rows_in_col_labels )
{
  ostream& report = ReportStream();
  
  // Validation (check dimensions.)
  int n_row = (int)table.size();
  if ( n_row < 2 )
    return;
  
  if ( 1 != (int)table[0].size() )
    return;

  int n_col = (int)table[1].size();
  for (int ii=1; ii<n_row; ++ii)
    if ( (int)table[ii].size() != n_col )
      return;

  PageBreakIfPicasBelow( rows_in_col_labels + 5 );

  // Start printing report.
  report << "\\begin{center}\n";
  report << "\\begin{small}\n";
  // supertabular environment: allow tables to split across pages.
  // (Should col labels appear at the top of each page?  No, for now...)
  report << "\\tablefirsthead{\\hline \\multicolumn{" << n_col << "}{|c|}{" 
	 << this->LatexSafeString( table[0][0] ) 
	 << "} \\\\ \\hline}\n";
  report << "\\tablehead{\\hline \\multicolumn{" << n_col << "}{|c|}{" 
	 << this->LatexSafeString( table[0][0] ) 
	 << ", {\\em continued}} \\\\ \\hline}\n";
  report << "\\tabletail{\\hline \\multicolumn{" << n_col 
	 << "}{|r|}{\\em continued on next page\\ldots} \\\\ }\n";
  report << "\\tablelasttail{\\hline}\n";

  report << "\\begin{supertabular}{|";
  for (int ii=0; ii<n_col; ++ii)
    report << justifications[ii] << "|";
  report << "} \\hline\n";

  // For each column, compute the maximum width of the heading, and the maximum
  // width of the body.  Then computed a recommended kern for the right side
  // of the column, so as to center it if it is right-justified.

  vec<int> head_width(n_col, 0), body_width(n_col, 0);
  vec<String> recommended_kern(n_col);
  for ( int ii = 0; ii < n_col; ii++ )
  {     for ( int jj = 1; jj <= rows_in_col_labels; jj++ )
             if ( !table[jj][ii].Contains( "$" ) )
                  head_width[ii] 
                       = max( head_width[ii], (int) table[jj][ii].size( ) );
        for ( int jj = rows_in_col_labels + 1; jj < n_row; jj++ )
             body_width[ii] = max( body_width[ii], (int) table[jj][ii].size( ) );
        if ( body_width[ii] < head_width[ii] )
        {    float kern = 2.1 * ( head_width[ii] - body_width[ii] );
             recommended_kern[ii] = "{\\kern" + ToString(kern, 1) + "pt}";    }    }
    
  // Title -- this is now taken care of by the supertabular environment instead.
  //  report << "\\multicolumn{" << n_col << "}{|c|}{"
  //	 << this->LatexSafeString( table[0][0] ) << "} \\\\ \\hline\n";

  // Cells.
  for (int ii=1; ii<n_row; ++ii) {
    for (int jj=0; jj<n_col; ++jj ) {
      // Column labels must be always centered.

      if ( ii <= rows_in_col_labels ) 
      { if ( !table[ii][jj].empty() ) report << "\\multicolumn{1}{|c|}{";      
          report << LatexSafeString( table[ii][jj] );
	  if ( !table[ii][jj].empty() ) report << "}"; }
      else
      { if ( justifications[jj] == "r" ) report << recommended_kern[jj];
        report << LatexSafeString( table[ii][jj] );
        if ( justifications[jj] == "r" ) report << recommended_kern[jj]; }

      if ( jj < n_col - 1 )
	report << " & ";
      else {
	report << " \\\\ ";
	if ( ii >= rows_in_col_labels )
	  report << "\\hline\n";
	else
	  report << "\n";
      }
    }
  }

  report << "\\end{supertabular}\n"
	 << "\\end{small}\n"
	 << "\\end{center}\n"
	 << endl;
}

void
latex_report_formatter::PrintFigureReference( const String& figure_location )
{
  ostream& report = ReportStream();
  
  report << "\\begin{center}\n"
	 << "\\epsfig{file=" << figure_location << ".eps}\n"
	 << "\\end{center}\n"
	 << endl;
}

void
latex_report_formatter::PrintSideBySideFiguresReference( const String& figure_location_1,
							 const String& figure_location_2 )
{
  ostream& report = ReportStream();
  
  report << "\\begin{center}\n"
	 << "\\mbox{\n"
	 << "\\epsfig{file=" << figure_location_1 << ".eps}\n"
	 << "\\epsfig{file=" << figure_location_2 << ".eps}\n"
	 << "}\n"
    	 << "\\end{center}\n"
	 << endl;
}

void
latex_report_formatter::PrintSequenceWithArrows( const vec<int>& sequence )
{
  ostream& report = ReportStream();
  
  report << "Links Between Alignments: ";

  unsigned int num_items = sequence.size();
  for ( unsigned int item_idx = 0; item_idx < num_items; ++item_idx )
  {
    if ( ( item_idx % 18 == 3 && item_idx > 20 ) ||
	 item_idx == 21 )
      report << "\n\\par" << endl;

    if ( sequence[item_idx] == -1 )
      report << "$\\cdots";
    else
      report << sequence[item_idx];
    
    if ( item_idx != num_items - 1 )
      report << "$\\kern2pt$\\rightarrow$\\kern2pt";
  }
  report << "\n\\par\\noindent" << endl;
}

void
latex_report_formatter::PrintKnownContigAlignment( const String& title,
						   const vec<String>& alignments,
						   const vec< vec<int> >& link_lists,
						   const vec< vec<int> >& control_lists)
{
  ostream& report = ReportStream();

  // Title.
  report << "\\vspace{0.1in}"
	 << "\\par\\noindent\\verb!"<< title << "!\n";

  // Alignments.
  for (int ii=0; ii<(int)alignments.size(); ii++)
    report << "\\par\\noindent\\verb|" << " " << alignments[ii] << "|\n";
 
  // Links. The double change of size is done to get the baselinestretch
  //  command to kick in. More or less the same for the extra blank lines.
  report << "\n"
	 << "\\normalsize\n"
	 << "\\renewcommand{\\baselinestretch}{1.25}\n"
	 << "\\scriptsize\n"
	 << "\n"
	 << "\\par\\noindent\\verb|Links:  |\n";
    
  for (int ii=0; ii<(int)link_lists.size(); ii++) {
    const vec<int> &the_list = link_lists[ii];
    const vec<int> &the_control = control_lists[ii];

    report << "{\\bf [}$";
    unsigned int num_items = the_list.size();
    for ( unsigned int item_idx = 0; item_idx < num_items; ++item_idx ) {
      if ( the_list[item_idx] == -1 )
	report << "\\cdots";
      else
	report << the_list[item_idx];
      
      if ( item_idx != num_items - 1 ) {
	// Replace right arrow with a slithering arrow, if the link
	//  is suspect.
	bool suspect_link = false;
	
	if ( the_control[item_idx] > 0 )
	  suspect_link = true;

	// Print.
 	if ( suspect_link ) 
          report << " \\kern3pt\\lower2.3pt\\hbox{\\LARGE$\\leadsto$}\\kern2pt ";
 	else
 	  report << " \\rightarrow ";
      }
    }
    
    if ( ii < (int)link_lists.size()-1 )
      report << "${\\bf ]}\n$\\;\\;\\;\\;\\;$\n";
    else
      report << "${\\bf ]}.\n";
  }
  
  // Back to regular linespace (same trick as before to get linestretch to
  //  kick in).
  report << "\n"
	 << "\\normalsize\n"
	 << "\\renewcommand{\\baselinestretch}{1}\n"
	 << "\\scriptsize\n"
	 << "\n";

}

void 
latex_report_formatter::PrintDocumentEnd()
{
  ostream& report = ReportStream();
  
  report << "\\end{document}" << endl;
}

void
latex_report_formatter::GetFileNames(String& dir_name, String& report_name)
{
  dir_name = report_dir_;
  report_name = report_name_;
}

void
latex_report_formatter::GetPrettyName(String& report_pretty_name)
{
  report_pretty_name = report_pretty_name_;
}

void 
latex_report_formatter::PrintDocumentLayout( String& pre_dir,
					     String& data_dir,
					     String& run_dir,
					     String& sub_dir )
{
  // Generate latex output (pre_dir and report_filename are full path names.)
  String full_run = pre_dir + "/" + data_dir + "/" + run_dir;
  String full_sub = full_run;
  if ( "" != sub_dir )
    full_sub += "/" + sub_dir;

  String full_report = full_sub + "/" + report_dir_;

  Remove( full_report + "/" + report_name_.Before( ".tex" )  + ".dvi" );
  String command = "cd " + full_report + "; latex "
    + report_name_ + " > /dev/null";
  char exit_value = System(command);
  if ( exit_value != 0 ) {
    cout << "\n"
	 << "Warning! System call to latex failed! There will be no report."
	 << "\n"
	 << endl;
    return;
  }

  // If latex fails, there is nothing to do.
  String temp_name = report_name_;
  String name_root = temp_name.Before( ".tex" );
  if ( !IsRegularFile( full_report + "/" + name_root + ".dvi" ) )
       FatalErr( "Latex did not generate a dvi file." );

  // Run dvips, and ps2pdf...
  command = "cd " + full_report + "; dvips " + name_root + ".dvi"
    + " -o " + report_pretty_name_ + " > " + name_root + ".ps.log 2>&1";
  if ( System(command) != 0 ) FatalErr( "System call to dvips failed." );
  
  command = "cd " + full_report + "; ps2pdf " + report_pretty_name_;
  System( command );
  
  // ...and move only .ps and .pdf files in the run directory (only
  //  if sub_dir is empty.)
  if ( "" == sub_dir ) {
    String start_ps = full_report + "/" + report_pretty_name_;
    String end_ps = full_run + "/" + report_pretty_name_;
    if ( IsRegularFile ( start_ps ) ) Mv( start_ps, end_ps );
    else FatalErr( "postscript file was not generated." );
    
    if ( ! report_pretty_name_.Contains( ".ps", -1 ) )
      cout << Date() << ": Warning! Could not generate pdf file!" << endl;
    else {
      String pdf_name = report_pretty_name_.Before( ".ps" ) + ".pdf";
      String start_pdf = full_report + "/" + pdf_name;
      String end_pdf = full_run + "/" + pdf_name;
      if ( IsRegularFile ( start_pdf ) )
	Mv( start_pdf, end_pdf );
      else
	cout << Date() << ": Warning! Could not generate pdf file!" << endl;
    }
  }

}

const String 
latex_report_formatter::LatexSafeString( const String& a_string)
{
  // Sante --- Tue Jul 10 11:20:29 EDT 2001
  // Can deal with almost all Latex reserved characters (with the exception
  // of: '$', '&', '~', '^', and '\',) and with simple math in-text formulae.

  // If every character were a reserved character, the new string would be
  // twice as long as the old.
  const int max_length = a_string.size() * 2;
  String safe_string;
  safe_string.resize(max_length);
  unsigned int safe_string_idx = 0;

  bool math_mode = false;
  for ( unsigned int i = 0; i < a_string.size(); i++ ) {
    char ch = a_string[i];

    // If find a '$', switch math mode on/off.
    if ( '$' == ch )
      math_mode = !math_mode;

    // Replace reserved characters, if not in math mode.
    if ( !math_mode )    
      if ( '%' == a_string[i] ||
	   '#' == a_string[i] ||
	   '_' == a_string[i] ||
	   '{' == a_string[i] ||
	   '}' == a_string[i] )
	safe_string[safe_string_idx++] = '\\';
    
    safe_string[safe_string_idx++] = a_string[i];
  }

  safe_string.resize(safe_string_idx);
  return safe_string;
}
