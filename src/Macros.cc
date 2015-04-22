// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include "Macros.h"
#include "system/Assert.h"
#include <iostream>

// Perform macro substitutions in expression.

void MacroSet::expandInPlace( String& expression ) const 
{
  int pos = 0;
  while ( pos < (int) expression.size() ) 
  {
    if ( expression[pos] != '$' ) 
    {
      ++pos;
      continue;
    }

    if ( pos > 0 && expression[pos-1] == '\\' ) 
    {
      expression.erase( (unsigned int) pos-1, 1 );
      continue;
    }
    
    int macro_begin;
    int macro_end;
    int delimiter_size;
    int expression_size = expression.size();
    
    if ( pos+1 == expression_size ) 
      return;

    String found_macro;
    int found_macro_size;

    if ( expression[pos+1] == '{' ) 
    {
      macro_begin = pos+2; // the macro name begins after '${'
      for ( macro_end = macro_begin; 
	    macro_end < expression_size;
	    ++macro_end )
	if ( expression[macro_end] == '}' )
	  break;
      
      if ( macro_end == expression_size )
	return;

      // extract macro name from expression
      found_macro_size = macro_end - macro_begin;
      found_macro = expression.substr(macro_begin, found_macro_size);
      
      this->expandInPlace( found_macro );

      delimiter_size = 3; // '$' '{' and '}'
    }

    else
    {
      macro_begin = pos+1; // the macro name begins after '$'
      for ( macro_end = macro_begin; 
	    macro_end < expression_size;
	    ++macro_end )
	if ( ! isalnum( expression[ macro_end ] ) &&
	     ! ( '_' == expression[ macro_end ] ) )
	  break;

      found_macro_size = macro_end - macro_begin;
      found_macro = expression.substr(macro_begin, found_macro_size);

      delimiter_size = 1; // '$'
    }

    // search through list of known macros
    std::map<String,String>::const_iterator macro_iter;
    macro_iter = macros_.find( found_macro );

    if ( macro_iter != macros_.end() )
    {
      const String& replacement = macro_iter->second;
      expression.replace( (unsigned int) pos, 
                          found_macro_size+delimiter_size, 
                          replacement.c_str() ); 	// replace it
      pos += replacement.size();
    }

    else
      ++pos;
  }
}

String MacroSet::expand( const String& expression ) const {
  String expansion(expression) ;
  expandInPlace(expansion);
  return expansion;
}

void MacroSet::add(const String& macro_name, String expansion) {
  Assert( ! macro_name.empty() );
  macros_[macro_name] = expand(expansion);
}

void MacroSet::print(std::ostream& out) const {
  std::map<String,String>::const_iterator macro_iter;
  for ( macro_iter = macros_.begin(); 
	macro_iter != macros_.end(); 
	macro_iter++ ) {
    out << "'" 
	<<  macro_iter->first << "' expands to '" 
	<< macro_iter->second << "'." << std::endl;
  }
}
  
