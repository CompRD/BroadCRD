// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef MACROSET_H
#define MACROSET_H

#include "String.h"
#include <map>
#include <ostream>

// A MacroSet contains macros, which are a name-replacement pair, in
// the sense that when the MacroSet is passed an expression to expand,
// any occurance of $name is replaced by replacement.  For example, if
// the macro "name" is assigned the expression "Bob" (via the MacroSet
// method add( "name", "Bob" ) ), the result of calling expand on "My
// name is $name" would be the string "My name is Bob".

// Both $name and ${name} refer to the same macro.  The expand method
// will treat anything that starts with a $ as a macro, with the name
// extending to the first non-letter, non-underscore character.  With
// brackets, the name extends from the opening bracket to the closing
// bracket.

// Using the bracket notation, macros can be used to construct other
// macros' names, e.g. given $first_name = 'Bob' and $select = 'first',
// an expansion of "${${select}_name}", would return "Bob".

// Macros that are undefined in the MacroSet are left as is, dollar
// sign intact.

// Dollar signs may be escaped with a backslash, preventing the
// MacroSet from attempting to find the associated macro.  The
// backslash is removed in this process.  For example, "I spent \$2.50
// on that hot dog." becomes "I spent $2.50 on that hot dog."

class MacroSet {
 public:
  void add(const String& macro_name, String macro_expression);

  String expand(const String& expression) const;
  void expandInPlace(String& expression) const;

  bool defined( const String &macro_name ) const
  {
    return ( macros_.count( macro_name ) > 0 );
  }

  void print(std::ostream& out) const;

 private:
  std::map<String,String> macros_;
};

#endif
