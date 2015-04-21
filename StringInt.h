// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// string_int: a class with a string and an int in it, which is easy to sort
// and search by the string.

#ifndef STRING_INT_H
#define STRING_INT_H

#include "CoreTools.h"

class string_int 
{
  public:
     String s;
     int i;

     string_int( ) { }
     string_int( const String& s_arg, int i_arg )
       : s(s_arg), i(i_arg)
     { }

     friend bool operator<( const string_int& si1, const string_int& si2 )
     {    return si1.s < si2.s;    }

  friend ostream & operator<<(ostream & os, const string_int & s) {
    os << s.s << " " << s.i << " ";
    return os;
  }

};

int BinPosition( const vec<string_int>& v, const String& x );

#endif
