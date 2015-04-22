// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include "CoreTools.h"
#include "StringInt.h"

int BinPosition( const vec<string_int>& v, const String& x )
{    if ( v.size( ) == 0 ) return -1;
     int first = 0, last = (int) v.size( ) - 1, next;
     while (1)
     {    if (first == last) 
	  return ( !(x < v[last].s) && !(v[last].s < x) ) ? last : -1;
          next = first + (last - first) / 2;
          if ( x < v[next].s ) last = next;
          else if ( v[next].s < x ) first = next + 1;
          else return next;    }    }
