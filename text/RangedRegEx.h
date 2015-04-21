// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef RANGEDREGEX_H
#define RANGEDREGEX_H

#include "String.h"

// SubstituteIntegerRanges() looks for substrings of the given
// expression of the form "[#number-number]" and expands them in place
// to extended regular expressions that would match that range of
// numbers.  For example, "[#6-12]" would be expanded to
// "(([6-9])|(1[0-2]))".

void SubstituteIntegerRanges( String& expression );

#endif
