// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef PARSECONFIG
#define PARSECONFIG

#include <fstream>

#include "String.h"
#include "system/Types.h"
#include "Vec.h"

void TranslateReadName( vec<String>& statements, String read_name,
     vec<String>& variables, vec<String>& values, Bool& excluded, Bool& unmatched,
     vec<String>& global_vars, vec<String>& global_vals );

void ParseConfigFile( istream& in, vec<String>& statements );

#endif
