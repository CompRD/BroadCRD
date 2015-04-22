// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef READTYPES
#define READTYPES

#ifdef NOTNOTNOTNOTNOTNOTNOTNOTNOTNOTNOTNOTNOTNOTNOTNOT

#include "String.h"
#include "Vec.h"

// See also ParseConfig.cc, ParseConfig.h.

// First we define the variables which are allowed in the reads.config file,
// and the values they are allowed to take on.  

// The following definition provides the list of variables which are allowed
// in the reads.config file, and the values which they are allowed to take on.

// Here is how to read it:
//
// An EnumVariablei can take on any of the i listed values, i = 1, 2, 3.   [1 byte]
// An NonnegVariable can take on any nonnegative integer value.            [4 bytes]
// A FourCharVariable is a string of 1 to 4 alphanumeric characters.       [4 bytes]
// A StringVariable is an arbitrary string.                                [2 bytes]

// These should be ordered so that word boundaries are maintained.

#define READNAMEVARIABLES                                        \
     NonnegVariable(   Plate )                                   \
     NonnegVariable(   InsertLength )                            \
     NonnegVariable(   Instance )                                \
     NonnegVariable(   InsertLengthDev )                         \
     FourCharVariable( Well )                                    \
     StringVariable(   Library )                                 \
     EnumVariable2(    Chemistry, dye_primer, dye_terminator )   \
     EnumVariable3(    Orientation, forward, reverse, none )     \
     EnumVariable1(    Class, PairedWholeGenome )

// Current size: 4 * 7 bytes.

// Then we define the legal combinations of variables which a read can parse to.

// This file defines a routine which -- upon being passed a variable and a
// value -- decides whether they are legal, and translates both to integers.

class string_pointer {
     public:

     static vec<String> contents;

     unsigned short loc;

     string_pointer( ) { }

     string_pointer( const String& s )
     {    if ( contents.size( ) == 65535 )
          {    cout << "\nToo many strings have been encountered as values "
                    << "in the reads.config file.  I though this would never "
                    << "happen.\nTo fix, class string_pointer will need to be "
                    << "redesigned (which should not be hard).\n";
               exit(1);    }
          loc = contents.size( );
          contents.push_back(s);    }

     String AsString( ) const { return contents[loc]; }

};

//////////////////// PASS ONE DEFINITIONS /////////////////////

// These definitions are used to build the class readinfo.

#define EnumVariable1( NAME, VAL1 ) char NAME;
#define EnumVariable2( NAME, VAL1, VAL2 ) char NAME;
#define EnumVariable3( NAME, VAL1, VAL2, VAL3 ) char NAME;
#define NonnegVariable( NAME ) int NAME;
#define FourCharVariable( NAME ) int NAME;
#define StringVariable( NAME ) string_pointer NAME;

class readinfo {
     public:

     READNAMEVARIABLES

     readinfo( )
     {

          #undef EnumVariable1
          #undef EnumVariable2
          #undef EnumVariable3
          #undef NonnegVariable
          #undef FourCharVariable
          #undef StringVariable

          #define EnumVariable1( NAME, VAL1 ) NAME = -1;
          #define EnumVariable2( NAME, VAL1, VAL2 ) NAME = -1;
          #define EnumVariable3( NAME, VAL1, VAL2, VAL3 ) NAME = -1;
          #define NonnegVariable( NAME ) NAME = -1;
          #define FourCharVariable( NAME ) NAME = -1;
          #define StringVariable( NAME ) NAME = -1;

     };

};

#undef EnumVariable1
#undef EnumVariable2
#undef EnumVariable3
#undef NonnegVariable
#undef FourCharVariable
#undef StringVariable

//////////////////// PASS TWO DEFINITIONS /////////////////////

#define EnumVariable1( NAME, VAL1 )                                           \
     else if ( variable == ##NAME )                                           \
     {    if ( value == ##VAL1 ) read.NAME = 0;                               \
          else                                                                \
          {    cout << "\nThe variable " << variable << " was assigned the "  \
                    << "illegal value " << value << ".\n";                    \
               cout << "The only legal value for this variable is "           \
                    << ##VAL1 << ".\n";                                       \
               exit(1);    }    }

#define EnumVariable2( NAME, VAL1, VAL2 )                                     \
     else if ( variable == ##NAME )                                           \
     {    if ( value == ##VAL1 ) read.NAME = 0;                               \
          else if ( value == ##VAL2 ) read.NAME = 1;                          \
          else                                                                \
          {    cout << "\nThe variable " << variable << " was assigned the "  \
                    << "illegal value " << value << ".\n";                    \
               cout << "The legal values for this variable are " << ##VAL1    \
                    << " and " << ##VAL2 << ".\n";                            \
               exit(1);    }    }

#define EnumVariable3( NAME, VAL1, VAL2, VAL3 )                               \
     else if ( variable == ##NAME )                                           \
     {    if ( value == ##VAL1 ) read.NAME = 0;                               \
          else if ( value == ##VAL2 ) read.NAME = 1;                          \
          else if ( value == ##VAL3 ) read.NAME = 2;                          \
          else                                                                \
          {    cout << "\nThe variable " << variable << " was assigned the "  \
                    << "illegal value " << value << ".\n";                    \
               cout << "The legal values for this variable are " << ##VAL1    \
                    << ", " << ##VAL2 << ", and " << ##VAL3 << ".\n";         \
               exit(1);    }    }

#define NonnegVariable( NAME )                                                \
     else if ( variable == ##NAME ) read.NAME = value.Int( );

#define FourCharVariable( NAME )                                              \
     else if ( variable == ##NAME )                                           \
     {    if ( value.size( ) > 4 )                                            \
          {    cout << "\nThere is an upper limit of four characters on "     \
                    << "size of a variable called " << ##NAME                 \
                    << ", but you've assigned it the value " << value         \
                    << ".\n";                                                 \
               exit(1);    }                                                  \
          read.NAME = 0;                                                      \
          if ( value.size( ) >= 1 ) read.NAME = value[0];                     \
          if ( value.size( ) >= 2 ) read.NAME ^= (value[1] << 8);             \
          if ( value.size( ) >= 3 ) read.NAME ^= (value[2] << 16);            \
          if ( value.size( ) >= 4 ) read.NAME ^= (value[3] << 24);    }

#define StringVariable( NAME )                                                \
     else if ( variable == ##NAME ) read.NAME = value;

void ModifyReadInfo( String variable, String value, readinfo& read )
{    if ( 0 == 1 );
     READNAMEVARIABLES
     }

//////////////////////////////////////////////////////////////

/*
inline TranslateVariableAndValue( const String& variable, const String& value,
     int& var, int& val )
{
     if ( 0 == 1 );
     EnumVariable2( Chemistry, alt, terminator )

          ...

// Read type 1.

     if ( defined(read.Chemistry) 
          && defined(read.Instance)
          && defined(read.Orientation)
          && defined(read.Well)
          && defined(read.Plate)
          && defined(read.InsertLength)
          && defined(read.InsertLengthDev)
          && read.Class.AsString( ) == "PairedWholeGenome" )
*/

#endif

#endif
