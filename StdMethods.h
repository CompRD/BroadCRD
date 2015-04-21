/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header: StdMethods.h

   Defines macros that take a list of fields of a class, and define standard
   methods of the class (copy constructor, assignment operator, binary read and write
   operators).

   For a class with fields x, y and z, put the line

   > STD_METHODS3(x,y,z);

   into the definition of the class, to create the standard methods.
   Use the appropriate macro for the number of fields in your class.

   *IMPORTANT*: If the class has a parent, use an STD_METHODSPn instead of
   STD_METHODSn macro!  Otherwise, the copy constructor and assignment operator
   will be incorrect.
*/

#ifndef __INCLUDE_StdMethods_h
#define __INCLUDE_StdMethods_h

#define GARBAGE_SEMICOLON_EATER(T)                              \
  typedef int __ ## T ## _dummy_eat_semicolon_

#define STD_METHODS1(T,f1)                           \
   T( const T& t ): f1(t.f1) { }                     \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     return *this;                                   \
   }                                                 \
  GARBAGE_SEMICOLON_EATER(T)

#define STD_METHODS2(T,f1,f2)                        \
   T( const T& t ): f1(t.f1), f2(t.f2) { }           \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     return *this;                                   \
   }                                                 \
  GARBAGE_SEMICOLON_EATER(T)

#define STD_METHODS3(T,f1,f2,f3)                     \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3) { } \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     return *this;                                   \
   }                                                 \
  GARBAGE_SEMICOLON_EATER(T)

#define STD_METHODS4(T,f1,f2,f3,f4)                  \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4) { }                                  \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     return *this;                                   \
   }                                                 \
  GARBAGE_SEMICOLON_EATER(T)

#define STD_METHODS5(T,f1,f2,f3,f4,f5)               \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5) { }                        \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     return *this;                                   \
   }                                                 \
  GARBAGE_SEMICOLON_EATER(T)

#define STD_METHODS6(T,f1,f2,f3,f4,f5,f6)            \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5), f6(t.f6) { }              \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     f6 = t.f6;                                      \
     return *this;                                   \
   }                                                 \
  GARBAGE_SEMICOLON_EATER(T)

#define STD_METHODS7(T,f1,f2,f3,f4,f5,f6,f7)         \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5), f6(t.f6), f7(t.f7) { }    \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     f6 = t.f6;                                      \
     f7 = t.f7;                                      \
     return *this;                                   \
   }                                                 \
  GARBAGE_SEMICOLON_EATER(T)


#define STD_METHODS8(T,f1,f2,f3,f4,f5,f6,f7,f8)      \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5), f6(t.f6), f7(t.f7),       \
       f8(t.f8) { }                                  \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     f6 = t.f6;                                      \
     f7 = t.f7;                                      \
     f8 = t.f8;                                      \
     return *this;                                   \
   }                                                 \
  GARBAGE_SEMICOLON_EATER(T)


#define STD_METHODS9(T,f1,f2,f3,f4,f5,f6,f7,f8,f9)   \
   T( const T& t ): f1(t.f1), f2(t.f2), f3(t.f3),    \
       f4(t.f4), f5(t.f5), f6(t.f6), f7(t.f7),       \
       f8(t.f8), f9(t.f9) { }                        \
   T& operator=( const T& t ) {                      \
     f1 = t.f1;                                      \
     f2 = t.f2;                                      \
     f3 = t.f3;                                      \
     f4 = t.f4;                                      \
     f5 = t.f5;                                      \
     f6 = t.f6;                                      \
     f7 = t.f7;                                      \
     f8 = t.f8;                                      \
     f9 = t.f9;                                      \
     return *this;                                   \
   }                                                 \
  GARBAGE_SEMICOLON_EATER(T)


#endif
// #ifndef __INCLUDE_StdMethods_h
