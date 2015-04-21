///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//
//  There comes a time in every young assembly's life when you're going to want to 
//  assemble the unassembled reads as their own project and incorporate that back into
//  the original assembly via MergeHaplotypes. 
//
//  Previously this would not work as the imported reads were still valid in the
//  original assembly, but would have different data (i.e. lengths, etc.) as a 
//  consequence of the different trimming that occurred during the "unplaced assembly"
//   and the Assembly class would assert.  
//
//  ASSERT_IF_DIFFERENT (default behavior): reproduces current behavior that asserts
//    if you try to add a read to an assembly that is valid in that assembly, but 
///   has different data than in that assembly.
//   
//  KEEP_EXISTING : Do not assert, add the read, but use the read data as it appears
//    in the original assembly.
//
//  REPLACE_EXISTING:  DO NOT USE!!!!  Your code will barf and you will be unhappy.
//    This is not implemented and probably should not be, but that's under consideration.
//





#ifndef ADDBEHAVIOR
#define ADDBEHAVIOR


enum AddBehavior { ASSERT_IF_DIFFERENT, KEEP_EXISTING, REPLACE_EXISTING };


#endif 
