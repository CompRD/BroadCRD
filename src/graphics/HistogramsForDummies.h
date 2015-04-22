// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef HISTOGRAMS_FOR_DUMMIES_H
#define HISTOGRAMS_FOR_DUMMIES_H



/*
 * Cambridge, September 14, 2001
 *
 * HistogramsForDummies.h
 *
 * A simple wrapper for BuildHistogram plus BarGraph, which fully
 * automatize the process of getting an histogram.
 */
#include "graphics/BarGraph.h"
#include "String.h"
#include "Vec.h"


/*
 * GenerateHistogram
 *
 * It creates an .eps bar graph figure, with the histogram of
 * the lengths given in input.
 *
 * This is basically a wrapper for David's BarGraph class, plus the
 * BuildHistogram class.
 *
 * Legenda:
 *  labels: label for figure (as in BarGraph;)
 *  lengths: vector of lengths;
 *  file_name: full path name of the file where the figure will be saved onto;
 *  cutoff: cutoff parameter (see remark;)
 *  mean: given mean;
 *  stdev: given standard deviation.
 *
 * Return:
 *  it returns 1, if it can create the bar graph file; or else it
 *  return 0 if any problem occurs.
 *
 * Remark:
 *  file_name will be appended a ".eps".
 *
 * Remark:
 *  if cutoff is = 0, then no cut-off will be performed. If cut-off is
 *  a positive value, then GenerateHistogram will reject all lengths
 *  such that | length - mean | >= cutoff*stdev. GenerateHistogram will
 *  calculate mean and stdev if they are not passed in input. If file_data
 *  is not empty, GenerateHistogram saves also the numeric data.
 */
int GenerateHistogram( vec<ns_psplot::freetext> labels,
		       const vec<float>& lengths,
		       String file_name,
		       String file_data = "",
		       float cutoff = 0,
		       float mean = 0,
		       float stdev = 0,
                       const int digits_to_right = 4 );
int GenerateHistogram( vec<ns_psplot::freetext> labels,
		       const vec<int>& lengths,
		       String file_name,
		       String file_data = "",
		       float cutoff = 0,
		       float mean = 0,
		       float stdev = 0,
                       const int digits_to_right = 4 );



# endif
