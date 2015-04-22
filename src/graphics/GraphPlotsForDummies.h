// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef GRAPHPLOT_FOR_DUMMIES_H
#define GRAPHPLOT_FOR_DUMMIES_H



/*
 * A simple wrapper for BuildHistogram plus BarGraph, which fully
 * automatize the process of getting an histogram.
 */
#include "graphics/BarGraph.h"
#include "String.h"
#include "Vec.h"

using namespace ns_psplot;



/*
 * labels: label for figure (as in BarGraph);
 * x_points[ii], y_points[ii]: coordinates of point ii;
 * file_name: where to save to (it will be appended an extra ".eps").
 * 
 * Returns 1 if all right, 0 for failure.
 */
int GenerateGraphPlot( const vec<freetext> &labels,
		       const vec<float> &x_points,
		       const vec<float> &y_poimts,
		       const String &file_name );



# endif
