#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <iostream>
#include <string>
#include <stdexcept>

#include "Vec.h"
#include "MainTools.h"
#include "graphics/Whiteboard.h"

int main(int argc, char* argv[])
{
    RunTime();

    BeginCommandArguments;
    CommandArgument_String(IN);
    CommandArgument_String(OUT);
    CommandArgument_Bool_OrDefault(ROW_NAMES, true);
    CommandArgument_Bool_OrDefault(COL_NAMES, true);
    CommandArgument_Bool_OrDefault(LOG2,      false);
    EndCommandArguments;

    if (ROW_NAMES == false)
    {
        printf("Lack of row names not yet supported.\n");
        assert(0);
    }
    if (COL_NAMES == false)
    {
        printf("Lack of column names not yet supported.\n");
        assert(0);
    }

    // 1. Load the data and compute the color bins.
    int line_number = 0;
    int nrows = 0;
    int ncols = 0;

    Ifstream(input, IN);
    vec<String>        row_names;
    vec<String>        column_names;
    vec< vec<double> > data;
    while (! input.eof())
    {
        char buf[8192];
        memset(buf, 0x00, 8192);
        input.getline(buf, 8192);
        String line(buf);
        vec<String> tokens;
        Tokenize(line, tokens);

        if (tokens.size() == 0) { break; } // we're done.

        if (COL_NAMES && line_number == 0)
        {
            column_names = tokens;
            ncols        = column_names.size();
            line_number += 1;
            continue;
        }

        if (ROW_NAMES)
        {
            row_names.push_back(tokens[0]);
            tokens.erase(tokens.begin());
        }

        vec<double> row(tokens.size());
        assert(row.size() == column_names.size());
        for (unsigned int i = 0; i < row.size(); i++)
        {
            row[i] = tokens[i].Double();
        }
        data.push_back(row);

        line_number += 1;
        nrows += 1;
    }

    double min = data[0][0];
    double max = data[0][0];
    for (unsigned int row = 0; row < data.size(); row++)
    {
        for(unsigned int column = 0; column < data[row].size(); column++)
        {
            if (LOG2)
            {
                data[row][column] = log2(data[row][column]);
            }

            if (::isinf(data[row][column])) { continue; }

            if (data[row][column] < min) { min = data[row][column]; }
            if (data[row][column] > max) { max = data[row][column]; }
        }
    }
    double intensity_scale = 1.0/(max - min);

    cout << min << " " << max << " " << intensity_scale << endl;

    // 2. Draw the image.
    double d  = 10;
    double Sx = 0 * d;
    double Sy = 0 * d;

    double total_size_x = (2*Sx) + (ncols*d);
    double total_size_y = (2*Sx) + (nrows*d);
    
    ns_whiteboard::whiteboard wb; 
    for (unsigned int row = 0; row < data.size(); row++)
    {
        for(unsigned int column = 0; column < data[row].size(); column++)
        {
            double intensity;
            if (::isinf(data[row][column])) intensity = 0;
            else                          intensity = (data[row][column] - min) * intensity_scale;

            //cout << data[row][column] << " " << min << " " << max << " " << intensity_scale << " " << intensity << endl;

            ns_whiteboard::xy_coords top_left((Sx + column*d), 
                                              (Sy + row*d));
            ns_whiteboard::xy_coords bottom_right((Sx + (column+1)*d), 
                                                  (Sy + (row+1)*d));
            if (::isinf(data[row][column]))
            {
                color c(1, 0, 0);
                ns_whiteboard::rect *pixel = new ns_whiteboard::rect(top_left, bottom_right, c);
                wb.Add(pixel);
            }
            else
            {
                color c(intensity, intensity, intensity);
                ns_whiteboard::rect *pixel = new ns_whiteboard::rect(top_left, bottom_right, c);
                wb.Add(pixel);
            }
        }
    }

    Ofstream(out, OUT);
    ns_whiteboard::ps_display display(out, total_size_x, total_size_y);
    wb.DisplayOn(&display);    
    out.close();
}

