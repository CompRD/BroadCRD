// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// Col n1 ... nk: read from standard input, select columns n1 through nk,
// print them in order, separated by spaces.  The first column is column 1.
// Col t n1 ... nk: same as above, except that only tabs are counted as 
// separators.

#include "MainTools.h"
#include "math/Functions.h"

int main(int argc, char *argv[])
{
    bool tab_mode = false;

    vec<int> cols;
    for (int i = 0; i < argc - 1; i++)
    {
        String s(argv[i+1]);
        if (i == 0 && s == "t")
        {
            tab_mode = true;
        }
        else if (!s.IsInt() || s.Int() <= 0)
        {
            std::cerr << "column number doesn't make sense" << endl;
            return 1;
        }
        else
        {
            cols.push_back(s.Int());
        }
    }
    
    int n = Max(cols);
    vec<String> v(n);

    while (std::cin)
    {
        String line;
        getline(std::cin, line);
    
        if (line.size() > 0)
        {
            size_t col_start = 0;
            int col_num = 0;
            for (size_t i = 0; i < line.size() && col_num < n; i++)
            {
                if ((!tab_mode && isspace(line[i])) || line[i] == '\t')
                {
                    v[col_num] = line.substr(col_start, i - col_start);
                    col_num++;
                    col_start = i + 1;

                    // skip multiple separators unless we're in tab-mode
                    while (!tab_mode && isspace(line[col_start]))
                    {
                        col_start++;
                    }
                    i = col_start - 1;
                }                
            }
            
            if (col_num < n)
            {
                v[col_num] = line.substr(col_start, line.size() - col_start);
            }
            
            for (size_t j = 0; j < cols.size(); j++)
            {
                std::cout << (j > 0 ? " " : "" ) << v[cols[j] - 1];
            }
            std::cout << "\n";
        }
    }

    return 0;
}
