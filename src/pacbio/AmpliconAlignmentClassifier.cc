///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This program accepts a BAM/SAM of aligned sequences and a list of 
// amplicon intervals. It then classifies each sequence as being
// (a) from one of the amplicons, (b) from some other part of the genome, or
// (c) unaligned. The BAM/SAM must be sorted by queryname in case multiple
// alignments are given for each read (a la BWASW). Classification is made
// by choosing the amplicon interval with the greatest percent covered
// from the read in question. The read names are assumed to end with the
// string "FX" where X is the index of the amplicon that we expect the
// sequence to match. The output is a confusion matrix.

#include <algorithm>

#include "FastIfstream.h"
#include "MainTools.h"
#include "lookup/SAM.h"

void load_intervals(vec<String>& chr, vec<size_t>& start,
    vec<size_t>& stop, String intervals)
{
    fast_ifstream testin(intervals);
    String line;
    while (!testin.fail())
    {
        getline(testin, line);
        if (!testin.fail())
        {
            if (line.Contains("\t"))
            {
                line = line.Before("\t");
            }
            line.GlobalReplaceBy(",", "");
            chr.push_back(line.Before(":"));
            start.push_back(line.Between(":", "-").Int());
            stop.push_back(line.RevAfter("-").Int());
            if (start.back() > stop.back())
            {
                std::cerr << line << " is an illegal interval" << std::endl;
                CRD::exit(EXIT_FAILURE);
            }
        }
    }
        
    return;    
}

void score_result(vec<vec<int> >& confusion_matrix,
    const String& current_read, const size_t best_i, const bool print_correct,
    const bool print_incorrect)
{
    unsigned int barcode = current_read.RevAfter("F").Int();
    confusion_matrix[barcode][best_i]++;

    if ((print_correct && barcode == best_i) ||
        (print_incorrect && barcode != best_i &&
        best_i < confusion_matrix.size()))
    {
        std::cout << current_read.RevBefore("_F") << "\t" <<
            "barcode=" << barcode << "\t" << "amplicon=" << best_i << "\n";
    }
        
    if (best_i < confusion_matrix.size())
    {
        //std::cout << current_read << "\t" << "F" << best_i << "\n";
    }
    else if (best_i == confusion_matrix.size())
    {
        //std::cout << current_read << "\tnonamplicon\n";
    }
    else if (best_i == confusion_matrix.size() + 1)
    {
        //std::cout << current_read << "\tunaligned\n";
    }
    else
    {
        FatalErr("illegal choice of best amplicon");
    }
    
    return;
}

int main(int argc, char** argv)
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(ALIGNMENT, "BAM/SAM file of aligned amplicons, "
        "sorted by readname");
    CommandArgument_String_Doc(INTERVALS, "File of amplicon location ranges");
    CommandArgument_Bool_OrDefault_Doc(PRINT_CORRECT, False, "Print the names "
        " of correctly barcoded reads");
    CommandArgument_Bool_OrDefault_Doc(PRINT_INCORRECT, False, "Print the "
        " names of incorrectly barcoded reads");
    EndCommandArguments;

    SAM::BAMFile alignment_file(ALIGNMENT);
    
    if (!alignment_file.getHeader() ||
        alignment_file.getHeader()->getSortOrder() != "queryname")
    {
        std::cerr << "input BAM/SAM must be sorted by queryname and have "
            "a correctly set @HD field\n";
        return 1;
    }
    
    vec<String> intervals_chr;
    vec<size_t> intervals_start;
    vec<size_t> intervals_stop;
    load_intervals(intervals_chr, intervals_start, intervals_stop, INTERVALS);
    
    vec<vec<int> > confusion_matrix(intervals_chr.size());
    for (size_t c = 0; c < confusion_matrix.size(); c++)
    {
        confusion_matrix[c].resize(intervals_chr.size() + 2, 0);
    }
    
    SAM::Record rec;
    String current_read = "";
    size_t best_i = intervals_chr.size() + 1;
    double best_i_overlap = -1;
    
    while (alignment_file.nextRecord(rec))
    {
        String read_name = rec.getQueryName();
        
        if (current_read != read_name)
        {
            if (!current_read.empty())
            {
                score_result(confusion_matrix, current_read, best_i, 
                    PRINT_CORRECT, PRINT_INCORRECT);
            }

            best_i = intervals_chr.size() + 1; // unmapped
            best_i_overlap = -1;
            current_read = read_name;
        }
        
        if (rec.isMapped())
        {
            if (best_i == intervals_chr.size() + 1)
            {
                best_i = intervals_chr.size(); // not mapped to an amplicon
            }
            SAM::Alignment aln(rec, alignment_file.getLogger());
            String alignment_chr = rec.getRefName();
            size_t alignment_start = aln.getRefStart();
            size_t alignment_stop = aln.getRefEnd();
            
            for (size_t i = 0; i < intervals_chr.size(); i++)
            {
                if (alignment_chr == intervals_chr[i])
                {
                    size_t overlap_start = std::max(alignment_start, 
                        intervals_start[i]);
                    size_t overlap_stop = std::min(alignment_stop, 
                        intervals_stop[i]);
                    double overlap_fraction =
                        double(overlap_stop - overlap_start) /
                        (intervals_stop[i] - intervals_start[i]);
                    if (overlap_stop > overlap_start &&
                        overlap_fraction > best_i_overlap)
                    {
                        best_i_overlap = overlap_fraction;
                        best_i = i;
                    }
                }
            }
        }
    }

    if (!current_read.empty())
    {
        score_result(confusion_matrix, current_read, best_i, PRINT_CORRECT,
            PRINT_INCORRECT);
    }
    
    vec<int> totals(confusion_matrix.size(), 0);
    for (size_t b = 0; b < totals.size(); b++)
    {
        totals[b] = Sum(confusion_matrix[b]);
    }
    
    std::cout << "barcode\t";
    for (size_t c = 0; c < confusion_matrix[0].size() - 2; c++)
    {
        std::cout << c << "\t";
    }
    
    std::cout << "correct\tincorrect\tnonamplicon\tunaligned\ttotal\n";
    for (size_t b = 0; b < confusion_matrix.size(); b++)
    {
        std::cout << b << "\t";
        int incorrect_count = 0;
        // print confusion
        for (size_t c = 0; c < confusion_matrix[b].size() - 2; c++)
        {
            std::cout << (confusion_matrix[b][c] / 
                static_cast<double>(totals[b]))<< "\t";
            if (b != c)
            {
                incorrect_count += confusion_matrix[b][c];
            }
        }
        // print correct/incorrect percentages
        std::cout << (confusion_matrix[b][b] / static_cast<double>(totals[b]))
            << "\t";
        std::cout << (incorrect_count / static_cast<double>(totals[b])) << "\t";
        // print nonamplicon/nonaligned percentages
        for (size_t c = confusion_matrix[b].size() - 2;
            c < confusion_matrix[b].size(); c++)
        {
            std::cout << (confusion_matrix[b][c] /
                static_cast<double>(totals[b])) << "\t";
        }
        std::cout << totals[b] << "\n";
    }
    
    return 0;
}
