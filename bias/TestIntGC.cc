///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * A program to measure the GC content of a set of test intervals on
 * a reference (useful for understanding the test intervals supplied
 * to BadCoverage).
 *
 * Original author: Michael G. Ross <mgross@broadinstitute.org>
 *
 */

#include "Basevector.h"
#include "Bitvector.h"
#include "FetchReads.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "VecUtilities.h"

void load_test_intervals(vec<unsigned int>& test_chr,
                         vec<unsigned int>& test_start,
                         vec<unsigned int>& test_stop,
                         vec<String> names, String test_intervals,
                         String ref_fasta)
{
    fast_ifstream testin(test_intervals);
    String line;
    while (!testin.fail())
    {
        getline(testin, line);
        if (!testin.fail())
        {
            line.GlobalReplaceBy(",", "");
            test_start.push_back(line.Between(":", "-").Int());
            test_stop.push_back(line.RevAfter("-").Int());
            int pos = Position(names, line.Before(":"));
            if (pos < 0)
            {
                cerr << "Failure to process TEST_INTERVALS "
                     << test_intervals << endl;
                cerr << "unable to locate chromosome " << line.Before(":")
                     << " in " << ref_fasta << endl;
                CRD::exit(EXIT_FAILURE);
            }
            test_chr.push_back(pos);
        }
    }
    SortSync(test_chr, test_start, test_stop);
    return;    
}

void load_fasta_info(vecbitvector& refamb, vec<String>& names, String ref_fasta)
{
    fast_ifstream rin(ref_fasta);
    vec<char> bases;
    String line;
    while (!rin.fail())
    {
        getline(rin, line);
        if (rin.fail() || (bases.nonempty() && line.Contains(">", 0)))
        {
            bitvector b(bases.size(), False);
            for (size_t i = 0; i < bases.size(); i++)
            {
                if (!Base::isBase(bases[i]))
                {
                    b.Set(i,True);
                }
            }
            refamb.push_back_reserve(b);
            bases.clear();
        }
        if (!rin.fail())
        {
            if (line.Contains(">", 0))
            {
                String n = line.After(">");
                if (n.Contains(" "))
                {
                    n = n.Before(" ");
                }
                names.push_back(n);
            }
            else
            {
                for (size_t i = 0; i < line.size(); i++)
                {
                    bases.push_back(line[i]);
                }
            }
        }
    }

    return;
}

int main(int argc, char** argv)
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(TEST_INTERVALS, "file of intervals, one line "
                               "per interval, form is chr*:*-*");
    CommandArgument_String_Doc(REF_FASTA, "reference file to get GC info from");
    EndCommandArguments;

    vecbitvector refamb;
    vec<String> names;
    vecbasevector genome;
    load_fasta_info(refamb, names, REF_FASTA);
    FetchReads(genome, 0, REF_FASTA);
    
    vec<unsigned int> test_chr;
    vec<unsigned int> test_start;
    vec<unsigned int> test_stop;
    load_test_intervals(test_chr, test_start, test_stop, names, TEST_INTERVALS,
                        REF_FASTA);

    longlong total_bases = 0;
    longlong total_gc = 0;
    vec<int> gc_counts(test_chr.size());
    vec<int> base_counts(test_chr.size());
    for (size_t i = 0; i < test_chr.size(); i++)
    {
        longlong amb_bases = 0;
        for (unsigned int loc = test_start[i]; loc < test_stop[i]; loc++)
        {
            if (refamb[test_chr[i]][loc])
            {
                amb_bases++;
            }
        }
        gc_counts[i] = genome[test_chr[i]].GcBases(test_start[i], test_stop[i])
            - amb_bases;
        base_counts[i] = (test_stop[i] - test_start[i]) - amb_bases;
        total_gc += gc_counts[i];
        total_bases += base_counts[i];
    }

    cout << TEST_INTERVALS << " report" << endl;
    cout << "total bases: " << total_bases << endl;
    cout << "total GC: " << total_gc << " "
         << (100 * double(total_gc) / double(total_bases)) << "%" << endl;

    double interval_mean = 0;
    for (size_t i = 0; i < gc_counts.size(); i++)
    {
        interval_mean += double(gc_counts[i]) / double(base_counts[i]);
    }
    interval_mean /= gc_counts.size();

    double interval_var = 0;
    for (size_t i = 0; i < gc_counts.size(); i++)
    {
        interval_var += pow(double(gc_counts[i]) / double(base_counts[i])
                            - interval_mean, 2);
    }
    interval_var /= (gc_counts.size() - 1);

    cout << "interval GC% mean: " << interval_mean << endl;
    cout << "interval GC% var: " << interval_var << endl;
    cout << "interval GC% stddev: " << sqrt(interval_var) << endl;

    return 0;
}
