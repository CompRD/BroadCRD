///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Align two sets of sequences in fasta files F1 and F2.

#include "Alignment.h"
#include "Basevector.h"
#include "MainTools.h"
#include "FetchReads.h"
#include "pairwise_aligners/LocalAlign.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "lookup/LookAlign.h"

#include <set>

int main(int argc, char *argv[])
{
    RunTime();

    BeginCommandArguments;
    CommandArgument_String_Doc(F1, "A FASTA of sequences (reads)");
    CommandArgument_String_Doc(F2, "Another FASTA of sequences (reference "
        "contigs)");
    CommandArgument_String_OrDefault_Doc(ALIGNER, "SWFREE", "Choose the SWFREE " 
        "(SmithWatFree), SWAFFINE(SmithWatAffine) or LOCAL aligner");
    CommandArgument_UnsignedInt_OrDefault_Doc(CAP, 0, "If non-zero, cap all "
        "sequences at this length.");
    CommandArgument_Int_OrDefault_Doc(MATCH, 1, "Local aligner's match "
        "score");
    CommandArgument_UnsignedInt_OrDefault_Doc(INDEL, 3, "Penalty for an indel "
        "in a SmithWatFree or Local alignment");
    CommandArgument_UnsignedInt_OrDefault_Doc(MISMATCH, 2, "Penalty for a "
        "mismatch in a SmithWatFree, SmithWatAffine or Local alignment");
    CommandArgument_UnsignedInt_OrDefault_Doc(GAP_OPEN, 3, "Penalty for "
        "opening a gap in SmithWatAffine");
    CommandArgument_UnsignedInt_OrDefault_Doc(GAP_EXTEND, 1, "Penalty for "
        "extending a gap in SmithWatAffine");
    CommandArgument_Bool_OrDefault_Doc(BRIEF, False, "Output short parseable "
        "format");
    CommandArgument_Bool_OrDefault_Doc(VISUAL, True, "Output human-readable "
        "visualization of each alignment");
    CommandArgument_Bool_OrDefault_Doc(ABBREVIATE, True, "Abbreviate alignment "
        "visualization");
    CommandArgument_Bool_OrDefault_Doc(PARSEABLE, True, "Output a machine "
        "parseable representation of each alignment");
    CommandArgument_Bool_OrDefault_Doc(NAMES, False, "Use the sequence names "
        "in the F1 and F2 FASTA files in the output");
    CommandArgument_String_OrDefault_Doc(SELECT_NAMES,"", "A file of names to "
        "select from F1 and align - implies NAMES.");
    CommandArgument_Bool_OrDefault_Doc(ONE_TO_ONE, False, "Align each read "
        "from F1 to the corresponding read from F2 only");
    CommandArgument_Double_OrDefault_Doc(MAX_MISMATCH_RATE, 1.0, "If there "
        "are more mismatches than this per base, throw out the alignment");
    CommandArgument_Double_OrDefault_Doc(MAX_ERROR_RATE, 1.0, "If there "
        "are more errors(mis, ins or dels) than this per base, throw out "
        "the alignment");  
    CommandArgument_Int_OrDefault_Doc(MIN_MATCHES, 0, "Minimum number of "
        "matching bases for alignment to be reported");
    CommandArgument_Bool_OrDefault_Doc(FW_OR_RC, True, "If True, try aligning "
        "each F1 read and its reverse complement, and report only the better "
        "of the two alignments");
    CommandArgument_Bool_OrDefault_Doc(BEST_ONLY, False, "Report only the "
        "F2 pairing for each read in F1");
    CommandArgument_Bool_OrDefault_Doc(FULL2, False, 
        "Report that alignment extends to both ends of second sequence.");
    EndCommandArguments;

    ForceAssert(ALIGNER == "LOCAL" || ALIGNER == "SWFREE"
        || ALIGNER == "SWAFFINE");

    if (MISMATCH >= 2 * INDEL)
    {
        FatalErr("DoAlign requires that the MISMATCH penalty be "
            "less than twice the INDEL penalty.");
    }

    vecbasevector b1, b2;
    vecString * names1 = 0;
    vecString * names2 = 0;
    vecqualvector q;

    std::set<String> selectedNames;
    if (!SELECT_NAMES.empty())
    {
        NAMES=true;
        String name;
        Ifstream(is, SELECT_NAMES);
        while (is)
        { 
            is >> name;
            selectedNames.insert(name); 
        }
    }

    if (NAMES)
    {
        names1 = new vecString();
        names2 = new vecString();
    }

    FetchReads(b1, q, names1, 0, F1);
    FetchReads(b2, q, names2, 0, F2);

    if (CAP > 0)
    {    
        for (size_t i1 = 0; i1 < b1.size(); i1++)
        {
            b1[i1].Cap(CAP);
        }
        for (size_t i2 = 0; i2 < b2.size(); i2++)
        {
            b2[i2].Cap(CAP);
        }
    }

    if (ONE_TO_ONE && b1.size() != b2.size())
    {
        FatalErr("Sizes of fasta files must be the same for ONE_TO_ONE option");
    }

    align a;
    alignment al;
    if (BRIEF)
    {
        if (NAMES)
        {
            std::cout << "BRIEF\tname1\tname2\t";
        }
        std::cout << "id1\tid2\tstart1\tend1\tlength1\tstart2\tend2\tlength2\t"
             << "mismatches\tinsertions\tdeletions\n";
    }
    String name;
    for (size_t i1 = 0; i1 < b1.size(); i1++)
    {
        //check for correct names if we are selecting.
        if (!selectedNames.empty())
        {
            name = (*names1)[i1];
        }
        
        if (b1[i1].size() > 0 && (selectedNames.empty() || 
            selectedNames.find(name) != selectedNames.end()))
        {
            static basevector b1_fw_rc;
            const basevector& b1_or_rc = (!FW_OR_RC ? b1[i1] : b1_fw_rc);
            ostringstream alignout;
            int best_errs = numeric_limits<int>::max();
            
            for (size_t i2 = ONE_TO_ONE ? i1 : 0;
                (ONE_TO_ONE && i2 == i1) || (!ONE_TO_ONE && i2 < b2.size());
                i2++)
            {
                if (b1[i1].size() <= b2[i2].size()) //avoid refsize < readsize
                {
                    if (ALIGNER == "LOCAL") 
                    {
                        LocalAlign(b1[i1], b2[i2], a, MATCH, MISMATCH, INDEL);
                    }
                    else if (ALIGNER == "SWFREE")
                    {
                        int best_loc;
                        SmithWatFree(b1[i1], b2[i2], best_loc, al, false, 
                            false, MISMATCH, INDEL);
                        a.UnpackFrom(al);    
                    }
                    else if (ALIGNER == "SWAFFINE")
                    {
                        SmithWatAffine(b1[i1], b2[i2], al, false, false, 
                            MISMATCH, GAP_OPEN, GAP_EXTEND);
                        a.UnpackFrom(al);
                    }
                    else
                    {
                        FatalErr(ALIGNER + " is an unknown aligner.");
                    }

                    if ( FULL2 
                        && ( al.pos2( ) > 0 || al.Pos2( ) < b2[i2].isize( ) ) )
                        continue;
        
                    if (FW_OR_RC)
                    {    
                        static align arc;
                        static alignment alrc;
                        static basevector b1rc;
                        b1rc.ReverseComplement(b1[i1]);
                        if (ALIGNER == "LOCAL")
                        { 
                            LocalAlign(b1rc, b2[i2], arc, MATCH,
                                MISMATCH, INDEL);
                        }
                        else if (ALIGNER == "SWFREE")
                        {    
                            int best_loc;
                            SmithWatFree(b1rc, b2[i2], best_loc, alrc, false, 
                                false, MISMATCH, INDEL);
                            arc.UnpackFrom(alrc);    
                        }
                        else if (ALIGNER == "SWAFFINE")
                        {
                            SmithWatAffine(b1rc, b2[i2], alrc, false, false,
                                MISMATCH, GAP_OPEN, GAP_EXTEND);
                            arc.UnpackFrom(alrc);
                        }
                        else
                        {
                            FatalErr(ALIGNER + " is an unknown aligner");
                        }
                        
                        static vec<int> errsfw;
                        errsfw = a.MutationsGap1Gap2(b1[i1], b2[i2]);
                        long weighted_errsfw_sum = errsfw[0] * MISMATCH
                            + errsfw[1] * INDEL + errsfw[2] * INDEL;

                        static vec<int> errsrc;
                        errsrc = arc.MutationsGap1Gap2(b1rc, b2[i2]);
                        long weighted_errsrc_sum = errsrc[0] * MISMATCH
                            + errsrc[1] * INDEL + errsrc[2] * INDEL;
                        
                        if (weighted_errsfw_sum <= weighted_errsrc_sum)
                        {
                            b1_fw_rc = b1[i1];
                        }
                        else
                        {    
                            b1_fw_rc = b1rc;
                            a = arc;    
                        }    
                    }
        
                    int len1 = b1[i1].size();
                    int len2 = b2[i2].size();
                    static vec<int> errs;
                    errs = a.MutationsGap1Gap2(b1_or_rc, b2[i2]);
                    if (float(errs[0]) / float(a.extent1()) <=
                            MAX_MISMATCH_RATE &&
                        float(errs[0]+errs[1]+errs[2]) / a.extent1() <=
                            MAX_ERROR_RATE &&
                        a.MatchingBases(b1_or_rc, b2[i2]) > MIN_MATCHES)
                    {                                                    
                        long weighted_errs_sum = errs[0] * MISMATCH
                            + errs[1] * INDEL + errs[2] * INDEL;
                            
                        if (BEST_ONLY && weighted_errs_sum < best_errs)
                        {
                            alignout.clear();
                            alignout.str("");
                        }
                            
                        if (!BEST_ONLY || weighted_errs_sum < best_errs)
                        {                
                            if (NAMES)
                            {
                                alignout << (*names1)[i1] << "\t"
                                    << (*names2)[i2] << "\t";
                            }
                            if (BRIEF)
                            {
                                alignout << "BRIEF\t" <<i1 << "\t" << i2 << "\t"
                                    << a.pos1()<< "\t" <<  a.Pos1() << "\t"
                                    << len1 << "\t"
                                    << a.pos2()<< "\t" <<  a.Pos2() << "\t"
                                    << len2 << "\t"
                                    << errs[0] << "\t" << errs[1] << "\t"
                                    << errs[2] << "\n"; 
                            } 
                            if (VISUAL)
                            {
                                alignout << "sequence " << i1;
                                if (b1_or_rc != b1[i1])
                                {
                                    alignout << "rc";
                                }
                                else
                                {
                                    alignout << "fw";
                                }
                                alignout << " vs sequence " << i2 
                                     << " insertions: " << errs[2] 
                                     << " deletions: " << errs[1]
                                     << " mutations: " << errs[0] 
                                     << " extent1: " << a.extent1() << "\n";
                                alignout << "a.pos1( ) = " << a.pos1()
                                    << ", a.Pos1( ) = " << a.Pos1()
                                    << ", len1 = " << len1 << "\n";
                                alignout << "a.pos2( ) = " << a.pos2()
                                    << ", a.Pos2( ) = " << a.Pos2()
                                    << ", len2 = " << len2 << "\n";
                                PrintVisualAlignment(ABBREVIATE, alignout,
                                                     b1_or_rc, b2[i2], a);    
                            }    
                            if (PARSEABLE)
                            {
                                look_align la;
                                la.a = a;
                                la.query_id = i1;
                                la.target_id = i2;
                                la.query_length = b1[i1].size();
                                la.target_length = b2[i2].size();
                                la.mutations = errs[0];
                                la.indels = errs[1] + errs[2];
                                la.rc1 = (b1_or_rc != b1[i1]);
                                la.PrintParseable(alignout, &b1[i1], &b2[i2]);
                            }
                        }
            
                        if (BEST_ONLY && weighted_errs_sum < best_errs)
                        {
                            best_errs = weighted_errs_sum;
                        }
                    }
                }
            } //end of inner for loop
            std::cout << alignout.str();
        }
    } //end of outer for loop
    
    delete names1;
    delete names2;
    
    return 0;
}
