///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This program accepts a coordinate-sorted BAM that represents the alignment
// of assembly contigs to a reference by BWASW. BWASW represents "chimeric"
// alignments by reporting multiple alignments for a particular read. In the
// case of aligning assembly contigs to the reference, these chimeric alignments
// indicate that a contig is missing a large chunk of sequence that is in
// the reference, probably indicating a biological deletion in the sample(s)
// used in the assembly. This program reports all such gaps, provided that
// the segments are placed consecutively on the reference and they are
// consistently oriented.

#include "MainTools.h"
#include "lookup/SAM.h"

int main(int argc, char** argv)
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(ASSEMBLYBAM, "BAM file of aligned assembly "
        "contigs");
    EndCommandArguments;
    
    SAM::BAMFile assembly_bam(ASSEMBLYBAM);
    SAM::Record read_record;
    
    std::string prev_ref_contig;
    std::string prev_assembly_contig;
    uint prev_start = 0;
    uint prev_end = 0;
    bool prev_reversed = false;
    
    while (assembly_bam.nextRecord(read_record))
    {
        if (read_record.isMapped())
        {
            std::string current_ref_contig = read_record.getRefName();
            std::string current_assembly_contig = read_record.getQueryName();
            SAM::Alignment read_alignment(read_record,assembly_bam.getLogger());
            uint current_start = read_alignment.getRefStart();
            uint current_end = read_alignment.getRefEnd();
            bool current_reversed = read_record.isReversed();
            
            if (prev_ref_contig == current_ref_contig &&
                prev_assembly_contig == current_assembly_contig)
            {
                if (prev_reversed == current_reversed &&
                    prev_end < current_start)
                {
                    std::cout << current_ref_contig << ":" << prev_end << "-" <<
                        current_start << "\n";
                }
                else if (prev_start > current_start)
                {
                    std::cerr << ASSEMBLYBAM << " is not coordinate sorted.\n";
                    exit(1);
                }
            }
            
            prev_ref_contig = read_record.getRefName();
            prev_assembly_contig = read_record.getQueryName();
            prev_start = current_start;
            prev_end = current_end;
            prev_reversed = current_reversed;
        }
    }
    
    return 0;
}
