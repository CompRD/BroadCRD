///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BuildJobs.
// Build job list for Human whole genome DISCOVAR runs.

#include "MainTools.h"
#include "system/RunTime.h"
#include "Basevector.h"
#include "feudal/VirtualMasterVec.h"
#include "math/Functions.h"


typedef VirtualMasterVec<BaseVec> VBaseVecVec;


uint32_t ComputeRegionCount(const uint32_t contig_size, const uint32_t region_size, const uint32_t min_overlap) {

  if (region_size >= contig_size) 
    return 1;
  else 
    return Max(ceil( (0.0 + contig_size - min_overlap)/(region_size - min_overlap) ), 1.0 );
}


int main(int argc, char *argv[])
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(JOB_LIST, "Job list file");
    CommandArgument_String_OrDefault_Doc(JOB_COMMAND, "run_region.sh", "Command to use in job list file");
    CommandArgument_String_OrDefault_Doc(JOB_ARGS, "", "Additional args to be passed when running JOB_COMMAND");
    CommandArgument_UnsignedInt_OrDefault_Doc(REGION_SIZE, 200, "Split genome into regions of this size (kb)");
    CommandArgument_UnsignedInt_OrDefault_Doc(REGION_OVERLAP, 10, "Minimum overlap between regions (kb)");
    CommandArgument_String_OrDefault_Doc(REF_GENOME, "/wga/scr4/bigrefs/human19", "Target genome");
    CommandArgument_Bool_OrDefault_Doc(HUMAN, true, "Genome is human - use X and Y for 23 and 24")
    EndCommandArguments;


    const String ref_fastb = REF_GENOME + "/genome.fastb";
   
    const uint32_t region_size = REGION_SIZE * 1000;
    const uint32_t region_overlap = REGION_OVERLAP * 1000;

    // Create job control file
    Ofstream( out, JOB_LIST );

    // Parse contigs/chromosomes from reference genome fastb file
    VBaseVecVec vbvv(ref_fastb.c_str());
    uint32_t contigs = 0;
    uint32_t regions = 0;
    String contig_name;
    for (VBaseVecVec::const_iterator it = vbvv.begin(); 
	 it != vbvv.end(); it++) {
      contigs++;
      uint32_t contig_size = it->size();
      uint32_t region_count = ComputeRegionCount(contig_size, region_size, region_overlap);
      contig_name = ToString(contigs);
      if (HUMAN) {
	if (contigs == 23)
	  contig_name = "X";
	else if (contigs == 24)
	  contig_name = "Y";
	else if (contigs > 24)
	  break;
      } 
      cout << "Contig: " << contig_name << " Size: " << contig_size << " Regions: " << region_count << endl;
      for (uint32_t region = 0; region < region_count; region++) {
	uint32_t start = (region_size - region_overlap) * region;
	uint32_t end = start + region_size;
	// Treat last region as special - extend it back to make it full size
	if (region == region_count - 1) {
	  end = contig_size - 1;
	  start = (region == 0 ? 0 : end - region_size);
	} 
	out << JOB_COMMAND << " " << regions << " " << contig_name <<  ":" << start << "-" << end 
	    << " " << JOB_ARGS << endl;
	regions++;
      }
    }
    cout << regions << " regions divided between " << contigs << " contigs." << endl;
    cout << "Done." << endl;
    
}
