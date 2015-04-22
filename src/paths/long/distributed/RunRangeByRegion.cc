///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "MainTools.h"
#include "paths/long/fosmid/Fosmids.h"
#include "system/ForkManager.h"

uint32_t find_leftmost_region(uint32_t start, const uint32_t region_size, const uint32_t overlap) {

  return Max(floor((0.0 + start - overlap)/(region_size - overlap)),0.0);
}

uint32_t find_rightmost_region(uint32_t stop, const uint32_t region_size, const uint32_t overlap) {

  return floor((0.0 + stop)/(region_size - overlap));
}


typedef std::tuple<String, uint32_t, uint32_t>  Region;


vec<Region> convert_range_to_region_list( const String& chr_name, const uint32_t range_start, uint32_t range_stop,
					  const uint32_t region_size, const uint32_t region_overlap) {
 
  vec<Region> region_list;
  
  uint32_t left_region = find_leftmost_region(range_start, region_size, region_overlap);
  uint32_t right_region = find_rightmost_region(range_stop, region_size, region_overlap);
  uint32_t region_count = (right_region - left_region) + 1;
  for (uint32_t region = left_region; region <= right_region; region++) {
    uint32_t region_start = (region_size - region_overlap) * region;
    uint32_t region_stop = region_start + region_size;
    region_list.push_back(make_tuple(chr_name, region_start, region_stop));
  }
  return region_list;
}


String build_command( const String& work_dir, const String& chr_name, 
		      const uint32_t region_start, const uint32_t region_stop) {

  String range = chr_name + ":" + ToString(region_start) + "-" + ToString(region_stop);
  String out_head = work_dir + "/" + "aaa";
  String tmp_dir = work_dir + "/" + "tmp";
  String log = work_dir + "/" + "LongProto.log";
  String cmd = "LongProto SAMPLE=human READS=#picard DATASET=1 LOGGING=REFTRACE_VARIANTS=True" 
    " X=" + range +
    " OUT_INT_HEAD=" + out_head +
    " TMP=" + tmp_dir +
    " > " + log;
  return cmd;
}
 


int main(int argc, char *argv[])
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String_Doc(ROOT, "Working directory");
  CommandArgument_String_OrDefault_Doc(TARGET, "fosmids", "The target region. Can be:\n"
				       "fosmids - all fosmid regions\n"
				       "<fosmid_id> - single fosmid region\n"
				       "<range> - range in the form chr:start-stop");
  CommandArgument_UnsignedInt_OrDefault_Doc(PARALLEL, 3, "Maximum number of jobs to run in parallel");
  CommandArgument_UnsignedInt_OrDefault_Doc(REGION_SIZE, 50, "Split genome into regions of this size (kb)");
  CommandArgument_UnsignedInt_OrDefault_Doc(REGION_OVERLAP, 10, "Minimum overlap between regions (kb)");
  EndCommandArguments;

  const uint32_t region_size = REGION_SIZE * 1000;
  const uint32_t region_overlap = REGION_OVERLAP * 1000;
  
  const String root_dir = ROOT;
  const String regions_dir = root_dir + "/" + "regions";
  
  Mkpath(root_dir);
  Mkpath(regions_dir);

  ForkManager pool(PARALLEL);
  
  vec<int> fosmid_ids = AllFosmids( );
  if (TARGET != "fosmids") 
    if ( TARGET.IsInt() )
      if (std::find(fosmid_ids.begin(), fosmid_ids.end(), TARGET.Int()) != fosmid_ids.end() )
	fosmid_ids = vec<int>(1,TARGET.Int());
      else
	FatalErr("Unknown Fosmid ID: " + TARGET);
  
  // Obtain list of regions to run
  vec<Region> region_list;
  if (TARGET == "fosmids" ||  TARGET.IsInt() ) {
    // User requested fosmid(s)

    if ( TARGET.IsInt() )
      // Single fosmid
      if (std::find(fosmid_ids.begin(), fosmid_ids.end(), TARGET.Int()) != fosmid_ids.end() )
	fosmid_ids = vec<int>(1,TARGET.Int());
      else
	FatalErr("Unknown Fosmid ID: " + TARGET);

    // for each fomsid in list
    for (int id : fosmid_ids) {
      // Get fosmid range info
      int g, rstart, rstop;
      String gid, loc;
      GetRegionInfo( ToString(id), g, rstart, rstop, gid, loc );
      String chr_name = (g == 22 ? "X": ToString(g+1) );
      // Add regions to list
      vec<Region> r = convert_range_to_region_list( chr_name, rstart, rstop, region_size, region_overlap);
      region_list.append(r);
      cout << "Fosmid: " << ToString(id) << "  Range: " << loc << "  Regions: " << ToString(r.size())
	   << " [" << ToString(std::get<1> (r[0]) ) << "-" << ToString(std::get<2> (r[r.size() -1 ]) ) << "]" << endl;
    }

  } else { 
    // User supplied range 
    String chr_name = TARGET.Before( ":" );
    int rstart = TARGET.Between( ":", "-" ).Int( );
    int rstop = TARGET.After( "-" ).Int( );
    // Add regions to list
    vec<Region> r = convert_range_to_region_list( chr_name, rstart, rstop, region_size, region_overlap);
    region_list.append(r);
    cout << "Range: " << TARGET << "  Regions: " << ToString(r.size())
	   << " [" << ToString(std::get<1> (r[0]) ) << "-" << ToString(std::get<2> (r[r.size() -1 ]) ) << "]" << endl;
  }


  // Build job list for regions
  for (Region region : region_list) {
    String chr_name;
    uint32_t region_start, region_stop;
    std::tie (chr_name, region_start, region_stop) = region;
    String work_dir = regions_dir + "/" + chr_name + "/" +  ToString(region_start) + "-" + ToString(region_stop);
    String cmd = build_command(work_dir, chr_name, region_start, region_stop);
    pool.add(cmd);
    Mkpath(work_dir);
  }
  
  // Run jobs
  int failed_jobs = pool.run();
  if (failed_jobs != 0) {
    cout << "Warning: " <<  failed_jobs << " jobs failed." << endl;
    cout << pool;
    return -1;
  }  else
    cout << "All jobs completed successfully." << endl;

}

