///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "A simple utility to print out source information for read IDs.  Will "
  "iterate back through as much historical read-renumbering information it "
  "has (from modules which use ReadTracker).";

/* ReadTrack
 *
 * Input files required:
 *
 * <READS>.readtrack
 * 
 * Other Arguments:
 * 
 * READ_IDS
 * A set of read IDs in braces (e.g., "{2345,123132,3454363}").
 *
 * Bruce Walker
 * 16 Dec 09
 ******************************************************************************/

#include <map>
#include <set>

#include "MainTools.h"
#include "util/ReadTracker.h"
#include "util/RunCommand.h"

// Return filename portion of path (after last slash)
String basename(const String path)
{
  size_t n = path.rfind("/");
  return path.substr(n+1);
}


int main( int argc, char **argv )
{
  RunTime( );

  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(READS,
    "Track read IDS in this file (requires a <READS>.readtrack file)");
  CommandArgument_String_OrDefault_Doc(READ_IDS, "",
    "Read IDs to track {id1,id2,id3,...} or track all by default.");
  CommandArgument_Bool_OrDefault_Doc(TRACK_BACK, True,
    "Track read IDs back to their origin, and not just the most recent change.")
  EndCommandArguments;

  // Trim off file extension if one was supplied
  READS = READS.SafeBefore(".fastb");
  READS = READS.SafeBefore(".readtrack");
  if (READS.EndsWith("."))
    READS = READS.RevBefore(".");

  // Make sure there is a readtrack file
  if ( !IsRegularFile( READS + ".readtrack" ) )
    FatalErr( "Unable to find the readtracker file: " + READS  + ".readtrack" );

  // Map of parent read trackers
  typedef map<String, ReadTracker> TrackMap;
  TrackMap tracks;

  // Set of missing read trackers that could not be found on disk
  typedef set<String> MissingMap;
  MissingMap missing;

  // Load the initial read tracker
  ReadTracker rt;
  rt.Load(READS);


  // Optionally only track a subset of the read IDs
  vec<longlong> read_ids;
  uint32_t n_reads;
  bool track_all_reads = (READ_IDS == "");
  if (track_all_reads) {
    n_reads = rt.size();
  } else {
    ParseLongLongSet(READ_IDS, read_ids, False);
    n_reads = read_ids.size();
  }

  
  // Display read ID history

  for (size_t i = 0; i < n_reads; ++i) {
    
    uint32_t r = (track_all_reads ? i : read_ids[i]);
    ForceAssertLt((uint32_t)r, rt.size());
    String source = rt.GetReadSource(r);
    uint32_t index = rt.GetReadIndex(r);

    cout << basename(READS) << ":" << r << " " << basename(source) << ":" << index;

    // Track read IDs back through previous read tracker files if possible
    while (TRACK_BACK) {

      // Get next read tracker if already loaded, or else look for it on disk
      TrackMap::iterator pos = tracks.find(source);
      if (pos == tracks.end()) { 

	// Only look on disk once, if we can't find it mark it in missing.
	if (missing.find(source) == missing.end()) { 

	  // Attempt to load the next read tracker (stripping off path if required)
	  ReadTracker rt_next;
	  if ( IsRegularFile( source + ".readtrack") )
	    rt_next.Load(source);
	  else
	    rt_next.Load(basename(source));

	  // Add read tracker to list or mark it as missing
	  if (rt_next.size() == 0) {
	    missing.insert(source); 
	  } else {
	    pos = tracks.insert(TrackMap::value_type(source, rt_next)).first;
	  }
	}
      }

      // If we found another tracker then write the ID history it contains
      if (pos != tracks.end()) {
	ReadTracker& rt_next = pos->second;
	ForceAssertLt((uint32_t)index, rt_next.size());
	source = rt_next.GetReadSource(index);
	index = rt_next.GetReadIndex(index);
	cout << " " << basename(source) << ":" << index;
      } else {
	// No more read trackers available - the history ends here
	break;
      }
    }

    cout << endl;
  }

  
}


