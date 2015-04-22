// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "tiled/ContigInformator.h"


//  Calculate the read coverage of a contig in consecutive, non-overlapping 
//  'window'-sized regions.  Display window (in bases along the contig) and coverage.
//  coverage calculated by #bases-in-window/window size
//
void ContigCoverage( vec<read_location> &vecLocs,
		     int contig,
		     int window )
{

     pair< vec<read_location>::iterator, vec<read_location>::iterator > range;

     read_location target;
     target.SetContig(contig);
   
     range = equal_range( vecLocs.begin(), 
			  vecLocs.end(),
			  target, 
			  order_read_locations_by_contig() );
     
     int contigLen( range.first->LengthOfContig() );
      
     int num_windows( ( (int) contigLen/window) + 1 ); 
     
     cout << "Contig " << contig << " has length " << contigLen 
	  << " and " << num_windows << " windows." << endl; 

     int start(0), end(window);
     while ( end < contigLen+window )
     {
       int num_bases(0);

       vec<read_location>::iterator rlIter = range.first;
       for ( ; rlIter != range.second; ++rlIter )
       {
	 if ( rlIter->StopOnContig() < start )
	   continue;

	 if ( rlIter->StartOnContig() > end )
	   break;

	 int start_base = Max( start, rlIter->StartOnContig() );
	 int end_base = Min( end, rlIter->StopOnContig() );

	 num_bases += (end_base-start_base);

       }

       cout << start <<"-" << end 
	    << " " << (float) num_bases/window << endl;

       start = end; 
       end += window;
     }
        
}


// For UNIQUELY PLACED reads only,
// print info about a read if its partner should be in the same contig but isn't
// info includes (for both the read and its partner): 
//          placement on human
//          read location from assembly
//          separation, stdev from pairs file
//
// a read partner should be in the same contig as a read if the distance from the
// read to the end of the contig (orientation dependent) is > sep + 4*std for the 
// insert.
void ReadAndPartnerInfo( vec<read_pairing> &pairs,
			 vec<int> &pairs_index,
			 vec<read_location> &locs,
			 vec<int> &locs_by_id,
			 vec<genome_pos> &uniq_human,
			 vec<int> &individual,
			 vecString &read_names,
			 vec<int> &read_lens,
			 int contig )
{

  
  pair< vec<read_location>::iterator, vec<read_location>::iterator > range;
  read_location target;
  target.SetContig(contig);
  
  range = equal_range( locs.begin(), 
		       locs.end(),
		       target, 
		       order_read_locations_by_contig() );
     
  int contigLen( range.first->LengthOfContig() );
  cout << "Contig " << contig << " (length=" << contigLen << ") "
       << " has " << distance(range.first,range.second) << " reads\n" << endl;

  // for each read on the contig, find out if the partner is NOT placed on the 
  // contig but should be.  Also mark unpaired reads.
  for ( ; range.first != range.second; ++range.first )
    {
      int loc_indx = distance( locs.begin(), range.first );
 
      int rdid = locs[loc_indx].ReadId();

      // sanity check for correct individual (Clint)
      if ( individual[rdid] != 0 )
      {
	cout << "Read: " << read_names[rdid] << " ("
	     << rdid << ") is not from Clint. " << endl;
	continue;
      }

      int partner_indx = pairs_index[rdid];
       
       if ( partner_indx < 0 )
	 continue;
       
       int rd_partner = pairs[partner_indx].Partner(rdid);
       int should_be_in_contig = IsPartnerPlaced( locs[loc_indx], 
						  locs_by_id,
						  pairs, 
						  pairs_index );

       // partner is not in contig, but this makes sense from where the read is 
       // in the contig and the insert info, so continue
       if ( should_be_in_contig == 1 )
	 continue;


       // at this point, we've got a read whose partner should be in the contig
       // but isn't.
       // 
       // if they are both uniquely placed, print out
       if ( uniq_human[rdid].t < 0 || uniq_human[rd_partner].t < 0 )
	 continue;
       
       read_pairing this_pair = pairs[ pairs_index[ rdid ] ];

       cout << "Partner to " << read_names[rdid] << " ("
	    << rdid << ") - " << read_names[rd_partner] << " (" 
	    << rd_partner << ") - is not in contig " << contig 
	    << " but should be." << endl;
       cout << locs[loc_indx];
       
       int partner_loc_indx = locs_by_id[ rd_partner ];
       if ( partner_loc_indx < 0 )
	 cout << "rd "<< rd_partner << " (len="
	      << read_lens[rd_partner] << ") --> not in assembly."<< endl;
       else
	 cout << locs[ locs_by_id[ rd_partner ] ];

       cout << "Sep/std: " << this_pair.sep <<"/" << this_pair.sd << endl;

       cout << "Placement (both reads uniquely placed) on human:" << endl;
       cout << "\t" << read_names[rdid] << ": "<< uniq_human[rdid].pos2 << " " 
	    << uniq_human[rdid].t << endl;

       cout << "\t" << read_names[rd_partner] << ": "<< uniq_human[rd_partner].pos2
	    << " " << uniq_human[rd_partner].t << endl;

       cout << "\n";
       
    }
}




// return -1 if partner not placed (error)
// return 0 if partner not placed in same contig but should be
// return 1 if partner placed in same contig
//			 
int IsPartnerPlaced( read_location &loc,
		     vec<int> &locs_by_id,
		     vec<read_pairing> &pairs,
		     vec<int> &pairs_index )
{
	
  int rdid = loc.ReadId();
  int partner_indx = pairs_index[rdid];
  
  if ( partner_indx < 0 )
  {
    cout << "Partner not placed." << endl;
    return -1;
  }

  int contigLen = loc.LengthOfContig();

  
  //  decision:  should the partner be in this contig?
  int start_on_contig = loc.StartOnContig();
  int distance_to_contig_end(-1);
  if ( loc.OrientationOnContig() == ForwardOr )
  {
    distance_to_contig_end = contigLen-loc.StopOnContig();
  }
  else
  {
    distance_to_contig_end = loc.StartOnContig();
  }

  read_pairing &this_pair = pairs[partner_indx];
  float max_stretch = this_pair.sep + 4.0*this_pair.sd;

  int rd_partner = pairs[partner_indx].Partner(rdid);       
  int rd_partner_loc_indx = locs_by_id[ rd_partner ];
  
  //  partner not in this contig...
  if ( rd_partner_loc_indx < 0 )
  {
    // ...but it should be
    if ( distance_to_contig_end > max_stretch )
      return 0;
  }

  return 1;
}





