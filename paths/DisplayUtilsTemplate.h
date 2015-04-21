/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header file: DisplayUtilsTemplate.h

   Utilities to simplify the display of assembly information
   (both the final assembly and various intermediate stages)
   as HTML.
   
   @file
*/

#ifndef __INCLUDE_paths_DisplayUtilsTemplate_h
#define __INCLUDE_paths_DisplayUtilsTemplate_h

#include <sstream>
#include <algorithm>
#include "paths/DisplayUtils.h"

// Type: begend_t
// A genome placement id, together with "begin" or "end" marker.
// Represents either the start or the end of a placement
// on the genome.
typedef pair< align_id_t, Bool > begend_t;

template <class GENOME_PLACEMENT>
genome_part_pos_t BegEndPos( const begend_t& be, const vec< GENOME_PLACEMENT >& placements ) {
  const GENOME_PLACEMENT& gp = placements[ be.first ];
  return be.second ? gp.StartOnTarget() : gp.EndOnTarget()-1;
}

/**
   Functor: begend_cmp

   Used to sort begin/end markers.
*/
template <class GENOME_PLACEMENT>
struct begend_cmp {
private:
  const vec< GENOME_PLACEMENT >& placements;
  
public:
  begend_cmp( const vec<GENOME_PLACEMENT>& _genomePlacements ):
    placements( _genomePlacements ) { }
  
  bool operator() ( const begend_t& be1, const begend_t& be2 ) const {
    Bool be1isBeg = be1.second;
    Bool be2isBeg = be2.second;
    genome_part_pos_t be1pos = BegEndPos( be1, placements );
    genome_part_pos_t be2pos = BegEndPos( be2, placements );

    // first of all, sort by position on reference
    if ( be1pos < be2pos ) return True;
    if ( be1pos > be2pos ) return False;
    // if position is the same:
    //   ... put begins first, so that the begin marker of a placement is always before
    // the end marker even if they're at the same base
    if ( be1isBeg && !be2isBeg ) return True;
    if ( !be1isBeg && be2isBeg ) return False;
    // if they're both begin markers at the same position, put the marker
    // whose end is further to the right first.
    if ( be1isBeg ) {
      const GENOME_PLACEMENT& p1 = placements[ be1.first ], p2 = placements[ be2.first ];
      if ( p1.EndOnTarget() > p2.EndOnTarget() ) return True;
      if ( p1.EndOnTarget() < p2.EndOnTarget() ) return False;
    } 
    // if all else fails, sort by the placement id
    return be1.first < be2.first;
  }
};  // struct begend_cmp

/**
   Function: ShowGenomePlacements

   Given a set of possibly overlapping ranges on the reference,
   construct an HTML representation of these ranges.
*/
template <class GENOME_PLACEMENT>
html_t ShowGenomePlacements( const vec< GENOME_PLACEMENT >& placements,
			     const vec< html_t >& placement_td_contents,
			     const vec< html_attrs_t >& placement_td_attrs,
			     genome_part_id_t GPART, nbases_t genomePartSize,
			     const vec< placement_link_t > *placement_links ) {

  ForceAssertEq( placements.size(), placement_td_contents.size() );
  ForceAssertEq( placements.size(), placement_td_attrs.size() );

  genome_part_id_t gpartId = GPART;
  genome_part_id_t actualGenomePartId =
    (gpartId >= 0 ? gpartId : -gpartId-1);


  for ( int i = 0; i < placements.isize(); i++ ) {
    const GENOME_PLACEMENT& p = placements[ i ];
    //ForceAssertLe( 0, p.QueryId() );
    //ForceAssertLe( 0, p.TargetId() );
    
    //ForceAssertLe( 0, p.StartOnQuery() );
    //ForceAssertLe( p.StartOnQuery(), p.EndOnQuery() );
    
    ForceAssertLe( 0, p.StartOnTarget() );
    ForceAssertLt( p.StartOnTarget(), p.EndOnTarget() );
    ForceAssertLe( p.EndOnTarget(), genomePartSize );
    
    //ForceAssertEq( p.TargetId(), actualGenomePartId );
  }

  //
  // Make a list of records where each record denotes either the beginning
  // or the end of a range.
  //
  vec< begend_t > begEnds;
  for ( int i = 0; i < placements.isize(); i++ ) {
    //ForceAssertEq( placements[i].TargetId(), actualGenomePartId );
    begEnds.push( i, True );
    begEnds.push( i, False );
  }
  
  Sort( begEnds, begend_cmp<GENOME_PLACEMENT>( placements ) );
  
  vec< align_id_t > activePlacements;
  
  vec< vec< align_id_t > > columns;
  columns.reserve( begEnds.size() );
  
  vec< nbases_t > columnWidths;
  int maxActivePlacements = 0;
  
  typedef int table_column_num_t;

  // Local vars: first and last table column occupied by each genome_placement
  //    placementFirstCol - first column
  //    placementLastCol - last column
  vec< table_column_num_t > placementFirstCol( placements.size(), -1 ),
    placementLastCol( placements.size(), -1 );
  
  for ( int i = 0; i < begEnds.isize(); i++ ) {
    
    genome_part_pos_t lastPos = -1;
    Bool gotColWidth = False;
    int j; 
    for ( j = i; j < begEnds.isize(); j++ ) {
      
      const begend_t& be = begEnds[j];
      align_id_t gpid = be.first;
      Bool isBeg = be.second;
      
      genome_part_pos_t pos = BegEndPos( be, placements );
      
      if ( j > i  &&  pos > lastPos ) {
	columnWidths.push_back( pos - lastPos );
	gotColWidth = True;
	break;
      }
      lastPos = pos;
      
      if ( isBeg ) {
	// find the first free spot in the placements array
	int freeSpot = -1;
	for ( int j = 0; j < activePlacements.isize() && freeSpot==-1; j++ )
	  if ( activePlacements[j] == -1 )
	    freeSpot = j;
	
	if ( freeSpot == -1 ) {
	  activePlacements.resize( activePlacements.size() + 1 );
	  freeSpot = activePlacements.size()-1;
	}
	
	activePlacements[ freeSpot ] = gpid;
	placementFirstCol[ gpid ] = columns.isize();
      } else {
	replace( activePlacements.begin(), activePlacements.end(), gpid, -1 );
	ForceAssertGe( placementFirstCol[ gpid ], 0 );
	ForceAssertGt( columns.isize(), 0 );
	placementLastCol[ gpid ] = columns.isize()-1;
      }
    }  // for all begEnds at the same pos
    
    if ( !gotColWidth ) {
      for ( int k = 0; k < activePlacements.isize(); k++ )
	ForceAssertEq( activePlacements[k], -1 );
    } else {
      columns.push_back( activePlacements );
      update_max( maxActivePlacements, activePlacements.isize() );
    }
    
    i = j - 1;
  }  // for each begend
  
  ForceAssertEq( columnWidths.size(), columns.size() );

  // check that we have recorded the start and end column of each placement
  for ( align_id_t i = 0; i < placements.isize(); i++ ) {
    ForceAssertGe( placementFirstCol[ i ], 0 );
    ForceAssertGe( placementLastCol[ i ], 0 );
    ForceAssertLe( placementFirstCol[ i ], placementLastCol[ i ] );
  }
  ostringstream cleg;

  genome_part_pos_t columnStart = 1;
  
  cleg << "<table border=\"1\">" << endl;
  
  cleg << "<tr>";
  cleg << "<th valign=\"middle\" rowspan=\"" << (maxActivePlacements + 1) << "\">" <<
    ( gpartId < 0 ? "-" : "" ) << (gpartId >= 0 ? gpartId : -gpartId-1) << "</th>" << endl;
  
  nbases_t leftMargin = BegEndPos( begEnds.front(), placements );
  if ( leftMargin > 0 ) {
    cleg << "<th title=\"1 - " << ToStringAddCommas( leftMargin ) << "\" >" <<
      ToStringAddCommas( leftMargin ) << "</th>";
    columnStart += leftMargin;
  }
  
  for ( int col = 0; col < columnWidths.isize(); col++ ) {
    cleg << "<th title=\"" << ToStringAddCommas( columnStart ) << " - "
	 << ToStringAddCommas( columnStart + columnWidths[col] - 1 )  
	 << "\">" << ToStringAddCommas( columnWidths[col]) << "</th>";
    columnStart += columnWidths[ col ];
  }
  
  nbases_t rightMargin = genomePartSize - BegEndPos( begEnds.back(),  placements ) - 1;
  if ( rightMargin > 0 ) {
    cleg << "<th title=\"" << ToStringAddCommas( columnStart ) << " - "
	 << ToStringAddCommas( columnStart + rightMargin - 1) << "\" >" <<
      ToStringAddCommas( rightMargin ) << "</th>";
  }
  
  cleg << "\n</tr>" << endl;
  
  for ( int row = 0; row < maxActivePlacements; row++ ) {
    cleg << "<tr align=\"center\" >";
    if ( leftMargin > 0 )
      cleg << "<td></td>" << endl;
    for ( int col = 0; col < columns.isize(); col++ )
      if ( row >= columns[col].isize() || columns[col][row] == -1 )
	cleg << "<td></td>" << endl;
      else {
	align_id_t gpid = columns[col][row];
	cleg << "<td " + placement_td_attrs[ gpid ] + ">" + placement_td_contents[ gpid ] + "</td>";
      }
    if ( rightMargin > 0 )
      cleg << "<td></td>" << endl;
    cleg << "</tr>" << endl;
  }

  if ( placement_links ) {

    // Local type: placement_link_id_t
    // Id number of a link between two genome placements; index of the link
    // in *placement_links.
    typedef int placement_link_id_t;
    typedef triple< table_column_num_t, table_column_num_t, placement_link_id_t > link_in_row_t;

    // Local var: row2links
    // For each row, (firstCol,lastCol,linkId) of links placed on that row
    vec< vec< link_in_row_t > > row2links;
  
    for ( placement_link_id_t linkId = 0; linkId < placement_links->isize(); linkId++ ) {
      align_id_t gpid1 = (*placement_links)[ linkId ].first;
      align_id_t gpid2 = (*placement_links)[ linkId ].second;

      table_column_num_t fromCol = min( placementFirstCol[ gpid1 ], placementFirstCol[ gpid2 ] );
      table_column_num_t toCol = max( placementLastCol[ gpid1 ], placementLastCol[ gpid2 ] );

      int row = 0;
      vec< link_in_row_t >::iterator newLinkPos;
      link_in_row_t newLink( fromCol, toCol, linkId );


      for ( ; row < row2links.isize(); row++ ) {
	vec< link_in_row_t >& linksThisRow = row2links[ row ];

	newLinkPos = lower_bound( linksThisRow.begin(), linksThisRow.end(), newLink );
	if ( ( newLinkPos == linksThisRow.begin()  ||
	       (newLinkPos-1)->second < newLink.first ) &&
	     ( newLinkPos == linksThisRow.end()  ||
	       newLink.second < newLinkPos->first ) )
	  break;
      }
      if ( row == row2links.isize() ) {
	row2links.resize( row2links.isize() + 1 );
	newLinkPos = row2links.back().end();
      }
      cout << " putting link " << (*placement_links)[newLink.third] << " in row " << row << endl;
      row2links[ row ].insert( newLinkPos, newLink );
    }  // for each link

    cleg << "<tr><td colspan=\"" <<  columnWidths.size() + ( leftMargin ? 1 : 0 ) + ( rightMargin ? 1 : 0 ) + 1 << "\"><hr></td></tr>\n";

    for ( int row = 0; row < row2links.isize(); row++ ) {
      cleg << "<tr>";
      const vec< link_in_row_t >& linksThisRow = row2links[ row ];
      ForceAssert( row2links[ row ].nonempty() );
      
      cleg << "<td></td>";
      if ( leftMargin )
	cleg << "<td></td>";
      
      if ( linksThisRow.front().first > 0 )
	cleg << "<td colspan=\"" << linksThisRow.front().first << "\"></td>";
      for ( int j = 0; j < linksThisRow.isize(); j ++ ) {
	if ( j > 0 ) {
	  int colsToSkip = linksThisRow[j].first - linksThisRow[j-1].second - 1;
	  if ( colsToSkip > 0 )
	    cleg << "<td colspan=\"" << colsToSkip << "\"></td>";
	}
	cleg << "<td colspan=\"" << linksThisRow[j].second - linksThisRow[j].first + 1 << "\">" <<
	  (*placement_links)[ linksThisRow[j].third ].third << "</td>";
      }
      if ( linksThisRow.back().second < columnWidths.isize() - 1 )
	cleg << "<td colspan=\"" << columnWidths.isize() - 1 - linksThisRow.back().second << "\"></td>";
      
      if ( rightMargin )
	cleg << "<td></td>";
      
      cleg << "</tr>\n";
    }

    
  }  // if placement links were given
  
  cleg << "</table>" << endl;

  return cleg.str();
}

#define DEFINE_DISPLAY_UTILS(T) \
  template html_t ShowGenomePlacements<T>( const vec< T >& placements,   \
      const vec< html_t >& placement_td_contents,                                       \
			     const vec< html_attrs_t >& placement_td_attrs,              \
			     genome_part_id_t GPART, nbases_t genomePartSize,        \
		 const vec< placement_link_t > *placement_links ); \
  template \
   genome_part_pos_t BegEndPos<T>( const begend_t& be, const vec< T >& placements )

#endif
// #ifndef __INCLUDE_paths_DisplayUtilsTemplate_h
