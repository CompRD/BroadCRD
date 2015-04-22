// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research


#ifndef VALIDATE_CHIMP_ON_HUMAN_H
#define VALIDATE_CHIMP_ON_HUMAN_H

class contig_on_human {

     public:

     int t;
     int pos;
     int Pos;
     Bool human_placement_inconsistent;
     Bool no_unique_placements;
     Bool excessive_coverage;
     Bool haplotype_inconsistent;
     Bool link_inconsistent;
     Bool overlapping_contigs;

     contig_on_human( ) 
     {    t = -1;
          pos = -1;
	  Pos = -1;
	  human_placement_inconsistent = False;
	  no_unique_placements = False;
	  excessive_coverage = False;
          haplotype_inconsistent = False;
          link_inconsistent = False;
          overlapping_contigs = False;     }

     contig_on_human( int t_arg, int pos_arg, int Pos_arg )
       : t(t_arg), pos(pos_arg), Pos(Pos_arg) 
     {    human_placement_inconsistent = False;
          no_unique_placements = False;
          excessive_coverage = False;
          haplotype_inconsistent = False;
	  link_inconsistent = False;
          overlapping_contigs = False;     }

     Bool has_problem( ) const
     {    return human_placement_inconsistent ||
	         no_unique_placements || 
	         excessive_coverage || 
	         haplotype_inconsistent ||
	         link_inconsistent ||
                 overlapping_contigs;    } 

     friend istream& operator>>( istream& in, contig_on_human& x )
     {    in >> x.t >> x.pos >> x.Pos;
          int y;
	  in >> y;
	  x.human_placement_inconsistent = y ? True : False;   
	  in >> y;
	  x.no_unique_placements = y ? True : False;
	  in >> y;
	  x.excessive_coverage = y ? True : False;
	  in >> y;
	  x.haplotype_inconsistent = y ? True : False;
	  in >> y;
	  x.link_inconsistent = y ? True : False;
	  in >> y;
	  x.overlapping_contigs = y ? True : False;
          return in;     }

     friend ostream& operator<<( ostream& out, const contig_on_human& x )
     {    return out << x.t << " " << x.pos << " " << x.Pos
		     << " " << (int) x.human_placement_inconsistent
		     << " " << (int) x.no_unique_placements
		     << " " << (int) x.excessive_coverage
		     << " " << (int) x.haplotype_inconsistent
		     << " " << (int) x.link_inconsistent
		     << " " << (int) x.overlapping_contigs << "\n";    }

};



class contig_on_human_old {

     public:

     int t;
     int pos;
     int Pos;
     Bool human_placement_inconsistent;
     Bool no_unique_placements;
     Bool excessive_coverage;
     Bool haplotype_inconsistent;
     Bool link_inconsistent;

     contig_on_human_old( ) 
     {    t = -1;
          pos = -1;
	  Pos = -1;
	  human_placement_inconsistent = False;
	  no_unique_placements = False;
	  excessive_coverage = False;
          haplotype_inconsistent = False;
          link_inconsistent = False;     }

     contig_on_human_old( int t_arg, int pos_arg, int Pos_arg )
       : t(t_arg), pos(pos_arg), Pos(Pos_arg) 
     {    human_placement_inconsistent = False;
          no_unique_placements = False;
          excessive_coverage = False;
          haplotype_inconsistent = False;
	  link_inconsistent = False;     }

     Bool has_problem( ) const
     {    return human_placement_inconsistent ||
	         no_unique_placements || 
	         excessive_coverage || 
	         haplotype_inconsistent ||
	         link_inconsistent;     } 

     friend istream& operator>>( istream& in, contig_on_human_old& x )
     {    in >> x.t >> x.pos >> x.Pos;
          int y;
	  in >> y;
	  x.human_placement_inconsistent = y ? True : False;   
	  in >> y;
	  x.no_unique_placements = y ? True : False;
	  in >> y;
	  x.excessive_coverage = y ? True : False;
	  in >> y;
	  x.haplotype_inconsistent = y ? True : False;
	  in >> y;
	  x.link_inconsistent = y ? True : False;
          return in;     }

     friend ostream& operator<<( ostream& out, const contig_on_human_old& x )
     {    return out << x.t << " " << x.pos << " " << x.Pos
		     << " " << (int) x.human_placement_inconsistent
		     << " " << (int) x.no_unique_placements
		     << " " << (int) x.excessive_coverage
		     << " " << (int) x.haplotype_inconsistent
		     << " " << (int) x.link_inconsistent << "\n";    }

};


#endif
