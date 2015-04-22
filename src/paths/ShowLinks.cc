///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ShowLinks.  Display the links from scaffold S to other scaffolds.  The argument
// MAX_OVERLAP determines how large the implies overlap between S and the other
// scaffolds can be - in order to display the links.

#include "Equiv.h"
#include "MainTools.h"
#include "Superb.h"
#include "math/Functions.h"
#include "paths/ReadLoc.h"
#include "system/ParsedArgs.h"

class Link {
     public:

     Link( ) { }
     Link( const int m1, const int s1, const int m2, const int s2, const int sep, 
          const int dev, const int read_class, const int dist_to_end1, 
          const int dist_to_end2, const Bool rc1, const Bool rc2 ) : m1(m1), s1(s1),
          m2(m2), s2(s2), sep(sep), dev(dev), read_class(read_class), 
          dist_to_end1(dist_to_end1), dist_to_end2(dist_to_end2), rc1(rc1), 
          rc2(rc2) { }

     friend Bool operator<( const Link& l1, const Link& l2 )
     {    return l1.s2 < l2.s2;     }

     int m1;
     int s1;
     int m2;
     int s2;
     int sep;
     int dev;
     int read_class;
     int dist_to_end1;
     int dist_to_end2;
     Bool rc1;
     Bool rc2;

};

int main(int argc, char *argv[])
{
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
     CommandArgument_Int(S);
     CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
          "if specified, file extension is .READLOCS_PREFIX.readlocs "
          "instead of .readlocs");
     CommandArgument_Int_OrDefault(MAX_OVERLAP, 10000);
     EndCommandArguments;

     // Define directories.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

     // Load scaffolds.

     int ntigs 
          = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY + ".contigs.fastb" );
     vec<superb> scaffolds;
     ReadSuperbs( sub_dir + "/" + ASSEMBLY + ".superb", scaffolds );
     vec<int> to_super( ntigs, -1 ), to_super_pos( ntigs, -1 );
     vec<int> to_super_posr( ntigs, - 1 );
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    int n = scaffolds[i].Ntigs( );
          for ( int j = 0; j < n; j++ )
          {    to_super[ scaffolds[i].Tig(j) ] = i;
               to_super_pos[ scaffolds[i].Tig(j) ] = j;   
               to_super_posr[ scaffolds[i].Tig(j) ] = n - j - 1;    }    }

     // Process contigs.

     int s1 = S;
     const superb& S1 = scaffolds[s1];
     String head = sub_dir + "/" + ASSEMBLY;
     if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
     read_locs_on_disk locs_file( head, run_dir );
     vec<Link> Links;
     for ( int p1 = 0; p1 < scaffolds[s1].Ntigs( ); p1++ )
     {    int m1 = scaffolds[s1].Tig(p1);
          vec<read_loc> locs;
          locs_file.LoadContig( m1, locs );
          for ( int j = 0; j < locs.isize( ); j++ )
          {    const read_loc& rl = locs[j];
               if ( !rl.PartnerPlaced( ) ) continue;
               int x1 = rl.Start( );
               int y1 = scaffolds[s1].SubSuperLength( 0, p1 - 1 ) + x1;
               int m2 = rl.PartnerContigId( );
               int s2 = to_super[m2], p2 = to_super_pos[m2];
	       if ( s2 < 0 ) continue;
               const superb& S2 = scaffolds[s2];
               if ( s2 == s1 ) continue;
               int x2 = rl.PartnerStart( );
               int y2 = scaffolds[s2].SubSuperLength( 0, p2 - 1 ) + x2;

               // For now just consider 1/2 of the cases:

               if ( rl.Fw( ) && rl.PartnerRc( ) )
               {    int dist_to_end1 = S1.Len(p1) - rl.Stop( ) 
                         + S1.SubSuperLength( p1 + 1, S1.Ntigs( ) - 1 );
                    int dist_to_end2 = rl.PartnerStart( ) 
                         + S2.SubSuperLength( 0, p2 - 1 );
                    int sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
                    int dev = rl.Dev( );
                    if ( sep < -MAX_OVERLAP ) continue;
                    Links.push( m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
                         dist_to_end1, dist_to_end2, False, False );    }
               if ( rl.Fw( ) && rl.PartnerFw( ) )
               {    int dist_to_end1 = S1.Len(p1) - rl.Stop( ) 
                         + S1.SubSuperLength( p1 + 1, S1.Ntigs( ) - 1 );
                    int dist_to_end2 = S2.Len(p2) - rl.PartnerStop( ) 
                         + S2.SubSuperLength( p2 + 1, S2.Ntigs( ) - 1 );
                    int sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
                    int dev = rl.Dev( );
                    if ( sep < -MAX_OVERLAP ) continue;
                    Links.push( m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
                         dist_to_end1, dist_to_end2, False, True );    }
               if ( rl.Rc( ) && rl.PartnerRc( ) )
               {    int dist_to_end1 = rl.Start( ) + S1.SubSuperLength( 0, p1 - 1 );
                    int dist_to_end2 = rl.PartnerStart( ) 
                         + S2.SubSuperLength( 0, p2 - 1 );
                    int sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
                    int dev = rl.Dev( );
                    if ( sep < -MAX_OVERLAP ) continue;
                    Links.push( m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
                         dist_to_end1, dist_to_end2, True, False );    }    }    }

     Sort(Links);
     for ( int i = 0; i < Links.isize( ); i++ )
     {    const Link& L = Links[i];
          int j;
          for ( j = i + 1; j < Links.isize( ); j++ )
               if ( Links[j].s2 != Links[i].s2 ) break;
          equiv_rel e( j - i );
          for ( int k1 = i; k1 < j; k1++ )
          {    for ( int k2 = i + 1; k2 < j; k2++ )
               {    int delta = 
                         Abs( Links[k1].dist_to_end1 - Links[k2].dist_to_end1 )
                         + Abs( Links[k1].dist_to_end2 - Links[k2].dist_to_end2 );
                    int min_delta = 10;
                    if ( delta < min_delta && Links[k1].rc1 == Links[k2].rc1 
                         && Links[k1].rc2 == Links[k2].rc2 ) 
                    {    e.Join( k1 - i, k2 - i );    }    }    }
          vec<int> reps;
          e.OrbitReps(reps);
          if ( reps.size( ) > 1 )
          {    cout << "\n";
               for ( int r = 0; r < reps.isize( ); r++ )
               {    int k = i + reps[r];
                    const Link& L = Links[k];
                    if ( L.read_class == 0 ) cout << "frag ";
                    if ( L.read_class == 1 ) cout << "jump ";
                    if ( L.read_class == 2 ) cout << "long ";
                    int s2 = L.s2, sep = L.sep, dev = L.dev;
                    int m1 = L.m1, m2 = L.m2;
                    cout << "s" << s1 << ( L.rc1 ? ".rc" : ".fw" ) 
                         << "[" << m1 << "] --> ";
                    cout << "s" << s2 << ( L.rc2 ? ".rc" : ".fw" ) 
                         << "[" << m2 << "], ";
                    int dist_to_end1 = L.dist_to_end1;
                    int dist_to_end2 = L.dist_to_end2;
                    String dist_to_end = ToString(dist_to_end1) + "/"
                         + ToString(dist_to_end2);
                    cout << "sep: " << sep << " +/- " << dev << ", ->end: "
                         << dist_to_end << "\n";    }    }
          i = j - 1;    }    }
