// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// PhysicalSupers: given alignments of reads in a project to reference genome,
// find those read pairs which are validated in the sense that both ends can be
// placed logically on the genome, roughly the right distance apart.  For each
// such pair, mark all the bases in the corresponding interval on the reference 
// genome.  Call a connected marked component a "physical super".  Compute the 
// N50 of the physical super sizes.  This gives an estimate for the theoretical
// upper limit of the supercontig N50 in an assembly.

// If INDIVIDUAL is defined, use only reads from the given individual.  (For this,
// reads.individual must exist.)

#include "MainTools.h"
#include "FastIfstream.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "ReadPairing.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTable.h"

const unsigned int UNDEFINED = 1000000000;

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Double_OrDefault(MAX_DEV, 4.5);
     CommandArgument_String(ALIGNMENT_FILES);
     CommandArgument_UnsignedInt(FIRST_CONTIG);
     CommandArgument_UnsignedInt(LAST_CONTIG);
     CommandArgument_UnsignedInt_OrDefault(INDIVIDUAL, UNDEFINED);
     EndCommandArguments;

     String run_dir = PRE + "/" + DATA + "/" + RUN;

     // Load some assembly data.

     READ( run_dir + "/reads.pairto", vec<read_pairing>, pairs );
     int nreads = MastervecFileObjectCount( run_dir + "/reads.ids" );

     vec<int> individual;
     if (INDIVIDUAL != UNDEFINED) 
          BinaryReader::readFile( run_dir + "/reads.individual", &individual );

     // Read alignments.  Note that we may have aligns_count < aligns.size( ),
     // with the difference explained by junk at the end of aligns.

     vec< vec<int> > aligns_index(nreads);
     vec<String> aligns_files = AllFilesInSource(ALIGNMENT_FILES);
     int aligns_count;
     String line;
     vec<look_align_plus> aligns;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) aligns.resize(aligns_count);
          aligns_count = 0;
          for ( int i = 0; i < aligns_files.isize( ); i++ )
          {    fast_ifstream in( aligns_files[i] );
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    if ( !line.Contains( "QUERY", 0 ) ) continue;
                    if ( pass == 2 ) 
                    {    aligns[aligns_count].ReadParseable(line);
                         int t = aligns[aligns_count].target_id;
                         if ( t < (int) FIRST_CONTIG || t > (int) LAST_CONTIG ) 
                              continue;
                         int id = aligns[aligns_count].query_id;
                         aligns_index[id].push_back(aligns_count);    }
                    ++aligns_count;    }    }    }

     // Look for validated pairs.

     vec< vec<ho_interval> > covered( LAST_CONTIG + 1 );
     for ( int i = 0; i < pairs.isize( ); i++ )
     {    int id1 = pairs[i].id1, id2 = pairs[i].id2;
          if ( INDIVIDUAL != UNDEFINED && individual[id1] != (int) INDIVIDUAL )
               continue;

          // If one end or the other isn't uniquely placed on the genome, try to
          // validate the pair anyway by looking for some consistent logical
          // full-length placement (of which there may be more than one).

          if( aligns_index[id1].size( ) > 5 ) continue;
          if( aligns_index[id2].size( ) > 5 ) continue;
          for ( int j1 = 0; j1 < aligns_index[id1].isize( ); j1++ )
          {    const look_align_plus& l1 = aligns[ aligns_index[id1][j1] ];
               if ( !l1.FullLength( ) ) continue;
               for ( int j2 = 0; j2 < aligns_index[id2].isize( ); j2++ )
               {    const look_align_plus& l2 = aligns[ aligns_index[id2][j2] ];
                    if ( !l2.FullLength( ) ) continue;
                    if ( l1.target_id != l2.target_id ) continue;
                    if ( l1.rc1 == l2.rc1 ) continue;
                    int sep;
                    if ( !l1.rc1 ) sep = l2.a.pos2( ) - l1.a.Pos2( );
                    else sep = l1.a.pos2( ) - l2.a.Pos2( );
                    float dev = float( Abs( sep - pairs[i].sep ) ) 
                         / float( pairs[i].sd );
                    if ( ( sep >= 0 && sep <= 2 * pairs[i].sep ) || dev <= MAX_DEV )
                    {    if ( !l1.rc1 )
                         {    ho_interval h( l1.a.pos2( ), l2.a.Pos2( ) );
                              covered[ l1.target_id ].push_back(h);    }
                         else
                         {    ho_interval h( l2.a.pos2( ), l1.a.Pos2( ) );
                              covered[ l1.target_id ].
                                   push_back(h);    }    }    }    }    }

     // Compute coverage.

     vec< vec<ho_interval> > cov( LAST_CONTIG + 1 );

     vec<int> physical_super_sizes;
     for ( int i = (int) FIRST_CONTIG; i <= (int) LAST_CONTIG; i++ )
     {    int n = 0;
          for ( int j = 0; j < covered[i].isize( ); j++ )
               n = Max( n, covered[i][j].Stop( ) );
          ExtractGivenCoverage( n, 1, covered[i], cov[i] );
          for ( int j = 0; j < cov[i].isize( ); j++ )
               physical_super_sizes.push_back( cov[i][j].Length( ) );    }
     Sort(physical_super_sizes);
     PRINT( N50(physical_super_sizes) );    }
