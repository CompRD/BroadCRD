///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// GrabReads.  Dump the reads that are in some contigs, along with their partners.
// This returns the .filt reads.

#include "Basevector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "Superb.h"
#include "paths/ReadLoc.h"

int main(int argc, char *argv[])
{
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
     CommandArgument_String_OrDefault_Doc(TIGS, "", 
          "if unspecified, process all contigs;"
          " otherwise it is one of the following: \n"
          "(a) a list of contig ids (in ParseIntSet format) or \n"
          "(b) the letter s followed by a list of scaffolds or \n"
          "(c) s<scaffold id>.<list of indices of contigs in the scaffold");
     CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
          "if specified, file extension is .READLOCS_PREFIX.readlocs "
          "instead of .readlocs");
     CommandArgument_Bool_OrDefault(PARTNERS, True);
     CommandArgument_Bool_OrDefault(FRAG_EDIT, False);
     CommandArgument_String_OrDefault(QUALB, "");
     EndCommandArguments;

     // Define directories.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

     // Parse TIGS.

     vec<int> tigs;
     {    if (TIGS == "") 
          {    int n_tigs = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY
                    + ".contigs.fastb" );
               for ( int j = 0; j < n_tigs; j++ )
                    tigs.push_back(j);    }
          else if (TIGS.Contains("s", 0)) 
          {    TIGS = TIGS.After("s");
               vec<superb> scaffolds;
               ReadSuperbs(sub_dir + "/" + ASSEMBLY + ".superb", scaffolds);
               if (TIGS.Contains(".")) 
               {    int scaffold = TIGS.Before(".").Int();
                    ForceAssertLt(scaffold, scaffolds.isize());
                    vec<int> spos;
                    ParseIntSet(TIGS.After("."), spos);
                    for (int j = 0; j < spos.isize(); j++)
                         tigs.push_back(scaffolds[scaffold].Tig(spos[j]));    }
               else 
               {    vec<int> s;
                    ParseIntSet(TIGS, s);
                    for (int i = 0; i < s.isize(); i++) 
                    {    int scaffold = s[i];
                         ForceAssertLt(scaffold, scaffolds.isize());
                         for (int j = 0; j < scaffolds[scaffold].Ntigs(); j++)
                              tigs.push_back(scaffolds[scaffold].Tig(j));
                              }    }    }
          else ParseIntSet(TIGS, tigs);    }

     // Process contigs.

     String head = sub_dir + "/" + ASSEMBLY;
     if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
     read_locs_on_disk locs_file( head, run_dir );
     vec<int> frag_ids, jump_ids, long_ids;
     for ( int it = 0; it < tigs.isize( ); it++ )
     {    vec<read_loc> locs;
          locs_file.LoadContig( tigs[it], locs );
          for ( int j = 0; j < locs.isize( ); j++ )
          {    const read_loc& rl = locs[j];
               if ( rl.Frag( ) )
               {    frag_ids.push_back( rl.ReadId( ) );
                    if (PARTNERS) frag_ids.push_back( rl.PartnerReadId( ) );    }
               if ( rl.Jump( ) )
               {    jump_ids.push_back( rl.ReadId( ) );
                    if (PARTNERS) jump_ids.push_back( rl.PartnerReadId( ) );    }
               if ( rl.LongJump( ) )
               {    long_ids.push_back( rl.ReadId( ) );
                    if (PARTNERS)
                         long_ids.push_back( rl.PartnerReadId( ) );    }    }    }
     UniqueSort(frag_ids);
     UniqueSort(jump_ids);
     UniqueSort(long_ids);
     vecbasevector frag_reads, jump_reads, long_reads;
     vecqualvector frag_quals, jump_quals, long_quals;
     if ( frag_ids.size( ) > 0 )
     {    if ( !FRAG_EDIT )
          {    frag_reads.Read( run_dir + "/frag_reads_filt_cpd.fastb", frag_ids );
               if ( QUALB != "" )
               {    frag_quals.Read( run_dir 
                         + "/frag_reads_filt_cpd.qualb", frag_ids );    }    }
          else 
          {    frag_reads.Read( run_dir + "/frag_reads_edit.fastb", frag_ids );
               if ( QUALB != "" )
               {    frag_quals.Read( run_dir + "/frag_reads_edit.qualb",
                         frag_ids );    }    }    }
     if ( jump_ids.size( ) > 0 )
     {    jump_reads.Read( run_dir + "/jump_reads_filt_cpd.fastb", jump_ids );
          if ( QUALB != "" )
          {    jump_quals.Read( run_dir + "/jump_reads_filt_cpd.qualb", 
                    jump_ids );    }    }
     if ( long_ids.size( ) > 0 )
     {    long_reads.Read( run_dir + "/long_jump_reads_filt.fastb", long_ids );
          if ( QUALB != "" )
          {    long_quals.Read( 
                    run_dir + "/long_jump_reads_filt.qualb", long_ids );    }    }
     for ( size_t i = 0; i < frag_reads.size( ); i++ )
          frag_reads[i].Print( cout, "frag_read_" + ToString( frag_ids[i] ) );
     for ( size_t i = 0; i < jump_reads.size( ); i++ )
          jump_reads[i].Print( cout, "jump_read_" + ToString( jump_ids[i] ) );
     for ( size_t i = 0; i < long_reads.size( ); i++ )
          long_reads[i].Print( cout, "long_jump_read_" + ToString( long_ids[i] ) );
     if ( QUALB != "" )
     {    vecqualvector allq(frag_quals);
          allq.Append(jump_quals);
          allq.Append(long_quals);
          allq.WriteAll(QUALB);    }    }
