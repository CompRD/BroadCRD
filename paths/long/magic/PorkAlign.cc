///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "FetchReads.h"
#include "paths/long/CreateGenome.h"
#include "PrintAlignment.h"
#include "util/Fastq.h"

int main( int argc, char* argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault(FASTA,"");
     CommandArgument_String_OrDefault(FASTQ,"");
     CommandArgument_Bool_OrDefault(RC, False);
     CommandArgument_Bool_OrDefault(VISUAL,False);
     CommandArgument_String_OrDefault_Doc(OUT_HEAD, "",
             "optionally write .align and .reference");
     EndCommandArguments;

     // read the read
     vecbasevector Q;
     if ( FASTA != "" ) {
         FetchReads(Q, 0, FASTA);
         if ( Q.size() == 0 ) FatalErr("read data didn't load");
         else if ( Q.size() > 1 ) cout << "warning: only first read will be used" << endl;
     } else if ( FASTQ != "" ) {
         vecqualvector dummy;
         fastq::ReadFastq(FASTQ,Q,dummy);
     } else
         FatalErr("must specify FASTA or FASTQ");

     basevector q(Q[0]);

     if ( RC ) q.ReverseComplement();


     // read the genome
     vecbasevector GG( "/wga/dev/references/Escherichia_coli/genome.fastb" );
     for ( int g = 0; g < (int) GG.size( ); g++ )
     {    basevector x;
          for ( int j = 0; j < GG[g].isize( ); j++ )
          {    if ( j >= 5 && GG[g][j] == GG[g][j-1] && GG[g][j] == GG[g][j-2]
                    && GG[g][j] == GG[g][j-3] && GG[g][j] == GG[g][j-4]
                    && GG[g][j] == GG[g][j-5] )
               {    continue;     }
               x.push_back( GG[g][j] );    }
          GG[g] = x;    }

     vecbasevector GGR(GG);
     for ( int g = 0; g < (int) GG.size( ); g++ )
          GGR[g].ReverseComplement( );
     GG.Append(GGR);
     const int K = 12;
     VecIntPairVec Glocs;
     CreateGlocs( GG, K, Glocs );



     // do the alignment
     vec<look_align> aligns;
     int BW_ADD = 300;
     const int MIN_CLUSTER = 20;
     const int MAX_OFFSET_DIFF = 1000;
     const int MISMATCH_PENALTY = 2;
     const int GAP_PENALTY = 3;
     Bool FW_ONLY = True;
     ClusterAligner( q, GG, K, Glocs, aligns, FW_ONLY, BW_ADD, MIN_CLUSTER,
          MAX_OFFSET_DIFF, MISMATCH_PENALTY, GAP_PENALTY );
     if ( aligns.nonempty( ) )
     {    cout << "LOC: " << " at " << aligns[0].target_id
               << "." << aligns[0].a.pos2( ) << "-" << aligns[0].a.Pos2( )
               << endl;    }
     if (VISUAL)
     {    for ( int i = 0; i < aligns.isize( ); i++ )
          {    cout << "\nalignment " << i+1 << endl;
               int g = aligns[i].target_id;
               const align& a = aligns[i].a;
               cout << "at " << g << "." << a.pos2( ) << "-" << a.Pos2( ) << endl;
               cout << "nhits=" << aligns[i].nhits << ", subs="
                       << aligns[i].Mutations()
                       << ", indels=" << aligns[i].Indels() << endl;
               PrintVisualAlignment( True, cout, q, GG[g], a );    }    }

     vec<Bool> to_remove( aligns.size( ), False );
     for ( int i = 0; i < aligns.isize( ); i++ )
          if ( !aligns[i].FullLength( ) ) to_remove[i] = True;
     EraseIf( aligns, to_remove );

     PRINT_TO( cout, aligns.size( ) );
     if ( aligns.empty( ) )
     {    cout << "0% of bases in matching 10-mers" << endl;
          return 0;
     }

     if ( OUT_HEAD != "" ) {
         if ( aligns.size() > 1 )
             cout << "WARNING: will only use first of " << aligns.size()
             << " alignments" << endl;
         basevector ref;
         align a = aligns[0].a;
         int g = aligns[0].target_id;
         ForceAssertGe(a.pos2(), 0);
         ref.SetToSubOf(GG[g],a.pos2(),a.extent2());
         Ofstream(of, OUT_HEAD+".reference");
         String id = ToString(g) + "." + ToString(a.pos2());
         ref.PrintCol(of, id, 80);

         a.Setpos2(0);           // we wrote out just a portion of the reference
         PRINT3(a.pos1(), a.pos2(), a.Nblocks());
         BinaryWriter::writeFile(OUT_HEAD+".align", a);
     }

     return 0;
}
