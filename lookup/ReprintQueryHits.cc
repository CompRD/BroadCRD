/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// ReprintQueryHits: read in parseable output ALIGNS of QueryLookupTable and
// reformat.  The other arguments to this program are a subset of those for
// QueryLookupTable.

// INPUT PARAMETERS:
//
//      LOOKUP_TABLE (or L): name of lookup table file
//
//      SEQS: name of input file containing query sequence bases
//      (if the name ends with fastb, this is treated as an Arachne-style fastb
//      file; otherwise it must be in fasta format)
//
//      SEQ_NAMES: optional file of names for the query sequences
//
//      QUALS: name of an optional input file containing query sequence quality
//      scores (if the name ends with qualb, this is treated as an Arachne-style
//      qualb file; otherwise it must be in fasta format)
//
//      START, STOP: specify range [START,STOP) of query sequences to process;
//      if unspecified, all query sequences are processed (only implemented for
//      fastb and qualb files)
//
//      SEQS_TO_PROCESS: if provided, this is a file of numeric query sequence ids
//      which are to be processed.  This cannot be used with START or STOP.
//
//      SEQ_NAMES_TO_EXCLUDE: if provided, this is a file of query names which are
//      to be excluded.

// PARSABLE OUTPUT PARAMETERS:
//
//      One or more of the following may be set to True to cause alignments to be
//      printed in the given style.  The default is to have only READABLE_BRIEF
//      set to True.  See README for descriptions.
//
//      PARSEABLE
//
//      PARSEABLE_BRIEF
//
//      VISUAL
//
//      REVERSE_DISPLAY: modifies VISUAL, causing rc alignments to be printed with
//                       the query sequence forward, the reference reversed.
//
//      READABLE_BRIEF
//
//      Naming of target genome contigs (fasta records) is governed by the parameter
//      TARGET_NAMING, which may be set to one of the following values:
//
//           numeric (default): the numeric index of the fasta record, according
//                              to the order encountered by MakeLookupTable
//
//           from_file: take the part of the file name (from which the fasta record
//                      came), after the last slash, and take the part of that
//                      before the first dot, and add [n] to it, where n is the
//                      index of the fasta record in the file
//
//           from_record: the part of the fasta record header which comes after
//           the ">"
//
//      Naming of query contigs (fasta records) is governed by the parameter
//      QUERY_NAMING, which may be set to one of the following values:
//
//           numeric (default): the numeric index of the fasta record
//
//           from_record: the part of the fasta record header which comes after
//           the ">" and before the first white space
//
//           from_names_file: from a file specified by the SEQ_NAMES argument

// LOGGING PARAMETERS:
//
//      DUMP_WINDOW: if positive, for each alignment, print the query sequence
//      and the corresponding window on the target sequence, extending by
//      the given amount in both directions.  (This facilitates experimentation
//      with alternate aligners.)
//
//      NO_HEADER: if True, omit the initial header showing the name of this
//      executable, etc.

// =================================================================================


#include <strstream>

#include "Basevector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTable.h"
#include "math/Functions.h"

#define ABORT(MSG)                                  \
{    cout << MSG << "  Abort." << endl << endl;     \
     exit(1);    }

// The maximum query size is 200,000,000.  The maximum number of bases
// in the target fasta files is 4,000,000,000.  These two numbers add up
// to less than 2^32, for a reason!

const unsigned int max_query_size = 200 * 1000 * 1000;

void arachne_signal_handler( int signal_number );

const unsigned int undefined = 1000000000;

int main( int argc, char *argv[] )
{
     NoDump( );
     ArachneInterruptHandler( &arachne_signal_handler_no_ctrlc_traceback );

     BeginCommandArguments;
     CommandArgument_String(ALIGNS);
     CommandArgument_String(SEQS);
     CommandArgument_String_OrDefault(SEQ_NAMES, "");
     CommandArgument_String_OrDefault(QUALS, "");
     CommandArgument_UnsignedInt_OrDefault(START, 0);
     CommandArgument_UnsignedInt_OrDefault(STOP, undefined);
     CommandArgument_String_OrDefault(LOOKUP_TABLE, "");
     CommandArgument_String_OrDefault(L, "");
     CommandArgument_UnsignedInt_OrDefault(DUMP_WINDOW, 0);
     CommandArgument_Bool_OrDefault(PARSEABLE, False);
     CommandArgument_Bool_OrDefault(PARSEABLE_BRIEF, False);
     CommandArgument_Bool_OrDefault(VISUAL, False);
     CommandArgument_Bool_OrDefault(REVERSE_DISPLAY, False);
     CommandArgument_Bool_OrDefault(READABLE_BRIEF, True);
     CommandArgument_String_OrDefault(TARGET_NAMING, "numeric");
     CommandArgument_String_OrDefault(QUERY_NAMING, "numeric");
     CommandArgument_Bool_OrDefault(PRINT_MEMORY_USAGE, False);
     CommandArgument_String_OrDefault(SEQS_TO_PROCESS, "");
     CommandArgument_String_OrDefault(SEQ_NAMES_TO_EXCLUDE, "");
     EndCommandArguments;

     // Process abbreviations and defaults for abbreviated parameters.

     if ( LOOKUP_TABLE == "" && L == "" )
          ABORT( "You must specify LOOKUP_TABLE (or L)." );
     if ( LOOKUP_TABLE != "" && L != "" && LOOKUP_TABLE != L )
          ABORT( "You can't specify both LOOKUP_TABLE and L." );
     if ( L != "" ) LOOKUP_TABLE = L;

     // Do some sanity checks on command line arguments.

     if ( TARGET_NAMING != "numeric" && TARGET_NAMING != "from_file"
          && TARGET_NAMING != "from_record"
          && TARGET_NAMING != "from_record_quoted" )
          ABORT( "TARGET_NAMING must either be \"numeric\" or \"from_file\""
               << " or \"from_record\"" << " or \"from_record_quoted\"." );
     if ( QUERY_NAMING != "numeric" && QUERY_NAMING != "from_record"
          && QUERY_NAMING != "from_names_file" )
     {    ABORT( "QUERY_NAMING must either be \"numeric\" or \"from_record\" "
               << " or \"from_names_file\"." );    }
     if ( START != 0 || STOP != undefined )
     {    if ( !SEQS.Contains( ".fastb", -1 )
               || ( QUALS != "" && !QUALS.Contains( ".qualb", -1 ) ) )
          {    ABORT( "START and STOP may only be specified if you use "
                    << "fastb and qualb files." );    }
          if ( SEQS_TO_PROCESS != "" )
          {    ABORT( "START and STOP may not be specified if SEQS_TO_PROCESS "
                    << "is also specified." );    }    }
     if ( !IsRegularFile(SEQS) ) ABORT( "I can't find your SEQS file." );
     if ( !IsRegularFile(LOOKUP_TABLE) )
          ABORT( "I can't find your LOOKUP_TABLE file." );
     if ( SEQS.Contains( ".fastb", -1 ) )
     {    if ( QUERY_NAMING != "numeric" && QUERY_NAMING != "from_names_file" )
          {    ABORT( "If SEQS is a fastb file, QUERY_NAMING must be \"numeric\" "
                    << " or \"from_names_file\"." );    }    }
     if ( QUERY_NAMING == "from_names_file" )
     {    if ( SEQ_NAMES == "" )
          {    ABORT( "You've specified QUERY_NAMING=from_names_file but "
                    << "haven't given a value for SEQ_NAMES." );    }    }

     // Read in lists of sequences to process and to exclude, if provided.

     vec<int> to_process;
     vec<String> to_exclude_names;
     if ( SEQS_TO_PROCESS.size( ) > 0 )
     {    Ifstream( in, SEQS_TO_PROCESS );
          while(1)
          {    int n;
               in >> n;
               if ( !in ) break;
               to_process.push_back(n);    }
          UniqueSort(to_process);
          if ( to_process.size( ) == 0 )
               ABORT( "No entries found in SEQS_TO_PROCESS." )
          if ( SEQS.Contains( ".fastb", -1 ) )
          {    START = to_process.front( );
               STOP = to_process.back( ) + 1;    }    }
     if ( SEQ_NAMES_TO_EXCLUDE.size( ) > 0 )
     {    Ifstream( in, SEQ_NAMES_TO_EXCLUDE );
          while(1)
          {    static String n;
               in >> n;
               if ( !in ) break;
               to_exclude_names.push_back(n);    }
          UniqueSort(to_exclude_names);    }

     // Read in sequences.

     vecbasevector seq;
     if ( SEQS.Contains( ".fastb", -1 ) )
     {    if ( START == 0 && STOP == undefined ) seq.ReadAll(SEQS);
          else
          {    if ( STOP > MastervecFileObjectCount(SEQS) )
               {    ABORT( "You've requested query sequence number " << STOP - 1
                         << ".  There aren't that many\nsequences." );    }
               seq.ReadRange( SEQS, START, STOP, 0 );    }    }
     else FetchReads( seq, 0, SEQS );
     vec<String> seq_names;
     seq_names.reserve( seq.size( ) );
     if ( QUERY_NAMING == "numeric" )
     {    for ( size_t i = 0; i < seq.size( ); i++ )
               seq_names.push_back( ToString(i+START) );    }
     else if ( QUERY_NAMING == "from_record" )
     {    fast_ifstream in(SEQS);
          String line, read_name;
          while(1)
          {    if ( in.fail( ) ) break;
               getline( in, line );
               if ( in.fail( ) ) break;
               ForceAssert( line.size( ) > 0 && line[0] == '>' );
               read_name.resize(0);
               for ( unsigned int i = 1; i < line.size( ); i++ )
               {    if ( isspace( line[i] ) ) break;
                    read_name += line[i];    }
               seq_names.push_back(read_name);
               while(1)
               {    char c;
                    in.peek(c);
                    if ( in.fail( ) || c == '>' ) break;
                    getline( in, line );    }    }    }
     else if ( QUERY_NAMING == "from_names_file" )
     {    fast_ifstream in(SEQ_NAMES);
          String line;
          unsigned int count = 0;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( count >= START && count < STOP ) seq_names.push_back(line);
               ++count;    }    }
     if ( seq.size( ) != seq_names.size( ) )
     {    cout << "The number of query sequences (" << seq.size( ) << ") does not "
               << "agree with the number of query\nsequence names "
               << "(" << seq_names.size( ) << ").  Abort." << endl;
          exit(1);    }
     for ( size_t i = 0; i < seq.size( ); i++ )
          if ( seq[i].size( ) > max_query_size )
          {    ABORT( "One of your query sequences has length " << seq[i].size( )
                    << ".  The maximum allowed length is 200 Mb." );    }

     // Read in quality scores for sequences, if given.

     vecqualvector qual;
     if ( !IsRegularFile(QUALS) );
     else if ( QUALS.Contains( ".qualb", -1 ) )
     {    if ( START == 0 && STOP == undefined ) qual.ReadAll(QUALS);
          else qual.ReadRange( QUALS, START, STOP, 0 );    }
     else
     {    int total_seqs = 0;
          longlong total_bases = 0;
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 ) qual.Reserve( total_bases, total_seqs );
               fast_ifstream quals(QUALS);
               String line;
               qualvector q;
               Bool first = True;
               while(1)
               {    getline( quals, line );
                    if ( quals.fail( ) )
                    {    if ( pass == 1 )
                         {    total_bases += q.size( );
                              ++total_seqs;    }
                         if ( pass == 2 ) qual.push_back(q);
                         q.clear( );
                         break;    }
                    ForceAssert( line.size( ) > 0 );
                    if ( line[0] == '>' )
                    {    if ( !first )
                         {    if ( pass == 1 )
                              {    total_bases += q.size( );
                                   ++total_seqs;    }
                              if ( pass == 2 ) qual.push_back(q);    }
                         first = False;
                         q.clear( );
                         while(1)
                         {    char c;
                              quals.peek(c);
                              if ( quals.fail( ) || c == '>' ) break;
                              getline( quals, line );
                              for ( int j = 0; j < (int) line.size( ); j++ )
                                   ForceAssert( isspace(line[j])
                                        || isdigit(line[j]) );
                              istrstream i( line.c_str( ) );
                              while(1)
                              {    int n;
                                   i >> n;
                                   if ( i.fail( ) ) break;
                                   q.push_back(n);    }    }    }
                    if ( quals.fail( ) )
                    {    if ( pass == 1 )
                         {    total_bases += q.size( );
                              ++total_seqs;    }
                         if ( pass == 2 ) qual.push_back(q);
                         q.clear( );
                         break;    }    }    }    }

     // Do sanity check on consistency of seqs and quals.

     if ( qual.size( ) )
     {    if ( seq.size( ) != qual.size( ) )
          {    ABORT( "You've supplied " << seq.size( ) << " query sequences "
                    << "but " << qual.size( ) << " quality score sequences." );    }
          for ( int i = 0; i < (int) seq.size( ); i++ )
          {    if ( seq[i].size( ) != qual[i].size( ) )
               {    ABORT( "Query sequence " << i << " has length "
                         << seq[i].size( ) << ", but you've supplied "
                         << qual[i].size( ) << " quality scores." );    }    }    }

     // Read header information from lookup table.
     lookup_table look(LOOKUP_TABLE);

     // Define target names.

     vec<String> target_names( look.NContigs( ) );
     for ( unsigned int i = 0; i < look.NContigs( ); i++ )
     {    if ( TARGET_NAMING == "numeric" ) target_names[i] = ToString(i);
          else if ( TARGET_NAMING == "from_file" )
               target_names[i] = look.ContigNameBasic(i);
          else if ( TARGET_NAMING == "from_recordd" )
               target_names[i] = look.ContigName(i);
          else
          {    static String targ;
               targ = look.ContigName(i);
               targ.GlobalReplaceBy( "\"", "" );
               target_names[i] = "\"" + targ + "\"";    }    }

     // Process QueryLookupTable output.

     String line;
     fast_ifstream in(ALIGNS);
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.size( ) < 5 || line[0] != 'Q' || line[1] != 'U' ) continue;
          if ( !line.Contains( "QUERY", 0 ) ) continue;
          static look_align_plus q;
          q.ReadParseable(line);

          q.query_id -= START;
          if ( BinMember( to_exclude_names, seq_names[q.query_id] ) ) continue;

          // Fetch genome bases.

          unsigned int start = (unsigned int) -1;
          static basevector t, trc;
          if ( PARSEABLE || VISUAL )
          {    unsigned int begin = look.ContigStart( q.target_id );
               start = 0;
               if ( 80 < q.a.pos2( ) ) start = q.a.pos2( ) - 80;
               unsigned int stop = q.a.Pos2( ) + 80;
               if ( !look.CanFetchBasesFromDisk(
                    begin + start, begin + stop ) )
               {    start = q.a.pos2( );
                    stop = q.a.Pos2( );    }
               t.Setsize( stop - start );
               look.FetchBasesFromDisk( begin + start, begin + stop, t );
               trc = t;
               trc.ReverseComplement( );    }

          // Print.

          if (PARSEABLE)
               q.PrintParseable( cout, seq[ q.query_id ],
                    ( qual.size( ) ? qual[ q.query_id ]
                         : qualvector(0) ),
                    t, start, seq_names, target_names );
          if (PARSEABLE_BRIEF)
               q.PrintParseableBrief( cout, seq[ q.query_id ],
                    ( qual.size( ) ? qual[ q.query_id ]
                         : qualvector(0) ),
                    t, start, seq_names, target_names );
          if (READABLE_BRIEF)
               q.PrintReadableBrief( cout, seq[ q.query_id ],
                    ( qual.size( ) ? qual[ q.query_id ]
                         : qualvector(0) ),
                    t, start, seq_names, target_names );
          if (VISUAL)
               q.PrintVisual( cout, seq[ q.query_id ],
                    ( qual.size( ) ? qual[ q.query_id ]
                         : qualvector(0) ), t, trc, start, False, REVERSE_DISPLAY );

          if ( DUMP_WINDOW > 0 )
          {    int id = q.query_id;
               seq[id].Print( cout, "query_" + ToString(id+START) );
               unsigned int begin = look.ContigStart( q.target_id );
               static basevector w;
               unsigned int start = 0;
               if ( (int) DUMP_WINDOW < q.a.pos2( ) )
                    start = q.a.pos2( ) - DUMP_WINDOW;
               unsigned int stop = q.a.Pos2( ) + DUMP_WINDOW;
               if ( !look.CanFetchBasesFromDisk( begin + start,
                    begin + stop ) )
               {    cout << "(Bases from " << begin + start << " to "
                         << begin + stop
                         << " are unavailable on disk.  Can't dump "
                         << "target window.  Sorry.)\n";    }
               else
               {    w.Setsize( stop - start );
                    look.FetchBasesFromDisk( begin + start,
                         begin + stop, w );
                    w.Print( cout, "genome_" + target_names[ q.target_id ]
                         + "." + ToString(start) + "-"
                         + ToString(stop) );    }    }    }    }
