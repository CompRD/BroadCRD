/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** 
   Program: EvaluateStrongKmers

   Evaluate how well our definition of <strong kmers> corresponds to the true
   <genomic kmers>.  Strong kmers are kmers that we think are genomic; genomic
   kmers are kmers that actually occur in the reference genome.

   Strong kmers are found by the program <FindStrongKmers>.

   \file
*/

#include "MainTools.h"

#include "math/Functions.h"
#include "system/System.h"

#include "ParseSet.h"
#include "kmers/KmerRecord.h"
#include "kmers/KmerShape.h"
#include "Qualvector.h"
#include "feudal/BinaryStream.h"
#include "kmer_freq/KmerShortMap.h"
#include "kmer_freq/WriteKmerFrequencies.h"
#include "paths/KmerPathInterval.h"

template <int K>
struct kmerOnlyLess
  : public binary_function<kmer_with_count<K>,kmer_with_count<K>,bool>
{
  bool operator() ( const kmer_with_count<K>& lhs, const kmer_with_count<K>& rhs ) const {
    return lt_kmer( lhs, rhs );
  }
};

/**
   Local function: Compare

   Parameters:

      file1 - <KmerShortMap> for the reference genome; this tells in particular
          which kmers are truly genomic.
      file2 - <KmerShortMap> for the strong kmers; this tells which kmers we
          think are genomic (as identified by the program <FindStrongKmers>).
      file3 - <KmerShortMap> for the original (unedited) reads: this tells which kmers
          occur in the original reads.  the reason we need this: if there are genomic
	  kmers that do not occur even once in the original reads, then we should
	  be concerned.
*/
template <class KSHAPE> 
void Compare( filename_t file1, filename_t file2, filename_t file3, 
              Bool evaluateReads, Bool showMistakes, int MAX_MISTAKES_TO_SHOW,
	      Bool LIST_KMERS, filename_t readsFile, filename_t qualsFile,
	      filename_t trustedFile, dirname_t run_dir ) {

  const nbases_t K = KSHAPE::KSIZE;
  vec< kmer_with_count<K> > table1, table2;
  PRINT2( file1, file2 );
  BinaryReader::readFile( file1, &table1 );
  BinaryReader::readFile( file2, &table2 );

  PRINT2( table1.size(), table2.size() );
  if ( !is_sorted( table1.begin(), table1.end() ) ) {
    UniqueSort( table1 );
    cout << " sorted table1!" << endl;
  } else {
    cout << " table1 already sorted." << endl;
  }
  if ( !is_sorted( table2.begin(), table2.end() ) ) {
    UniqueSort( table2 );
    cout << " sorted table 2!" << endl;
  } else {
    cout << " table2 already sorted." << endl;
  }
  
  KmerShortMap* table3 = 0;
  if (evaluateReads)
    table3 = new KmerShortMap(KSHAPE::getId(), file3);

  typedef BinaryIteratingWriter< vec<kmer_with_count<K> > > BIW;
  BIW* nonGenomicWriter = 0;
  BIW* withoutNonGenomicWriter = 0;
  BIW* withMissingWriter = 0;

  if (LIST_KMERS) {
    nonGenomicWriter = new BIW( file2 +".nonGenomicList" );
    withoutNonGenomicWriter = new BIW( file2 + ".withoutNonGenomicList" );
    withMissingWriter = new BIW( file2 + ".withMissingList" );
  }
  vec< kmer_with_count<K> > falsePositives;
  vec< kmer_with_count<K> > falseNegatives;

  typename vec< kmer_with_count<K> >::const_iterator iKmer1, iKmer2;
  iKmer1 = table1.begin();
  iKmer2 = table2.begin();

  const int cap = 1000;

  vec<int> missingByGc(K+1,0);
  vec<int> nonGenomicByGc(K+1,0);
  vec<int> sharedByGc(K+1,0);
  vec<int> strongByGc(K+1,0);
  vec<int> genomicByGc(K+1,0);

  vec<int> inReadsByGc(K+1,0);

  basevector kmer;

  int mistakesShownPos = 0, mistakesShownNeg = 0;

  while ( iKmer1 != table1.end() || iKmer2 != table2.end() ) {
    if ( iKmer1 != table1.end() 
	 && (iKmer2 == table2.end() || 
             kmerOnlyLess<K>()( *iKmer1, *iKmer2 ) ) ) {
      iKmer1->GetBasevector(kmer);
      missingByGc[GcBases(kmer)]++;
      if (LIST_KMERS) withMissingWriter->write( kmer_with_count<K>( kmer, 1 ) );

      if ( mistakesShownNeg++ < MAX_MISTAKES_TO_SHOW ) {
	cout << " false negative: " << *iKmer1 << endl;
      }
      
      if (evaluateReads && table3->IsStrong(kmer)) {
	inReadsByGc[GcBases(kmer)]++;
        if ( showMistakes ) {
          falseNegatives.push_back( *iKmer1 );
	}
      }
      ++iKmer1;
    } else if ( iKmer2 != table2.end()
		&& (iKmer1 == table1.end() || 
                    kmerOnlyLess<K>()( *iKmer2, *iKmer1 ) ) ) {
      iKmer2->GetBasevector(kmer);
      nonGenomicByGc[GcBases(kmer)]++;
      if (LIST_KMERS) nonGenomicWriter->write( kmer_with_count<K> (kmer, GcBases(kmer)) );
      if (LIST_KMERS) withMissingWriter->write( kmer_with_count<K>( kmer, 1 ) );
      if ( mistakesShownPos++ < MAX_MISTAKES_TO_SHOW ) {
	cout << " false positive: " << *iKmer2 << endl;
      }
      if ( showMistakes ) {
        falsePositives.push_back( *iKmer2 );
      }
      ++iKmer2;
    } else {
      iKmer1->GetBasevector( kmer );
      if (LIST_KMERS) withoutNonGenomicWriter->write( kmer_with_count<K>( kmer, 1 ) );
      if (LIST_KMERS) withMissingWriter->write( kmer_with_count<K>( kmer, 1 ) );
      ForceAssert( eq_kmer( *iKmer1, *iKmer2 ) ); 
      iKmer1->GetBasevector(kmer);
      sharedByGc[GcBases(kmer)]++;
      ++iKmer1;
      ++iKmer2;
    }
  }

  delete table3;

  if (LIST_KMERS) {
    nonGenomicWriter->close(); delete nonGenomicWriter;
    withoutNonGenomicWriter->close(); delete withoutNonGenomicWriter;
    withMissingWriter->close(); delete withMissingWriter;
  }

  int missing = 0;
  int nonGenomic = 0;
  int shared = 0;
  int inReads =0;
  for (int i = 0; i < K+1; ++i) {
    missing += missingByGc[i];
    nonGenomic += nonGenomicByGc[i];
    inReads += inReadsByGc[i];
    shared += sharedByGc[i];
    genomicByGc[i] = missingByGc[i] + sharedByGc[i];
    strongByGc[i] = nonGenomicByGc[i] + sharedByGc[i];
  }

  cout << "Strong Kmer Table Evaluation for K=" << K << "\n";
  cout << "-------------------------------------\n";
    
  
  cout << "GC\t[1]\t[2]\t[3]\t[4]\t[5]" << (evaluateReads ? "\t[6]" : "") << "\n";
  for (int i = 0; i < K+1; ++i)
    cout << i << "\t"
	 << genomicByGc[i] << "\t"
	 << strongByGc[i] << "\t"
	 << sharedByGc[i] << "\t"
	 << nonGenomicByGc[i] << "\t"
	 << missingByGc[i] << "\t"
	 << (evaluateReads ? ToString(inReadsByGc[i]) : "") << "\n";

  cout << "\n";
  cout << "[1] Genomic Kmers  : " << table1.isize() << "\n";
  cout << "[2] Strong Kmers   : " << table2.isize() << "\n";
  cout << "[3] Genomic Kmers in Strong     : " << shared << "\n";
  cout << "[4] Non Genomic Kmers in Strong : " << nonGenomic << "\n";
  cout << "[5] Genomic Kmers that are Missing from Strong          : " << missing << "\n";
  if (evaluateReads)
    cout << "[6] Genomic Kmers in Reads that are Missing from Strong : " << inReads << "\n";

  cout << "\n\n";

  if ( showMistakes ) {
    vecbasevector reads( readsFile );
    vecqualvector quals;
    if ( ! qualsFile.empty() )
      quals.ReadAll( qualsFile );
    vecbitvector trusted;
    if ( ! trustedFile.empty() )
      trusted.ReadAll( trustedFile );
    
    for ( int round = 0; round < 2; ++round ) {
      vec< kmer_with_count<K> >& mistakes = ( round == 0 ? falsePositives : falseNegatives );
      Sort( mistakes );

      // TODO: potentially dangerous truncation of index by these vec<int>'s
      vec< vec< int > > instanceIds( mistakes.size() );
      vec< vec< vec<int> > > instanceQuals( mistakes.size() );
      vec< vec< vec<int> > > instanceTrusteds( mistakes.size() );
      basevector kmer(K);
      cout << "Processing reads in " << reads.size() / 100000 << " passes: " << endl;
      for ( size_t i = 0; i < reads.size(); ++i ) {
        if ( i % 100000 == 0 ) Dot( cout, i / 100000 );
        for ( int j = 0; j < (int) reads[i].size() - (int)KSHAPE::getSpan() + 1; ++j ) {
          KSHAPE::extractKmer( kmer, reads[i], j );
          CanonicalForm cf = kmer.Canonicalize();
          kmer_with_count<K> target( kmer, 0 );
          // is this kmer an instance of the mistake?
          typename vec< kmer_with_count<K> >::iterator found = 
            lower_bound( mistakes.begin(), mistakes.end(), target, kmerOnlyLess<K>() );
          if ( ! eq_kmer( *found, target ) ) continue;
          int idx = distance( mistakes.begin(), found );
          // extract the associated quals
          basevector found_bases(K);
          vec<int> found_quals(K);
          vec<int> found_trusteds(K);
          // TODO: potentially dangerous truncation of index to play negative number trick
          int found_id = i;
          for ( int k = 0; k < K; ++k ) {
            int basePos = j + KSHAPE::getShapeOffset(k);
            found_bases.Set( k, reads[i][basePos] );
            if ( ! quals.empty() )
              found_quals[k] = quals[i][basePos];
            if ( ! trusted.empty() )
              found_trusteds[k] = trusted[i][basePos];
          }
          if ( cf == CanonicalForm::REV ) {
            found_quals.ReverseMe();
            found_trusteds.ReverseMe();
            found_bases.ReverseComplement();
            found_id = -found_id-1;
          }
          ForceAssert( found_bases == kmer );
          instanceIds[idx].push_back(found_id);
          if ( ! quals.empty() )
            instanceQuals[idx].push_back(found_quals);
          if ( ! trusted.empty() )
            instanceTrusteds[idx].push_back(found_trusteds);
        }
      }
      cout << endl;
      for ( unsigned int i = 0; i < mistakes.size(); ++i ) {
        cout << "False " << ( round == 0 ? "positive" : "negative" ) << " #" << i << ":" << endl;
        basevector kmer;
        mistakes[i].GetBasevector( kmer );
        cout << setw(10) << "read id";
        for ( int k = 0; k < K; ++k ) 
          cout << setw(3) << kmer.at(k);
        cout << endl;
        for ( unsigned int j = 0; j < instanceQuals[i].size(); ++j ) {
          int id = instanceIds[i][j];
          bool rc = (id<0);
          if (rc) id = -id-1;
          String idString = ToString(id) + (rc?"rc":"fw");
          cout << setw(10) << idString;
          for ( int k = 0; k < K; ++k )
            cout << setw(3) << instanceQuals[i][j][k];
          if ( instanceTrusteds[i].size() >= j &&
               Sum( instanceTrusteds[i][j] ) == K )
            cout << " trusted";
          cout << endl;
        }
        /*
         * for ( unsigned int j = 0; j < instanceTrusteds[i].size(); ++j ) {
         *   for ( int k = 0; k < K; ++k )
         *     cout << setw(3) << ( instanceTrusteds[i][j][k] ? "!" : " " );
         *   cout << endl;
         * }
         */
        cout << endl;
      }
    }
  }
}



int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA); 
  CommandArgument_String(RUN); 
  CommandArgument_String_OrDefault_Doc(READS_IN, "reads",
    "Original reads from which strong kmers were identified"); 
  CommandArgument_String_OrDefault_Doc(GENOME_IN, "genome",
    "Genome, against which the strong kmers are to be evaluated"); 
  CommandArgument_String_OrDefault_Doc(STRONG_IN, "",
    "Strong kmers to evaluate against truth (the genome)");  
  CommandArgument_Bool_OrDefault(SHOW_MISTAKES, False);
  CommandArgument_Int_OrDefault(MAX_MISTAKES_TO_SHOW, 1);
  CommandArgument_Bool_OrDefault_Doc(EVALUATE_READS, True,
    "Evaluate kmers in original reads in addition to strong and genomic kmers");
  CommandArgument_String_OrDefault_Doc(USE_AS_GENOMIC_KMERS, "",
    "Take the set of genomic kmers from this file, rather than from the genome");
  CommandArgument_Bool_OrDefault_Doc(REBUILD_GENOME_KMER_FREQS, False,
    "Rebuild genome frequency table, overwriting existing table");
  CommandArgument_Bool_OrDefault_Doc(REBUILD_READS_KMER_FREQS, False,
    "Rebuild reads frequency table, overwriting existing table");
  CommandArgument_Bool_OrDefault_Doc(BUILD_MISSING_TABLES, True,
    "Build missing frequency tables for genome and reads");
  CommandArgument_Bool_OrDefault_Doc(LIST_KMERS, False,
    "Write out lists of problem kmers (non-genomic, missing, etc)");
  CommandArgument_String(K);

  EndCommandArguments; 

  vec<KmerShapeId> Ks;
  ParseKmerShapeIdSet( K, Ks, false );
  ForceAssertSupportedKShapes( Ks );
  
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;


  String refFile = data_dir + "/" + GENOME_IN + ".fastb";
  String readsFile = run_dir + "/" + READS_IN + ".fastb";
  String qualsFile = run_dir + "/" + READS_IN + ".qualb";
  String trustedFile = data_dir + "/" + READS_IN + ".trusted";

  if ( ! IsRegularFile( trustedFile ) )
    trustedFile = "";
  if (STRONG_IN == "")
    STRONG_IN = READS_IN + ".strong";

  String genomeFileBase = data_dir + "/" + GENOME_IN + ".freq_table.k";
  String strongFileBase = run_dir + "/"+ STRONG_IN + ".k";
  String readsFileBase = run_dir +"/" + READS_IN + ".freq_table.k";

  for ( unsigned int i = 0; i < Ks.size(); ++i ) {
    #define CASE(K) \
        String Kid = ToString(K::getId()); \
	if (!BUILD_MISSING_TABLES && !IsRegularFile(genomeFileBase + Kid) ) {  \
          FatalErr("Could not find genome kmer frequency table for K=" + Kid); \
        } else if ( REBUILD_GENOME_KMER_FREQS || !IsRegularFile(genomeFileBase + Kid) ) { \
          vecbasevector allReads( refFile ); \
          cout << "Finding " << Kid << "-mer frequencies in genome." << endl; \
          WriteKmerFrequencies<K>( allReads, genomeFileBase + Kid, true ); } \
	if (EVALUATE_READS) { \
	  if (!BUILD_MISSING_TABLES && !IsRegularFile(readsFileBase + Kid) ) {  \
            FatalErr("Could not find reads kmer frequency table for K=" + Kid); \
          } else if (REBUILD_READS_KMER_FREQS || !IsRegularFile(readsFileBase + Kid) ) { \
            vecbasevector allReads( readsFile ); \
            cout << "Finding " << Kid << "-mer frequencies in uncorrected reads." << endl; \
            WriteKmerFrequencies<K>( allReads, readsFileBase + Kid, true ); } \
	} \
        Compare<K>( USE_AS_GENOMIC_KMERS.empty() ? genomeFileBase  + ToString(Ks[i]) : ( run_dir + "/" + USE_AS_GENOMIC_KMERS ), \
                    strongFileBase + ToString(Ks[i]), \
                    readsFileBase + ToString(Ks[i]), \
                    EVALUATE_READS, SHOW_MISTAKES, MAX_MISTAKES_TO_SHOW, LIST_KMERS, \
		    readsFile, qualsFile, trustedFile, run_dir )
    DISPATCH_ON_KSHAPE( Ks[i], CASE );
  }

  return 0;
}

