///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// TODO: we should be dealing with potentially selecting a sample instead of hard-coding everything for index 0

#include "MainTools.h"
#include "feudal/BinaryStream.h"
#include "FastaFileset.h"
#include "paths/simulation/VCF.h"
#include "Basevector.h"
#include "efasta/EfastaTools.h"

#include "graph/Digraph.h"
#include "paths/simulation/GenomeVCFToHyperBaseVector.h"
#include "paths/HyperKmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"

#include "paths/long/fosmid/FosmidPool.h"
#include "paths/HyperBasevector.h"


namespace {

void eFastaToHyper(const int K, const int NUM_THREADS, const VecEFasta& ereads, HyperKmerPath& h, HyperBasevector& hbv, const string& READS_OUT = "" )
{
    const int LRP_VERBOSITY = 0;

    cout << Date() << ": expanding eFasta reads" << endl;

    vecbasevector correctedv;
    vec< triple<int,int,int> > origin;
    for ( size_t id = 0; id < ereads.size( ); id++ ) {
	vec<basevector> b;
	ereads[id].ExpandTo(b);
	for ( size_t j = 0; j < b.size( ); j++ )
	    correctedv.push_back_reserve( b[j] );
    }

    if ( READS_OUT != "" )
	correctedv.WriteAll(READS_OUT);

    // Make paths.
    cout << Date() << ": pathing" << endl;
    vecKmerPath paths, paths_rc, unipaths;
    vec<tagged_rpint> pathsdb, unipathsdb;

    ReadsToPathsCoreY( correctedv, K, paths );
    CreateDatabase( paths, paths_rc, pathsdb );
    Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );


    cout << Date() << ": building adjacency graph" << endl;
    digraph A;
    BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths,
	      unipathsdb, A );
    BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );

    KmerBaseBroker kbb(K, paths, paths_rc, pathsdb, correctedv);

    // Create HyperBasevector.
    cout << Date() << ": making hyperbasevector" << endl;
    hbv.Initialize(h, kbb);
}

//---- dumpVcfChromosome - pretty printer has mostly been subsumed into the VCF class itself as VCFChromosome::Print()
void dumpVcfChromosome( const VCFChromosome& vcfc )
{
    for ( auto i = vcfc.begin(); i != vcfc.end(); ++i ) {
        cout << "NAME " << vcfc.name() << endl;
        cout << "POS  " << i->POS1 << endl;
        cout << "REF  " << i->REF << endl;
        cout << "ALT  (" << i->ALT.size() << " els)" << endl;
        cout << "FILTER " << i->FILTER << endl;
//        cout << "use_this[]=" << use_this[vcfc.idx(i)] << endl;
        cout << "INFO " << i->INFO << endl;
        cout << "GENOTYPE_FORMAT " << i->GENOTYPE_FORMAT << endl;
        for ( size_t j = 0; j < i->GENOTYPE.size(); j++ ) {
	     cout << "GENOTYPE " << i->GENOTYPE[j] << endl;
	     cout << "GT " << i->GT[j] << endl;
	     // next line broken for haploid.
	     cout << "GT_IDX " << i->GT_IDX[j][0] << "," << i->GT_IDX[j][1] << endl;
	     cout << "GT_PHASED " << i->GT_PHASED[j] << endl;
        }
    }
}

//---- eFastaChunkUpBetweenVariants - takes a single efasta seq and breaks, between variants {},
// into a VecEFasta with K-based overlaps
//
void eFastaChunkUpBetweenVariants( const efasta& seq, const int K, VecEFasta& efastav_out)
{
    efastav_out.clear();

    bool processing_bracket = false;

    for ( efasta::const_iterator i = seq.begin(); i != seq.end(); /* i */ ) {

	if (!processing_bracket) {
	    // new read -- okay, lets find us some brackets.
	    efastav_out.push_back( efasta() );
	    while ( i != seq.end() && *i != '{' )
		efastav_out.back().push_back( *i++ );

	    processing_bracket = true;
	} else {
	    // okay, we saw a bracket, let's find the end and then find a K-base space

	    // find the end bracket
	    while ( i != seq.end() && *i != '}' )
		efastav_out.back().push_back(*i++);

	    // find a K-sized space
	    int space = -1;	// we start off down one because we must consume the }  first
	    while ( i != seq.end() && *i != '{' && space < K ) {
		efastav_out.back().push_back(*i++);
		space++;
	    }

	    if ( space == K ) {		// Did we exit because we found a suitable space?
		processing_bracket=false;	// start a new read, next cycle
		i -= K;		// we want a K-sized overlap with the next read, so back up
	    }
	}
    }
}


//---- pushUnique - takes a list of alternates (e.g. {AT,A}), removes duplicates, and adds to efasta
//
void pushUnique( const vector<string>& alts, efasta& seq )
{
    map<string,bool> flipflop;
    for ( const auto& alt : alts )
	flipflop[alt] = true; 			// true entry for each unique

    size_t count = flipflop.size();
    if ( count == 1 ) {				// everything is the same
	seq.append( alts[0] );			// ... so only append one copy
    } else {
	seq.push_back('{');
	for ( const auto& alt : alts ) {
	    if ( flipflop[alt] ) {		// not encountered yet
		flipflop[alt] = false;		// ... so mark it seen
		seq.append( alt );		// ... append
		if ( --count > 0 )		// ... add a comma if more to come
		    seq.push_back(',');
	    }
	}
	seq.push_back('}');
    }
}


//---- phasedVcfToEfasta - this converts genome+vcf to efasta for a specified range of a single chromosome.  Output in "seq."
//
// Here are the rules:
// - as with everything in this program, we only look at the genotype associated with the first sample
// - if a genotype is unspecified, then we output the reference base at that coordinate
// - we *always* output alternatives, even if they are equal, e.g. {C,C} -- this will be cleaned up later
// - phasing spans from the first phased variant (|) until a non-phased variant occurs (/).  In order to represent
//   this in efasta, we therefore all phased calls and interventing bases into alternatives -- e.g. {CGATA,TGATT}
//
void phasedVcfToEFasta( const vecbasevector& genome, const VCF& vcf, vector<bool>& use_this, const Range range, efasta& seq, const bool debug = false, const int sample = 0 )
{
    ForceAssert( vcf.size() == 1 );             // sanity check -- should have loaded only one chromosome
    const VCFChromosome& vcfc = *(vcf.begin());

    cout << Date() << ": using PHASED algorithm" << endl;

    //---- translate to a single efasta record
    const basevector& ref_seq = genome[ range.chr0() ];

    auto edit_iter = vcfc.begin();

    for (size_t ref_pos0 = range.start0(); ref_pos0 < range.end0(); ) {	// for each base

	// echo bases until an edit
	while ( ref_pos0 < range.end0() && ( edit_iter == vcfc.end() || ref_pos0 < edit_iter->POS1 - 1) )
	    seq.push_back( Base::val2Char( ref_seq[ref_pos0++] ) );


	// if we have an edit...
	if ( ref_pos0 < range.end0() ) {
	    ForceAssertLt( edit_iter->POS1-1, range.end0());

	    if (debug)
		cout << "checking edit at POS1=" << edit_iter->POS1 << ", ref_pos0=" << ref_pos0 << endl;

	    // if it's an unused edit outside of a phase-set, then skip it
	    if ( ! use_this[ vcfc.idx( edit_iter ) ] ) {
		if (debug) {
		    cout << "FILTERING (skipping) idx=" << vcfc.idx(edit_iter) << endl;
		    edit_iter->Print(cout);
		}
		edit_iter++;
		continue;
	    }

	    // double check that there are no phase sets (PS tags) as we do not support them, yet
	    if ( edit_iter->GT_PHASED[sample] && edit_iter->GENOTYPE_FORMAT.find("PS") != string::npos ) {
		cout << "WARNING: we do not yet support phase-sets, however this record has one:" << endl;
		cout << "\tCHROM " << vcfc.name() << ", POS1=" << edit_iter->POS1 << endl;
		cout << "\tFORMAT***** " << edit_iter->GENOTYPE_FORMAT << endl;
		cout << "\tGENOTYPE[" << sample << "] " << edit_iter->GENOTYPE[sample] << endl;
	    }

	    // print out the edit for debugging
	    if ( debug ) edit_iter->Print(cout);


	    // since we may have phase information across edits, we're going to be less smart about "compressing"
	    // identical alternatives (e.g. homozygous positions).  So we'll always open and close braces.  Anyway,
	    // the pathing later will clean it up.
	    vector<string> alts;
	    alts.push_back(string());

	    size_t nGt = edit_iter->GT_IDX[sample].size();		// number of alleles


	    // for each allele
	    for ( size_t igt = 0; igt < nGt; ++igt ) {
		size_t reset = ref_pos0;				// we're going to come back here repeatedly for each allele
		auto reset_iter = edit_iter;

		if (debug) cout << "ALLELE " << igt << endl;

		// for the entire phase-set
		while ( ref_pos0 < range.end0() ) {			// we keep going until we hit an unphased edit OR the end

		    // we may be walking in-between phased changes
		    if ( ref_pos0 < edit_iter->POS1 - 1 ) {
			while ( ref_pos0 < edit_iter->POS1-1 ) {
			    alts.back().push_back( Base::val2Char( ref_seq[ ref_pos0 ] ) );			// ... so output the reference base
			    ref_pos0++;
			}
		    }

		    // now we're at an edit... should we use it?
		    if ( !use_this[vcfc.idx(edit_iter)] ) {
			// no, echo base at that position.
			alts.back().push_back( Base::val2Char( ref_seq[ ref_pos0 ] ) );			// ... so output the reference base
			ref_pos0++;
		    } else {
			// yes, expand the edit
			char alt = edit_iter->GT_IDX[sample][igt]-1;
			if ( alt < -1 ) cout << "WARNING: missing alternative at reference 0-position " << ref_pos0 << endl;

			if ( alt <= -1 ) {					// if the alternative is missing, we're going to take the reference.
			    alts.back().append( edit_iter->REF.ToString() );
			    if (debug) cout << "output " << edit_iter->REF.ToString() << endl;
			} else if ( ! edit_iter->ALT[alt].isBaseVec() ) {		// otherwise, is it a literal sequence alternative?
			    cout << "WARNING: non-basevec alternative " << edit_iter->ALT[alt].sRep << " at 0-position " << ref_pos0 << endl;
			    alts.back().append( edit_iter->REF.ToString() );			// ... if not, warn but use the reference
			    if (debug) cout << "output " << edit_iter->REF.ToString() << endl;
			} else {
			    alts.back().append( edit_iter->ALT[ alt ].bRep.ToString() );	// literal sequence -- output the alternate
			    if (debug) cout << "output " << edit_iter->ALT[ alt ].bRep.ToString() << endl;
			}


			ref_pos0 += edit_iter->REF.size();

		    }

		    if ( !edit_iter->GT_PHASED[sample] ) {				// no phasing, so leave this "phase-set" loop
			++edit_iter;
			break;
		    }


		    if ( ++edit_iter == vcfc.end() || !edit_iter->GT_PHASED[sample] ) {	// next sample will break phasing, so leave.
			break;
		    }


		}	// while ref_pos0


		// if there are more alleles, then push a comma and reset our counters
		if ( igt+1 < nGt ) {
		    alts.push_back(string());
		    ref_pos0 = reset;
                    edit_iter = reset_iter;
		}

	    }		// finished all alleles for the singleton or phase-set

	    // de-dup alternatives and add to seq
            if ( debug ) cout << "about to push alts, ref_pos0=" << ref_pos0 << endl;
	    pushUnique(alts, seq);
	}	// if ref_pos0 < range.end0()
    }	// for ref_pos0
}





//---- validateVcfAgainstGenome -- check that REF bases in the VCF match the genome in the RANGE selected
//
void validateVcfAgainstGenome( const vecbasevector& genome, const VCF& vcf, const Range& range )
{


    ForceAssert( vcf.size() == 1 );              // sanity check -- should have loaded only one chromosome
    const VCFChromosome& vcfc = *(vcf.begin());

    cout << Date() << ": validating refence against " << vcfc.size() << " edits." << endl;

    const basevector& ref_seq = genome[ range.chr0() ];

    for ( auto iter = vcfc.begin(); iter != vcfc.end(); ++iter ) {

	for ( size_t off = 0; off < iter->REF.size(); ++off ) {
	    size_t ref_pos0 = iter->POS1 - 1 + off;

	    ForceAssertLt(ref_pos0, ref_seq.size() );

	    if ( iter->REF[off] != ref_seq[ ref_pos0 ] ) {
		cout << "Error reference does not match VCF file at ref seq-index=" << range.chr0() << ", 0-pos=" << ref_pos0 << endl;
		cout << "reference base is " << Base::val2Char(ref_seq[ ref_pos0 ]) << endl;
		cout << "VCF does not match at 0-pos=" << off << " of REF string=" << iter->REF.ToString() << endl;
		iter->Print(cout);
		FatalErr("exiting...");
	    }
	}
    }

    cout << Date() << ": VCF matches reference within the selected range." << endl;
}

//---- unphasedVcfToEFasta - this is the original, unphased algorithm to be replaced with one that respects phasing
//
void unphasedVcfToEFasta( const vecbasevector& genome, const VCF& vcf, vector<bool>& use_this, const Range range, efasta& seq, const bool debug = false, const int sample = 0 )
{
    ForceAssert( vcf.size() == 1 );              // sanity check -- should have loaded only one chromosome
    const VCFChromosome& vcfc = *(vcf.begin());

    cout << Date() << ": using UNPHASED algorithm" << endl;

    //---- translate to a single efasta record
    const basevector& ref_seq = genome[ range.chr0() ];
    size_t ref_pos0 = range.start0();

    if ( find( use_this.begin(), use_this.end(), true ) != use_this.end() ) {  // normal case -- there are still edits left over

	for ( auto i = vcfc.begin(); i != vcfc.end(); ++i ) {

	    ForceAssertLt( ref_pos0, range.end0() );	// shouldn't if we've picked from VCF properly

	    // skip if not used
	    if ( !use_this[vcfc.idx(i)] ) continue;

	    // emit bases up until the edit
	    size_t pos0 = i->POS1-1;
	    for ( /* */ ; ref_pos0 < pos0; ref_pos0++ )
		seq.push_back( Base::val2Char( ref_seq[ref_pos0] ) );

	    // double check that there are no phase sets as we do not support them, yet
	    if ( i->GT_PHASED[sample] && i->GENOTYPE_FORMAT.find("PS") != string::npos ) {
		cout << "WARNING: we do not yet support phase-sets, however this record has one:" << endl;
		cout << "\tCHROM " << vcfc.name() << ", POS1=" << i->POS1 << endl;
		cout << "\tFORMAT***** " << i->GENOTYPE_FORMAT << endl;
		cout << "\tGENOTYPE[" << sample << "] " << i->GENOTYPE[sample] << endl;
	    }


	    // expand out the edit
	    if (debug) i->Print(cout);

	    size_t nGt = i->GT_IDX[sample].size();
	    ForceAssert(nGt > 0);

	   if (debug) cout << "***** PROCESSING...." << endl;

	    // now we go through the alternatives and make a list of the valid,
	    // unique ones to instantiate later.
	    vector<char> alts;
	    for ( size_t j = 0; j < nGt; j++ ) {
		if (debug) cout << "ALLELE" << j << endl;
		char alt = i->GT_IDX[sample][j];

		if (debug) cout << "alt=" << alt << endl;

		if ( alt < 0 ) {
		    alt = 0;
		    if ( debug ) cout << "missing alternative, forcing reference" << endl;
		}

		if ( static_cast<size_t>(alt) > i->ALT.size() )
		    i->Print(cout);
		ForceAssertLe(static_cast<size_t> (alt), i->ALT.size() );

		// now we are free to look at the alt
		if ( alt == 0 || i->ALT[alt-1].isBaseVec() ) {
		    // push back uniques only
		    if ( std::find( alts.begin(), alts.end(), alt ) == alts.end() )
			alts.push_back(alt);
		} else if ( i->ALT[alt-1].sRep != "." ) {	 // just ignore '.', but bark about others
		    cout << Date() << ": WARNING - found <ID> format alternative in file." << endl;
		}
	    }


	    // okay, so now we add the sequence for the alternatives to the
	    // output efasta in seq.
	    if ( alts.size() > 1 ) seq.push_back('{');
	    for ( size_t j = 0; j < alts.size(); ++j ) {
		char alt = alts[j];

		if ( alt == 0 )
		    seq.append(i->REF.ToString() );		// use the reference sequence
		else
		    seq.append(i->ALT[alt-1].bRep.ToString() ); // supplied alternate sequence

		if ( j + 1 < alts.size() ) seq.push_back(',');
	    }
	    if ( alts.size() > 1 ) seq.push_back('}');


	    // and advance past the reference bases
	    ref_pos0 += i->REF.size();
	}

	// finish out past the last edit
	for ( /* ref_pos0 */ ; ref_pos0 < range.end0(); ref_pos0++ )
		seq.push_back( Base::val2Char( ref_seq[ref_pos0] ) );

    }
    else {
	// translate the whole range to an efasta record
	for ( size_t i = range.start0() ; i < range.end0(); i++ )
		seq.push_back( Base::val2Char( ref_seq[i] ) );
    }
}


void genomeVcfToEFasta( const vecbasevector& genome, const VCF& vcf, const vec<String>& filters, const int K, const Range range, VecEFasta& efastav_out, const bool DEBUG = false, const bool PHASING = false, const int sample_idx = 0 )
{

    //---- check whether there are any edits
    if ( vcf.size() > 0 ) {

	ForceAssert( vcf.size() == 1 );              // sanity check -- should have loaded only one chromosome (or none, if there were no matches)

	const VCFChromosome& vcfc = *(vcf.begin());

	//    dumpVcfChromosome(vcfc);

	vector<bool> use_this( vcfc.size(), true );

	//---- filter out non-passing elements
	size_t negatives0 = 0;		// keep track of how many we filter out

	// skip filters if filters array is empty (size 0) or if it contains a single
	// empty string (size 1, filters[0]=="").
	if ( filters.size() > 1 || ( filters.size() == 1 && filters[0] != "" ) ) {
	    for ( auto i = vcfc.begin(); i != vcfc.end(); ++i ) {
		if ( !Member( filters, String(i->FILTER) ) ) {
		    use_this[vcfc.idx(i)] = false;
		    negatives0++;
		}
	    }
	}

	cout << "excluded " << negatives0 << " of " << use_this.size() << " based on (not) matching VCF filters:" << endl;
	for ( size_t i = 0; i < filters.size(); ++i )
	    cout << filters[i] << ( i+1 < filters.size() ? ", " : "" );
	cout << endl;

	//---- find adjacent overlapping changes and favor left-most until no overlapping changes are found
	size_t negatives1 = 0;
	bool found;
	do {
	    size_t last_right = 0;	// end of HALF-OPEN interval, so first pos AFTER
	    found = false;

	    for ( auto i = vcfc.begin(); i != vcfc.end(); ++i ) {
		if (use_this[vcfc.idx(i)]) {
		    if ( i->POS1 < last_right ) { 	// items are sorted, so this means overlap
			use_this[vcfc.idx(i)] = false;
			found = true;
			negatives1++;
		    }
		    else
			last_right = i->POS1 + i->REF.size();
		}
	    }

	} while ( found );

	size_t negatives = negatives0 + negatives1;

	cout << "excluding " << negatives1 << " of " << use_this.size() << " based on overlap filtering." << endl;
	cout << "after filtering, " << use_this.size()  - negatives << " of " << use_this.size() << " remain." << endl;

	// TODO: how do we deal with N's in the FASTA since we're reading FASTB + other composite bases?

	//---- translate to a single efasta record
	efasta seq;
	if ( PHASING )
	    phasedVcfToEFasta(genome, vcf, use_this, range, seq, DEBUG, sample_idx);
	else
	    unphasedVcfToEFasta(genome, vcf, use_this, range, seq, DEBUG, sample_idx);

	//---- okay, now CHUNK the efasta into reads with K overlap
	eFastaChunkUpBetweenVariants( seq, K, efastav_out );

    } else {

	//---- No edits, so just extract the reference
	efastav_out.assign(1, efasta() );
	efasta& seq = efastav_out.back();
	const basevector& ref_seq = genome[ range.chr0() ];

	for ( size_t ref_pos0 = range.start0(); ref_pos0 < range.end0(); ++ref_pos0 )
	    seq.push_back( Base::val2Char( ref_seq[ref_pos0] ) );

    }
    cout << Date() << ": done making eFasta " << endl;

}


void dumpDot( const HyperBasevector& hbv, const String& name  )
{
         const Bool DOT_LABEL_CONTIGS = True;
         const Bool DOT_LABEL_VERTICES = False;
         vec<double> lengths( hbv.EdgeObjectCount( ) );
         for ( int i = 0; i < hbv.EdgeObjectCount( ); i++ )
              lengths[i] = hbv.EdgeLengthKmers(i);
         vec<String> edge_id_names( hbv.EdgeObjectCount( ) );
         for ( int i = 0; i < hbv.EdgeObjectCount( ); i++ )
             edge_id_names[i] = ToString(i);
         Ofstream( dout, name );
         hbv.PrettyDOT( dout, lengths, HyperBasevector::edge_label_info(
              HyperBasevector::edge_label_info::DIRECT, &edge_id_names ),
              DOT_LABEL_CONTIGS, DOT_LABEL_VERTICES, NULL, NULL, NULL, NULL,
              NULL, NULL );

}

};	//---- namespace local


int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;

     CommandArgument_String_Doc(VCF_FILE, "VCF file with variants called on GENOME");
     CommandArgument_String_Doc(RANGE, "e.g. chrX:0-10000 chr name into FASTA genome and zero-based range of bases to select.  "\
					 "Alternatively, you can specify a fosmid pool region number.");
     CommandArgument_String_Doc(HBV, "outputs HyperBasevector file HBV");
     CommandArgument_Int(K);
     CommandArgument_Int_OrDefault(NUM_THREADS, 0);
     CommandArgument_Bool_OrDefault_Doc(DEBUG, false, "turn on extra output");
     CommandArgument_Bool_OrDefault_Doc(PHASING, false, "used phased algorithm");
     CommandArgument_String_OrDefault_Doc(VCF_CHR_PREFIX, "", "HORRIBLE kludge for vcf file with a different chromosome name than GENOME_*");
     CommandArgument_StringSet_OrDefault_Doc(FILTERS, "PASS", "set of acceptable filter states (default: PASS). "\
					     "Specify FILTERS='' if you'd like to not filter.");
     CommandArgument_String_OrDefault_Doc(EFASTA_OUT, "", "dump the efasta file");
     CommandArgument_String_OrDefault_Doc(READS_OUT, "", "dump the expanded efasta as 'reads' to this FASTB");
     CommandArgument_String_OrDefault_Doc(FASTA_OUT, "", "dump the HyperBasevector to FASTA");
     CommandArgument_String_OrDefault_Doc(FASTB_OUT, "", "dump the HyperBasevector to FASTB");
     CommandArgument_String_OrDefault_Doc(VCF_OUT, "", "dump the matching lines to a VCF (plus the headers)");
     CommandArgument_String_OrDefault_Doc(GENOME_FASTA, "", "GENOME to read in -- FASTA format.");
     CommandArgument_String_OrDefault_Doc(GENOME_FASTB, "/wga/scr4/bigrefs/human19/genome.fastb", "GENOME to read in -- FASTB format.");
     CommandArgument_String_OrDefault_Doc(GENOME_NAMES, "/wga/scr4/bigrefs/human19/genome.fastb.names", "names file for GENOME read in.");
     CommandArgument_String_OrDefault_Doc(DOT,"","outputs DOT");
     CommandArgument_String_OrDefault_Doc(SAMPLE_NAME,"","if sample name is not selected, then the *first* sample is chosen.");
     CommandArgument_Int_OrDefault_Doc(EXTEND_BY,0,"extend ends of region by this amount");

     EndCommandArguments;

     NUM_THREADS=configNumThreads(NUM_THREADS);


     if ( PHASING ) cout << "**** WARNING: phased code has not been throughly vetted yet" << endl;

     vecbasevector genome;
     vecString chr_names;

     if ( GENOME_FASTB != "" ) {
	 cout << Date() << ": reading genome from " << GENOME_FASTB << endl;
	 genome.ReadAll(GENOME_FASTB);

	 if ( GENOME_NAMES != "" ) {
	     cout << Date() << ": reading genome names from " << GENOME_NAMES << endl;
	     chr_names.ReadAll(GENOME_NAMES);
	     // remove chromosome names from the first space onward.

	     for_each(chr_names.begin(), chr_names.end(),
		     [] (String& name) {
			 size_t pos = name.find_first_of(' ');
			 if ( pos != string::npos ) name.resize(pos);
		     }
	     );
	     cout << Date() << ": read " << chr_names.size() << " names" << endl;
	 }
     }
     else if ( GENOME_FASTA != "" ) {
	 FastFetchReads(genome, &chr_names, GENOME_FASTA );
     } else
	 FatalErr("you must specify one of GENOME_FASTB or GENOME_FASTA.");

     Range coords(RANGE, chr_names);
     coords.extend_by(EXTEND_BY);

//     cout << Date() << ": processing range " << coords.chr0() << ":" << coords.start0() << "-" << coords.end0() << endl;

     cout << Date() << ": processing chromosome name " << coords.chr() << ", range1 [" << coords.start1() << "," << coords.end1() << ")" << endl;

     VCF variants(VCF_FILE, VCF_CHR_PREFIX+coords.chr(), coords.start1(), coords.end1(), VCF_OUT);

     // Validate VCF against the GENOME
     validateVcfAgainstGenome(genome, variants, coords);

     // find the sample that we're looking for
     int sample_idx = 0;
     if ( SAMPLE_NAME != "" ) {
	 sample_idx = variants.getSampleIndex(SAMPLE_NAME);
	 if (sample_idx < 0) FatalErr("SAMPLE_NAME specified (" << SAMPLE_NAME << ") was not found in the VCF file");
     }

     // make Efasta
     VecEFasta efastav;
     genomeVcfToEFasta( genome, variants, FILTERS, K, coords, efastav, DEBUG, PHASING, sample_idx );

     if ( EFASTA_OUT != "" ) {
	 cout << Date() << ": writing efasta to " << EFASTA_OUT << endl;
	 ofstream out( EFASTA_OUT, ios::out );
	 for ( size_t i = 0; i < efastav.size(); ++i ) {
	     ostringstream s;
	     s << i;
	     efastav[i].Print(out, s.str().c_str());
	 }
	 out.close();
     }

     // make HBV -- maybe this process gets repeated for each contig?
     cout << Date() << ": making Hyper" << endl;
     HyperBasevector hbv;
     HyperKmerPath h;
     eFastaToHyper( K, NUM_THREADS, efastav, h, hbv, READS_OUT );
     if ( !hbv.Acyclic( ) ) cout << "Warning: HyperBasevector has a cycle." << endl;

     // write outputs
     cout << Date() << ": writing HyperBasevector to " << HBV << endl;
     BinaryWriter writer( HBV );
     hbv.writeBinary(writer);

     if ( DOT != "" ) {
	 cout << Date() << ": writing DOT to " << DOT << endl;
	 dumpDot( hbv, DOT );
     }

     if ( FASTB_OUT != "" ) {
	 cout << Date() << ": Writing FASTB to " << FASTB_OUT << endl;
	 hbv.DumpFastb(FASTB_OUT);
     }

     if ( FASTA_OUT != "" ) {
	 cout << Date() << ": Writing FASTA to " << FASTA_OUT << endl;
	 hbv.DumpFasta(FASTA_OUT, False);
     }

     cout << Date( ) << ": Done!" << endl;


     return 0;
}
