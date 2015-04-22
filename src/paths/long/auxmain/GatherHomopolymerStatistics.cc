///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Jun 26, 2013
//

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "Basevector.h"
#include "random/Random.h"
#include "simulation/BamReadGenerator.h"
#include "util/TextTable.h"
#include <array>


// strategy:
//
// 1. scan the genome for homopolymer runs, maybe going out to a max distance from random restarts
// 2. widen the range a bit for context and capture and find reads overlapping the entire range
//    (necessary to make an inference about total length, obviously)
// 3. calculate homopolymer length (maybe do alignment?  errors in flanking regions?)


class FindHomopolymer {
public:
    struct HRun {
	unsigned char base;
	size_t chr;
	size_t pos;
	size_t len;
	BaseVec run;
	BaseVec lcontext;
	BaseVec rcontext;

	friend std::ostream&  operator<<( std::ostream& out, HRun const& self ) {
	    out << self.len << " base " << Base::val2Char(self.base) << "-run at " << self.chr << ":" << self.pos << "-" << self.pos+self.len << "   ";
	    out << self.lcontext.ToString() << " | " << self.run.ToString() << " | " << self.rcontext.ToString();
	    return out;
	}

	std::string samRegionString() {
	    std::ostringstream s;
	    s << chrTranslate() << ":" << pos+1 << "-" << pos+1+len-1;
	    return s.str();
	}
    private:
	String chrTranslate() {
	    if ( chr < 22 )
		return ToString(chr+1);
	    else if ( chr == 22 )
		return ToString("X");
	    else if ( chr == 23 )
		return ToString("Y");
	    FatalErr("bug -- invalid chromosome number");
	    return ""; //NOTREACHED
	}
    };

    FindHomopolymer(BaseVecVec const& gvv, size_t min_len, size_t max_len, size_t context) :
	_gvv(gvv), _min_len(min_len), _max_len(max_len), _context(context) {};

    FindHomopolymer() = delete;

    void seed( int seed ) { srandomx( seed ); }

    HRun findHRun();

private:
    BaseVecVec const& _gvv;
    size_t _min_len;
    size_t _max_len;
    size_t _context;
    size_t const _scan = 10000;		// how many bases to scan from random starting point
    size_t const _max_tries = 100;	// how many random restarts before we give up
};



FindHomopolymer::HRun FindHomopolymer::findHRun()
{
    ForceAssertGe(_context, 1U);		// we always look 1 ahead, so context must be > 1

    for ( size_t tries = 0; tries < _max_tries; ++tries ) {

	// pick a random starting coordinate
	size_t chr = randint( _gvv.size() );
	BaseVec const& gv = _gvv[chr];
	size_t pos0 = big_random() % ( gv.size() - (2*_context + _max_len) ) + _context;

	// scan for homopolymer
	size_t limit = std::min( gv.size() - _context, pos0+_scan );

	if ( gv[pos0-1] == BASE_A || gv[pos0-1] == BASE_T ) continue;	// don't start in the middle of a homopolymer

	for ( size_t pos = pos0; pos < limit ; ++pos ) {
	    unsigned char b = gv[pos];

	    if ( b == BASE_A || b == BASE_T ) {

		size_t pos1 = pos+1;
		while ( pos1 < limit && gv[pos1] == b ) ++pos1;

		size_t len = pos1 - pos;
		if ( pos1 < limit && len <= _max_len && len >= _min_len ) {	// must have been found
		    BaseVec run(gv, pos, len);
		    BaseVec lcontext( gv, pos - _context, _context );
		    BaseVec rcontext( gv, pos + len, _context );
		    return HRun{ b, chr, pos, len, run, lcontext, rcontext };
		} else
		    pos = pos1-1;
	    }
	}
    }

    std::ostringstream s;
    s <<  "unlucky today -- didn't find a homopolymer after " << _max_tries << " tries";
    FatalErr(s.str());

    return HRun();	//NOTREACHED
}


class CentralRates {
public:
    CentralRates(size_t min, size_t max) : _min(min), _max(max) {
	for ( size_t i = _min; i <= _max; ++i )
	    _counts[i] = {{0,0,0}};
    };

    bool isSufficient() {
	for ( auto iter = _counts.cbegin(); iter != _counts.cend(); ++iter )
	    if ( iter->second[0] < _min_samples || iter->second[2] < _min_samples )
		return false;
	return true;
    }

    bool isNeeded(size_t sz) {
	return ( _counts[sz][0] < _min_samples || _counts[sz][2] < _min_samples );
    }

    bool addLine( size_t hlen, std::array<size_t, 5> const& sample ) {
	ForceAssertLe(hlen, _max);
	ForceAssertGe(hlen, _min);
	if ( isCentral( sample ) ) {
	    _counts[hlen][0] += sample[1];
	    _counts[hlen][1] += sample[2];
	    _counts[hlen][2] += sample[3];
	    return true;
	} else return false;
    }

    bool isCentral( std::array<size_t,5> const& sample ) {
	return ( sample[2] > _center_factor*(sample[0]+sample[1]+sample[3]+sample[4]) );
    }

    void Print(std::ostream& out ) {
	TextTable t;
	for ( auto iter = _counts.begin(); iter != _counts.end(); ++iter )
	    t << iter->first << Tab << iter->second[0] << Tab << iter->second[1] << Tab << iter->second[2] << EndRow;
	t.Print(out,2);
    }

    std::array<size_t,3> const& operator[](const size_t i) { return _counts[i]; };

    size_t Min() {return _min;}
    size_t Max() {return _max;}

private:
    size_t _min, _max;
    std::map<size_t, std::array<size_t,3> > _counts;
    const size_t _min_samples = 30;
    const double _center_factor = 5.;
};


int main(int argc, char* argv[])
{
    std::string const DEFAULT_GENOME = "/wga/scr4/bigrefs/human19/genome.fastb";
    size_t const UNIQUE_CONTEXT = 5;

    RunTime();
    BeginCommandArguments;
    CommandArgument_String_OrDefault_Doc(GENOME, DEFAULT_GENOME, "specify genome to scan (and to which reads are aligned)");
    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, "number of threads to use in parallel sections");
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_LEN, 10, "minimum length of reference homopolymer to consider");
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_LEN, 25, "maximum length of reference homopolymer to consider");
    CommandArgument_UnsignedInt_OrDefault_Doc(CONTEXT, 10, "context bases around homopolymer to print");
    CommandArgument_UnsignedInt_Doc(NUM_INSTANCES, "number of such instances to find");
    CommandArgument_StringSet_Doc(BAMS,"set list of BAM files to find reads in");
    CommandArgument_Bool_OrDefault_Doc(POOL, true, "whether to pool results across BAM files");
    CommandArgument_Bool_Doc(HUMAN, "must set HUMAN to true accepting that you understand that this is for human only");
    CommandArgument_Bool_OrDefault(VERBOSE, false);
    CommandArgument_String_OrDefault_Doc(TABLE,"", "file to which the final table can be written as space-delimited");
    CommandArgument_Bool_OrDefault_Doc(CCR, false, "compute central rates - whether to keep going until sufficient samples are generated")
    CommandArgument_Int_OrDefault_Doc(SEED,0,"if non-zero, seed for random number generator -- ensures a given output sequence");
    EndCommandArguments;

    if ( !HUMAN ) FatalErr("only implemented for human data for now -- issue is contig names and samtools");

    ForceAssertGe(MIN_LEN, 3U);			// we're going to remove 2 characters, so there better be at least 3.

    NUM_THREADS = configNumThreads(NUM_THREADS);
    omp_set_num_threads(NUM_THREADS);

    // read in the genome
    cout << Date() << ": reading in genome" << endl;
    BaseVecVec gvv( GENOME );

    cout << Date() << ": kludge for now -- removing contigs outside of [0,23] chr1-22,X,Y" << endl;
    gvv.resize(24);

    cout << Date() << ": genome read, " << gvv.size() << " contigs" << endl;

    cout << Date() << ": Searching for " << ( CCR ? "(at least) " : "" ) << NUM_INSTANCES << " homopolymers of length ["
	    << MIN_LEN << "," << MAX_LEN << "]" << endl;

    FindHomopolymer fh( gvv, MIN_LEN, MAX_LEN, CONTEXT );
    if ( SEED ) fh.seed(SEED);

    // gather homopolymer runs
    size_t nfound = 0;

    TextTable output;
    output << "left flank" << Tab << "right flank" << Tab << "len" << Tab << "I-2" << Tab;
    output << "I-1" << Tab << "I" << Tab << "I+1" << Tab << "I+2" << Tab << "coords";
    if ( !POOL ) output << Tab << "BAM";
    output << EndRow << DoubleLine;

    if (VERBOSE) cout << endl;
    CentralRates crates( MIN_LEN, MAX_LEN );
    while ( 1 ) {

	FindHomopolymer::HRun run = fh.findHRun();
	if ( CCR && !crates.isNeeded(run.len) ) continue;
	if ( VERBOSE ) {
	    cout << "found " << nfound << endl;
	    cout << run << endl;
	}
	Dot(cout, nfound);

	ForceAssertGe(run.lcontext.size(), UNIQUE_CONTEXT);
	ForceAssertGe(run.rcontext.size(), UNIQUE_CONTEXT);
	String lcontext = run.lcontext.ToString().substr(run.lcontext.size()-UNIQUE_CONTEXT, UNIQUE_CONTEXT);
	String rcontext = run.rcontext.ToString().substr(0, UNIQUE_CONTEXT);
	const char base = Base::val2Char(run.base);
	std::array<String,5> matches{ { lcontext+String(run.len-2, base)+rcontext,
				     lcontext+String(run.len-1, base)+rcontext,
				     lcontext+String(run.len,base) + rcontext,
				     lcontext+String(run.len+1,base)+rcontext,
				     lcontext+String(run.len+2, base)+rcontext } };
//	std::vector<size_t> counts, zero = {0, 0, 0, 0, 0};

	// now visit each BAM file, load the reads associated with this locus, and grep for the string of interest
	// not terribly efficient, but probably parsimonious
	const std::array<size_t,5> zero{{0,0,0,0,0}};
	std::array<size_t,5> aggregate = zero;

	ForceAssertEq(matches.size(), zero.size());

	for ( auto const& bam : BAMS ) {
	    std::array<size_t,5> counts = zero;

	    if ( VERBOSE ) cout << "visiting bam " << bam << endl;
	    if ( !POOL || &bam == &BAMS.front() ) counts = zero;
	    BamReadGenerator reader( CachedBAMFile(bam), run.samRegionString() );
	    BaseVecVec const& reads = reader.getReads();

	    for ( auto const& read : reads ) {
		for ( size_t i = 0; i < matches.size(); ++i ) {
		    if ( read.ToString().Contains( matches[i] ) || ReverseComplement(read).ToString().Contains( matches[i] ) ) {
			counts[i]++;
		    }
		}
	    }

	    bool central = crates.addLine(run.len, counts);

	    if ( POOL ) {
		if ( central )
		    for ( size_t i = 0; i < aggregate.size(); ++i )
			aggregate[i] += counts[i];
		if ( &bam == &BAMS.back() ) {
		    output << run.lcontext.ToString() << Tab << run.rcontext.ToString() << Tab << run.len << Tab;
		    for ( size_t i = 0; i < aggregate.size(); ++i )
			output << aggregate[i] << Tab;
		    output << run.samRegionString();
		    output << EndRow;
		}
	    } else {
		aggregate = counts;
		output << run.lcontext.ToString() << Tab << run.rcontext.ToString() << Tab << run.len << Tab;
		if ( !central ) output << "(";
		for ( size_t i = 0; i < aggregate.size(); ++i )
		    output << aggregate[i] << Tab;
		if ( !central ) output << ")";
		output << run.samRegionString();
		output << Tab << bam.SafeAfterLast("/");
		output << EndRow;
	    }
	    if ( !POOL && &bam == &BAMS.back() ) output << SingleLine;
	}


	if ( VERBOSE ) cout << endl;
	nfound++;

	if ( nfound % 20 == 0 ) {
	    std::cout << std::endl << "nfound=" << nfound << std::endl;
	    crates.Print(std::cout);
	}


	// termination criteria
	if ( !CCR && nfound >= NUM_INSTANCES ) break;
	else if ( CCR && nfound >= NUM_INSTANCES && crates.isSufficient() ) break;
    }

    std::cout << endl;
    output.Sort(2,true,2);
    output.Print(cout,2,"llrrrrrr");

    if ( TABLE != "" ) {
	ofstream outfile(TABLE);
	ForceAssert(outfile.good());
	auto table = output.GetTable();
	for ( size_t irow = 2; irow < table.size(); ++irow )
	    for ( auto itr = table[irow].begin(); itr != table[irow].end(); ++itr ) {
		outfile << *itr;
		if ( itr+1 != table[irow].end() ) outfile << " "; else outfile << endl;
	    }
    }

    if ( CCR ) {
	TextTable t;
	t << "h-length" << Tab << "del(1)" << Tab << "equal" << Tab << "ins(1)" << Tab << "del(1)-rate" << Tab << "ins(1)-rate" << EndRow;
	t << SingleLine;
	for ( size_t hlen = crates.Min(); hlen <= crates.Max(); ++hlen ) {
	    double sum = crates[hlen][0] + crates[hlen][1] + crates[hlen][2]; ForceAssertEq(crates[hlen].size(),3U);
	    t << hlen << Tab << crates[hlen][0] << Tab << crates[hlen][1] << Tab << crates[hlen][2];
	    t << Tab << crates[hlen][0]/ sum << Tab << crates[hlen][2] / sum << EndRow;
	}
	t.Print(std::cout, 1);
    }

    cout << Date() << ": about to exit" << endl;

    return 0;
}

