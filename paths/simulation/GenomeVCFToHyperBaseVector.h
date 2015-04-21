///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef GENOME_VCF_TO_HYPERBASEVECTOR_H
#define GENOME_VCF_TO_HYPERBASEVECTOR_H

#include "CoreTools.h"
#include "paths/long/fosmid/FosmidPool.h"

//---- Range class - parse ranges such as chrX:10-100000
//
// This class exists to cluster together some of the nastiness related to parsing ranges and Longproto
// fosmid numbers into one place.  I also wanted a consistent interface that made clear whether a
// particular index was zero-based or one-based, hence methods such as start0 and start1.

struct Range {
    explicit Range(string str, const vecString& names) {

	if ( str.find_first_of(':') == string::npos )
	    str = load_fosmid_region(str);

	std::istringstream s(str);

	if ( names.size() > 0 ) {
	    // names specified, so parse the front of the range as a string
	    // and find it in the names list.
	    char c = '\0';

	    s >> c;
	    while ( c != ':' ) {
		_chr.push_back(c);
		s >> c;
	    }

	    if ( c != ':' )
		FatalErr("no : in range");

	    for ( _chr0 = 0; _chr0 < names.size() && names[_chr0] != _chr; _chr0++ )
		/* do nothing */ ;

	    if ( _chr0 >= names.size() )
		FatalErr("Failed to find contig name " << _chr << " in file ");
	} else {
	    // names not specified, so parse the front of the range as an int
	    // and make it the zero-based contig index, say for a naked .FASTB with
	    // no associated .names file.
	    s >> _chr0;
	    if ( s.get() != ':' )  // if no :, assume fosmid pool region.
		FatalErr("no : in range.");

	}

	s >> _start0;
	if ( s.get() != '-' )
	    FatalErr("Failed to find dash in range string; peek showed: " << s.peek() );
	s >> _end0;
    }

    string load_fosmid_region(const string& s) {
	std::istringstream ss(s);
	ss >> _chr0;
	return load_fosmid_region(_chr0);
    }

    string load_fosmid_region( size_t chr0 ) {
	vec<String> regions;
	typedef vec< vec< pair<String,String> > > Stuff;
	Stuff stuff1, stuff2, stuff3;
	ParseFosmidPoolMetainfo( regions, stuff1, stuff2, stuff3 );
	ForceAssertLt( chr0, regions.size() );
	return regions[chr0];
    }

    const String dumb_chr_mapping() const {
	std::string pref = "";
	std::ostringstream s;
	if ( _chr0 < 22 )
	    s << pref << _chr0+1;
	else if ( _chr0 == 22 )
	    s << pref << "X";
	else if ( _chr0 == 23 )
	    s << pref << "Y";
	else if ( _chr0 == 24 )
	    s << pref << "M";
	return s.str();
    }

    const String chr() const {
	if ( _chr != "" ) return _chr;
	else return this->dumb_chr_mapping();
    };

    size_t chr0() const { return _chr0; };
    size_t start0() const { return _start0; };
    size_t end0() const { return _end0; };

    size_t start1() const { return _start0+1;};       // half-open interval [start1,end1)
    size_t end1() const { return _end0+1;};

    void extend_by( const int ext ) { _start0 -= ext; _end0 += ext; }

private:
    String _chr;
    size_t _chr0;       // chromosome index 0-based
    size_t _start0;     // half-open interval [start0,end0)
    size_t _end0;
};

#endif
