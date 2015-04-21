///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Jul 28, 2014 - <crdhelp@broadinstitute.org>
//

#ifndef LINEMAPREADER_H_
#define LINEMAPREADER_H_
#include "FastIfstream.h"
#include "feudal/CharString.h"


class LineMapReader {
public:
    struct LineLocus {
        String chr;
        int start;
        int stop;
        Bool ref_gap;
        int lineId;
        int edgeId1;
        int edgeId2;
        int len;
        double cov;
        friend std::ostream& operator<<(std::ostream& out, LineLocus const& self) {
            out << self.chr << " " <<
                    self.start << " " << self.stop << " ";
            if ( self.ref_gap ) out << "** REF GAP **";
            else {
                out << "line[" << self.lineId << "] " <<
                        self.edgeId1 << ".." << self.edgeId2 << " len=" << self.len;
                if ( self.cov >= 0 ) cout << " cov=" << self.cov;
            }
            return out;
        }
    };
    LineMapReader() = delete;
    explicit LineMapReader(String const& filename) {
        fast_ifstream zin(filename);
        String line;
        vec<String> tokens;
        while (1) {
            LineLocus tmp;
            getline( zin, line );
            if ( zin.fail() ) break;    // ick, ick
            if ( line == "" ) continue;
            tokens.clear();
            Tokenize(line, tokens);
            if ( tokens.size() < 3  ) {
                cout << "weird line (first kind) read in " << filename <<
                        ": " << line << endl;
                continue;
            }
            tmp.chr = tokens[0];
            tmp.start = RemoveCommas(tokens[1]).Int();
            tmp.stop = RemoveCommas(tokens[2]).Int();
            if (line.Contains("**REF GAP**")) {
                tmp.ref_gap = True;
            } else {
                if ( tokens.size() < 6u ) {
                    cout << "weird line read in " << filename <<
                            ": " << line << endl;
                    continue;
                }
                tmp.ref_gap = False;
                tmp.lineId = tokens[3].Between("line[","]").Int();
                tmp.edgeId1 = tokens[4].Before("..").Int();
                tmp.edgeId2 = tokens[4].After("..").Int();
                tmp.len = tokens[5].After("len=").Int();
                if ( tokens.size() > 6 ) tmp.cov = tokens[6].After("cov=").Double();
                else tmp.cov = -1.;
            }
            mLineLoci.push_back(tmp);
        }
    }

    void Dump() {
        for ( auto const& line: mLineLoci  )
            cout << line << endl;
    }

    vec<LineLocus> const& LineLoci() { return mLineLoci; }

private:
    vec<LineLocus> mLineLoci;
};


#endif /* LINEMAPREADER_H_ */
