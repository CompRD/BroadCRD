///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Nov 19, 2014 - <crdhelp@broadinstitute.org>
//

#include "paths/long/magic/SimpleDataFrame.h"
#include "feudal/CharString.h"
#include "system/System.h"
#include "Vec.h"

void SimpleDataFrame::ReadFromTextFile(String const& filename, String const& colformat)
{
    formats = StrSplit(colformat);
    std::ifstream inFile( filename, std::ios_base::in );
    if (inFile.bad()) { FatalErr("failed to open file " + filename); }

    std::string line;
    if ( !std::getline(inFile,line) || line.size() < 1 ) { FatalErr("error reading first line of "+filename); }
    if ( line[0] != '#' ) { FatalErr("no header in file " + filename); }

    if ( line == "# ATTRS" ) {
        this->PullAttributes(inFile);
        if ( !std::getline(inFile,line) || line.size() < 1 ) { FatalErr("error reading first line of head in " + filename); }
    }

    vec<String> header;
    if ( line == "# HEAD" ) {
        header = this->PullHead(inFile);
        if ( header.size() != formats.size() ) {
            FatalErr("header defined " + ToString(header.size()) +
                    " columns, but formats specified " +ToString(formats.size()) +
                    " in " + filename);
        }
        if ( !std::getline(inFile,line) || line.size() < 1 ) { FatalErr("error reading first line of head in " + filename); }
    } else FatalErr("no # HEAD in " + filename);

    if ( line == "# DATA" )
        this->PullData(inFile, header);
    else
        FatalErr("no # DATA in " + filename);
}

void SimpleDataFrame::PullAttributes(std::ifstream& inFile)
{
    std::string buf;

    while ( true ) {
        if ( inFile.peek() == '#' ) break;
        if ( !std::getline(inFile,buf) || buf.size() < 1 ) { FatalErr("Died during attributes"); }
        String line(buf);
        if ( !line.Contains(":") ) { FatalErr("attribute didn't have a colon"); }
        auto name = line.Before(":");
        double val = line.After(":").Double();

        this->attr[name] = PolyVal<double>::New();
        this->attr[name]->set_double(val);
        attr_names.push_back(name);
    }
}


vec<String> SimpleDataFrame::PullHead(std::ifstream& inFile)
{
    std::string buf;

    if ( !std::getline(inFile, buf) || buf.size() < 1 )  { FatalErr("Died during header"); }
    names = SimpleDataFrame::StrSplit(buf, ',');
    return names;
}

void SimpleDataFrame::PullData(std::ifstream& inFile, vec<String> const& header )
{
    std::string buf;
    ForceAssertEq(header.size(), formats.size());
    while ( std::getline(inFile, buf) && buf.size() ) {
        vec<String> values = StrSplit(buf, ',');
        if ( values.size() != header.size() ) {
            FatalErr("short line with only " + ToString(values.size()) +
                    " values, expected " + ToString( header.size() ) +
                    " at DATA line " + ToString(nrows));
        }
        for ( size_t i = 0; i < values.size(); ++i ) {
            auto const& name = header[i];
            auto format = ToLower(formats[i]);
            if ( nrows == 0 ) {
                if ( format == "s" )
                    this->cols[name] = PolyVal<vec<String>>::New();
                else if ( format == "f" )
                    this->cols[name] = PolyVal<vec<double>>::New();
                else if ( format == "i" )
                    this->cols[name] = PolyVal<vec<int>>::New();
                else
                    FatalErr("bad format");
            }
            this->cols[name]->push_back(values[i]);
        }
        nrows++;
    }
}

void SimpleDataFrame::Dump() const
{
    // output attributes
    if ( this->attr.size() )  {
        cout << "# ATTRS" << endl;
        for ( auto& kv : this->attr ) {
            cout << kv.first << ": " << kv.second->get_double() << endl;
        }
    }
    // output data columns
    if ( this->names.size() ) {
        cout << "# HEAD" << endl;
        for ( auto const& name : names ) {
            cout << name;
            if ( &name != &names.back() ) cout << ",";
        }
        cout << endl;
        cout << "# DATA" << endl;
        for ( size_t row = 0; row < nrows; ++row ) {
            cout << row << ": ";
            for ( size_t col = 0; col < names.size(); ++col ) {
                auto format = ToLower(formats[col]);
                auto const& name = names[col];
                if ( formats[col] == "f" )
                    cout << cols.at(name)->get_vec_double()[row];
                else if ( formats[col] == "s")
                    cout << cols.at(name)->get_vec_String()[row];
                else if ( formats[col] == "i")
                    cout << cols.at(name)->get_vec_int()[row];
                else
                    FatalErr("wrong format");
                if (&name != &names.back()) cout << ",";
            }
            cout << endl;
        }

    }
}

template <>
void PolyVal<vec<double>>::push_back(const String& s ) {
    double v = s.Double();
    val.push_back(v);
}

template <>
void PolyVal<vec<String>>::push_back(const String& s) {
    val.push_back(s);
}

template <>
void PolyVal<vec<double>>::push_back(const double v) {
    val.push_back(v);
}

template <>
void PolyVal<vec<int>>::push_back(const String& s) {
    val.push_back(s.Int());
}

template <>
void PolyVal<vec<int>>::push_back(const int v) {
    val.push_back(v);
}

template <class T>
void PolyVal<T>::push_back(const String& s)
{
    FatalErr("tried to push back a value on a non-vec type");
}

template <class T>
void PolyVal<T>::push_back(const double v)
{
    FatalErr("tried to push back a value on a non-vec type");
}

template <class T>
void PolyVal<T>::push_back(const int v)
{
    FatalErr("tried to push back a value on a non-vec type");
}
