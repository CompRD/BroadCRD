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

#ifndef SIMPLEDATAFRAME_H_
#define SIMPLEDATAFRAME_H_

#include <iosfwd>
#include <iostream>
#include <typeinfo>
#include "Vec.h"
#include "feudal/CharString.h"
#include "system/System.h"



class PolyValBase {
public:
    PolyValBase() { };
    virtual ~PolyValBase() {};

    virtual double& get_double() = 0;
    virtual String& get_String() = 0;
    virtual int& get_int() = 0;

    void set_double(const double d) { this->get_double() = d; }
    void set_String(const String& s) { this->get_String() = s; }
    void set_int(const int i) { this->get_int() = i; }

    virtual vec<double>& get_vec_double() = 0;
    virtual vec<String>& get_vec_String() = 0;
    virtual vec<int>& get_vec_int() = 0;

    virtual void push_back(const String& s) = 0;
    virtual void push_back(const double v) = 0;
    virtual void push_back(const int v) = 0;
    virtual void writeBinary(BinaryWriter& writer) const = 0;
    virtual void readBinary(BinaryReader& reader) = 0;
    static size_t externalSizeof() { return 0u; }

    virtual std::shared_ptr<PolyValBase> Clone() = 0;
};
SELF_SERIALIZABLE(PolyValBase);

template <class T> class PolyVal : public PolyValBase {
public:
    typedef PolyVal<T> MyType;
    PolyVal(T const& v) : val(v) {};
    PolyVal() {};
    ~PolyVal() {}

    T& get() {
        return val;
    }

    void writeBinary(BinaryWriter& writer) const {
        writer.write(val);
    }

    void readBinary(BinaryReader& reader) {
        reader.read(&val);
    }

    // this could be done more generically, but fine for now.
    double& get_double() {
        ForceAssertEq(typeid(T).hash_code(), typeid(double).hash_code());
        return reinterpret_cast<PolyVal<double>*>(this)->get();
    }

    String& get_String() {
        ForceAssertEq(typeid(T).hash_code(), typeid(String).hash_code());
        return reinterpret_cast<PolyVal<String>*>(this)->get();
    }

    int& get_int() {
        ForceAssertEq(typeid(T).hash_code(), typeid(int).hash_code());
        return reinterpret_cast<PolyVal<int>*>(this)->get();
    }

    vec<double>& get_vec_double() {
        ForceAssertEq(typeid(T).hash_code(), typeid(vec<double>).hash_code() );
        return reinterpret_cast<PolyVal<vec<double>>* >(this)->get();
    }

    vec<String>& get_vec_String() {
        ForceAssertEq(typeid(T).hash_code(), typeid(vec<String>).hash_code() );
        return reinterpret_cast<PolyVal<vec<String>>* >(this)->get();
    }

    vec<int>& get_vec_int() {
        ForceAssertEq(typeid(T).hash_code(), typeid(vec<int>).hash_code() );
        return reinterpret_cast<PolyVal<vec<int>>* >(this)->get();
    }

    void push_back(const String& s);
    void push_back(const double v);
    void push_back(const int v);

    static std::shared_ptr<MyType> New() {
        return std::shared_ptr<MyType>(new MyType); }

    std::shared_ptr<PolyValBase> Clone() {
        auto cl = this->New();
        *cl = *this;
        return cl;
    }

private:
    T val;
};


class SimpleDataFrame {
public:
    SimpleDataFrame(String const& filename, String const& colformat,
            Bool debugme = False) : nrows(0u), debug(debugme) {
        ReadFromTextFile(filename,colformat);
    }
    SimpleDataFrame(SimpleDataFrame const& src) { this->Duplicate(src); }
    SimpleDataFrame& operator=(SimpleDataFrame const& src) {
        this->Duplicate(src); return *this;
    };
    SimpleDataFrame() : nrows(0u), debug(False) {};
    ~SimpleDataFrame() {};

    void Dump() const;

    void Duplicate(SimpleDataFrame const& src) {
        for ( auto const& pair_attr : src.attr )
            attr[pair_attr.first] = pair_attr.second->Clone();
        for ( auto const& pair_col : src.cols )
            cols[pair_col.first] = pair_col.second->Clone();
        this->formats = src.formats;
        this->names = src.names;
        this->attr_names = src.attr_names;
        this->nrows = src.nrows;
        this->debug = src.debug;
    }

    void writeBinary(BinaryWriter& writer) const {
        writer.write(formats);
        writer.write(names);
        writer.write(attr_names);
        writer.write(nrows);
        for ( auto aname : attr_names )
            attr.at(aname)->writeBinary(writer);
        for ( auto name : names )
            cols.at(name)->writeBinary(writer);
    }
    void readBinary(BinaryReader& reader) {
        reader.read(&formats);
        reader.read(&names);
        reader.read(&attr_names);
        reader.read(&nrows);
        for ( auto aname : attr_names ) // allocate attributes
            this->attr[aname] = PolyVal<double>::New();
        for ( auto aname : attr_names ) // read
            this->attr.at(aname)->readBinary(reader);
        ForceAssertEq(names.size(), formats.size());
        for ( size_t i = 0; i < names.size(); ++i ) {
            auto format = ToLower(formats[i]);
            auto name = names[i];
            if ( format == "f" )
                this->cols[name] = PolyVal<vec<double>>::New();
            else if ( format == "s" )
                this->cols[name] = PolyVal<vec<String>>::New();
            else if ( format == "i" )
                this->cols[name] = PolyVal<vec<int>>::New();
            this->cols[name]->readBinary(reader);
        }
    }
    static size_t externalSizeOf() { return 0u; }

    // columns are vecs, so swapping rows is hard
    void SwapRows(size_t i, size_t j) {
        ForceAssertLt(i, nrows);
        ForceAssertLt(j, nrows);
        for ( size_t icol = 0; icol < names.size(); ++icol ) {
            auto format = ToLower(formats[icol]);
            auto name = names[icol];
            if ( format == "f" )
                std::swap(this->cols[name]->get_vec_double()[i],
                          this->cols[name]->get_vec_double()[j]);
            else if ( format == "s" )
                std::swap(this->cols[name]->get_vec_String()[i],
                          this->cols[name]->get_vec_String()[j]);
            else if ( format == "i" )
                std::swap(this->cols[name]->get_vec_int()[i],
                          this->cols[name]->get_vec_int()[j]);
        }
    }

    std::shared_ptr<PolyValBase> operator[](std::string const& s) {
        try {
            return cols.at(s);
        } catch ( std::out_of_range& e ) {
            cerr << "error retrieving column " << s << " from a DataFrame. Possible choices are:" << endl;
            for ( auto const& p : cols ) cout << p.first << endl;
            FatalErr("exiting...");
        }
    }
    std::shared_ptr<PolyValBase> operator()(std::string const& s) {
        try {
            return attr.at(s);
        } catch ( std::out_of_range& e ) {
            cerr << "error retrieving attr " << s << " from a DataFrame" << endl;
            for ( auto const& p : attr ) cout << p.first << endl;
            FatalErr("exiting...");
        }
    }

    vec<String> const& Names() const { return names; }
    size_t size() const {  return nrows; }

private:
    std::map<std::string,std::shared_ptr<PolyValBase>> attr;
    std::map<std::string,std::shared_ptr<PolyValBase>> cols;
    vec<String> formats;
    vec<String> names;
    vec<String> attr_names;
    size_t nrows;
    bool debug;

    void ReadFromTextFile(String const& filename, String const& colformat);

    void PullAttributes(std::ifstream& inFile);
    vec<String> PullHead(std::ifstream& inFile);
    void PullData(std::ifstream& inFile, vec<String> const& header);

    // stolen from my code in VCF.cc; I should put it in the string fcns
    // split - tokenize a string based on a separator and return a vector of strings representing the
    // characters between each instance of a separator.  If there is no separator, then the vector will
    // have one entry equal to the original string.
    static vec<String> StrSplit( const String& input, const char sep=',' )
    {
        vec<String> tokens;

        for (auto i = input.begin(); i != input.end(); /* */ ) {
            tokens.push_back(String());

            while ( i != input.end() && *i != sep )
                tokens.back().push_back(*i++);

            if ( i != input.end() ) i++;
        }

        return tokens;
    }
};
SELF_SERIALIZABLE(SimpleDataFrame);



#endif /* SIMPLEDATAFRAME_H_ */
