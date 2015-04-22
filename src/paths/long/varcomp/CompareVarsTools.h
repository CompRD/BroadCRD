///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef COMPARE_VARS_TOOLS_H
#define COMPARE_VARS_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "paths/long/fosmid/Fosmids.h"

class event {

     public:

     event( const String& alias, const String& chr, const int pos, 
          const String& ref, const String& alt, const String& genotype,
          const String& filter )
          : alias(alias), chr(chr), pos(pos), ref(ref), alt(alt), 
          genotype(genotype), filter(filter) { }

     String alias;
     String chr;
     int pos;
     String ref;
     String alt;
     String genotype;
     String filter;

     friend Bool operator<( const event& e1, const event& e2 )
     {    if ( e1.chr < e2.chr ) return True;
          if ( e1.chr > e2.chr ) return False;
          if ( e1.pos < e2.pos ) return True;
          if ( e1.pos > e2.pos ) return False;
          if ( e1.ref < e2.ref ) return True;
          if ( e1.ref > e2.ref ) return False;
          if ( e1.alt < e2.alt ) return True;
          if ( e1.alt > e2.alt ) return False;
          if ( e1.genotype < e2.genotype ) return True;
          if ( e1.genotype > e2.genotype ) return False;
          if ( e1.filter < e2.filter ) return True;
          if ( e1.filter > e2.filter ) return False;
          return e1.alias < e2.alias;    }

};

// A disk_map is a map from String to int that is maintained on disk and is in
// principle not corruptable.  It is specialized to the particular needs of
// CompareVars.  The String is not allowed to contain "=" or ";".  The map is
// stored on disk as x1=y1;x2=y2;...;xn=yn.

class disk_map {

     public:

     disk_map( const String& fn )
     {    fd_ = open( fn.c_str( ), O_RDWR | O_CREAT, 0775 );
          int64_t total = lseek( fd_, 0, SEEK_END );
          if ( total == 0 ) 
          {    nbytes_ = 0;
               write( fd_, &nbytes_, 8 );
               fsync(fd_);    }
          lseek( fd_, 0, SEEK_SET );
          read( fd_, &nbytes_, 8 );
          String contents;
          contents.resize(nbytes_);
          if ( nbytes_ > 0 ) read( fd_, &contents[0], nbytes_ );
          vec<String> fields;
          Tokenize( contents, {';'}, fields );
          for ( int j = 0; j < fields.isize( ); j++ )
               val_[ fields[j].Before( "=" ) ] = fields[j].After( "=" ).Int( );    }

     ~disk_map( )
     {    fsync(fd_);
          lseek( fd_, 0, SEEK_SET );
          write( fd_, &nbytes_, 8 );
          fsync(fd_);
          close(fd_);    }

     Bool Defined( const String& s )
     {    return val_.find(s) != val_.end( );    }

     int Value( const String& s )
     {    return val_[s];    }

     void Set( const String& s, const int x )
     {    ForceAssert( !Defined(s) );
          val_[s] = x;
          String eq = s + "=" + ToString(x);
          if ( nbytes_ > 0 ) eq = ";" + eq;
          lseek( fd_, nbytes_ + 8, SEEK_SET );
          write( fd_, eq.c_str( ), eq.size( ) );
          nbytes_ += eq.size( );    }

     private:

     int fd_;
     int64_t nbytes_;
     map<String,int> val_;

};

// Genotype map for the CEPH family.  The 17 member CEPH family has been deeply
// sequenced by Ilumina.  The family consists of NA12878 and her husband, their 
// four parents, and their eleven children (five daughters and six sons).  There
// are bams and vcfs for each of the 17 members.  From these data we deduce a 
// partial map of the family, which we describe.  In the map, the genome is divided 
// into intervals, with some gaps between them.  At a given point on the genome,
// a given child inherits one of mom's chromosomes and one of dad's chromosomes,
// so there are four possibilities.  The map asserts that on a given interval,
// for each child, the inheritance pattern is constant across the interval.  The
// map specifies a partition, consisting of four subsets of {d1,...,d5,s1,...,s6}, 
// however the map does not completely order these subsets.  These sets are 
// {{x1,x2},{y1,y2}}, where the inner sets and outer set are unordered.  

int GenotypeId( const vec<int>& s1, const vec<int>& s2, const vec<int>& s3 );

void GenotypeGroups( vec<vec<vec<int>>>& groups, vec<vec<int>>& compat );

void BuildGenotypeMap( const String& root );

void DefineGenotypeMap( const String& root,
     vec< pair< triple<String,int,int>, vec<int> > >& gmap,
     vec< vec< vec<int> > >& groups, vec< vec<int> >& compat );

int Search( disk_map& search_db, vec<vecbasevector>& reads,
     const vec<String>& reads_id, const String& home, const String& query,
     const String& dataset );

void RunFunctions( const String& FUN, const String& R, const String& home,
     const String& ref_dir, const String& loc, const String& ID,
     const basevector& query, const basevector& target, const String& fos_fasta, 
     const int rstart, const int rstop );

void RunSome( const String& root, const String& ID, String& com, 
     const Bool COUNT_ONLY, const Bool parallel );

#endif
