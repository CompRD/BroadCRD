///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"

#ifndef BNG_EDGE_H
#define BNG_EDGE_H

class bng_edge {

     public:

     bng_edge( ) { }
     bng_edge( const int len, const vec< pair<int,int> >& id_pos ) 
          : len_(len), id_pos_(id_pos) { }

     int Len( ) const { return len_; }
     void SetLen( const int len ) { len_ = len; }

     int Weight( ) const { return id_pos_.size( ); }

     const vec< pair<int,int> >& IdPos( ) const { return id_pos_; }
     void SetIdPos( const vec< pair<int,int> >& id_pos ) { id_pos_ = id_pos; }

     vec<int> Rids( ) const
     {    vec<int> x;
          for ( int i = 0; i < id_pos_.isize( ); i++ )
               x.push_back( id_pos_[i].first );
          return x;    }

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     private:

     int len_;
     vec< pair<int,int> > id_pos_;

};

SELF_SERIALIZABLE(bng_edge);

#endif
