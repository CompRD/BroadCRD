/*
 * BridgeTools.cc
 *
 *  Created on: Dec 20, 2013
 *      Author: blau
 */

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "pacbio/BridgeTools.h"


namespace pacbio_bridge_tools{
kmer_val_t KmerVal(const basevector&b, size_t K, size_t p) {
    ForceAssert(K <= 32);
    ForceAssert(p + K <= b.size());
    kmer_val_t out = 0;
    for (size_t ii = 0; ii < K; ++ii) {
        out <<= 2;
        out += b[p + ii];
    }
    return out;
}

unsigned int SmithWatAffine_loc( const basevector& S, const basevector& T, alignment& a
                               , SmithWatAffineResource_t& resource
                               , bool penalize_left_gap, bool penalize_right_gap, const int mismatch_penalty, const int gap_open_penalty, const int gap_extend_penalty )
{
     const int Infinity = 100000000;
     ForceAssertGt( S.size(), 0u );
     ForceAssertGt( T.size(), 0u );

     unsigned int n = S.size( ), N = T.size( );

     //     ForceAssertLe( n, N );

     avector<char> s, t;
     s.resize(n);
     for ( unsigned int i = 0; i < n; i++ )
          s(i) = S[i];
     t.resize(N);
     for ( unsigned int i = 0; i < N; i++ )
          t(i) = T[i];

     int best_score = Infinity;
     resource.reset(n+1,N+1);

     resource.score_x(0,0) = 0;
     resource.score_y(0,0) = Infinity;
     resource.score_z(0,0) = Infinity;
     resource.x_from(0,0) = 's';
     resource.y_from(0,0) = 's';
     resource.z_from(0,0) = 's';

     for ( unsigned int i = 1; i <= n; i++ )
     {    resource.score_x(i,0) = Infinity;
          resource.score_y(i,0) = Infinity;
          resource.score_z(i,0) = gap_open_penalty + gap_extend_penalty * i;
          resource.x_from(i,0) = 's';
          resource.y_from(i,0) = 's';
          resource.z_from(i,0) = 's';   }

     for ( unsigned int j = 1; j <= N; j++)
       {  resource.score_x(0,j) = Infinity;
          resource.score_y(0,j) = (penalize_left_gap ? gap_open_penalty + gap_extend_penalty * j : 0);
          resource.score_z(0,j) = Infinity;
          resource.x_from(0,j) = 's';
          resource.y_from(0,j) = 's';
          resource.z_from(0,j) = 's';   }

     for ( unsigned int i = 1; i <= n; i++ )
     {   for ( unsigned int j = 1; j <= N; j++ )
          {    unsigned int x_x = resource.score_x(i-1,j-1) + mismatch_penalty * ( s(i-1) != t(j-1) );
               unsigned int x_y = resource.score_y(i-1,j-1) + mismatch_penalty * ( s(i-1) != t(j-1) );
               unsigned int x_z = resource.score_z(i-1,j-1) + mismatch_penalty * ( s(i-1) != t(j-1) );
               unsigned int y_x = resource.score_x(i,j-1) + (i != n || penalize_right_gap ? gap_open_penalty : 0);
               unsigned int y_y = resource.score_y(i,j-1) + (i != n || penalize_right_gap ? gap_extend_penalty : 0);
               unsigned int y_z = Infinity; //score_z[i][j-1] + gap_open_penalty;
               unsigned int z_x = resource.score_x(i-1,j) + gap_open_penalty;
               unsigned int z_y = Infinity; //score_y[i-1][j] + gap_open_penalty;
               unsigned int z_z = resource.score_z(i-1,j) + gap_extend_penalty;

               resource.score_x(i,j) = Min( Min( x_x, x_y ), x_z );
               resource.score_y(i,j) = Min( Min( y_x, y_y ), y_z );
               resource.score_z(i,j) = Min( Min( z_x, z_y ), z_z );

               if ( x_x <= x_y )
               {    if ( x_x <= x_z ) resource.x_from(i,j) = 'x';
                    else resource.x_from(i,j) =  'z';    }
               else
               {    if ( x_y <= x_z ) resource.x_from(i,j) = 'y';
                    else resource.x_from(i,j) =  'z';    }

               if ( y_x <= y_y )
               {    if ( y_x <= y_z ) resource.y_from(i,j) = 'x';
                    else resource.y_from(i,j) =  'z';    }
               else
               {    if ( y_y <= y_z ) resource.y_from(i,j) = 'y';
                    else resource.y_from(i,j) =  'z';    }

               if ( z_x <= z_y )
               {    if ( z_x <= z_z ) resource.z_from(i,j) = 'x';
                    else resource.z_from(i,j) =  'z';    }
               else
               {    if ( z_y <= z_z ) resource.z_from(i,j) = 'y';
                    else resource.z_from(i,j) =  'z';    }    }    }

     best_score = Min( resource.score_x(n,N), Min( resource.score_y(n,N), resource.score_z(n,N) ) );

//     vec< vec<unsigned char> > *from;
     vec<unsigned char> *from;
     if ( resource.score_x(n,N) <= resource.score_y(n,N) )
     {    if ( resource.score_x(n,N) <= resource.score_z(n,N) ) from = &resource.x_from();
          else from = &resource.z_from();    }
     else
     {    if ( resource.score_y(n,N) <= resource.score_z(n,N) ) from = &resource.y_from();
          else from = &resource.z_from();    }

     int i = int(n);
     int j = int(N);
     int lcount = 0, g1count = 0, g2count = 0;
     int last_length = 0;
     avector<int> gaps(0), lengths(0);
     while(1)
     {
//          unsigned char dir = (*from)[i][j];
          unsigned char dir = (*from)[i*resource.n2()+j];
          //cout << dir;
          if ( from == &resource.x_from() )
          {    if ( g1count > 0 )
               {    if ( last_length > 0 )
                    {    gaps.Prepend( g1count );
                         lengths.Prepend( last_length );    }
                    g1count = 0;    }
               if ( g2count > 0 )
               {    if ( last_length > 0 )
                    {    gaps.Prepend( -g2count );
                         lengths.Prepend( last_length );    }
                    g2count = 0;    }
               ++lcount;
               --i;
               --j;   }
          else if ( from == &resource.z_from() )  // gap on long sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g1count == 0 );
               ++g2count;
               --i;    }
          else                           // gap on short sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g2count == 0 );
               ++g1count;
               --j;    }

          if ( dir == 'x') from = &resource.x_from();
          else if ( dir == 'y') from = &resource.y_from();
          else from = &resource.z_from();

          if( (*from)[i*resource.n2()+j] == 's' ) break;

    }

     //cout << "\n";

     if ( g1count != 0 ) gaps.Prepend( g1count );
     else if ( g2count != 0 ) gaps.Prepend( -g2count );
     else gaps.Prepend(0);

     lengths.Prepend( lcount );

     int pos1 = i;
     int pos2 = j;

     if ( gaps(0) < 0 )
     {   pos2 -= gaps(0);
         gaps(0) = 0;    }

     if ( gaps(0) > 0 )
     {   pos1 += gaps(0);
         gaps(0) = 0;    }

     int errors = best_score;
     a = alignment( pos1, pos2, errors, gaps, lengths );

     return best_score;    }

void pos1_to_pos2(vec<size_t>& map, const alignment&a){
    map.clear();
    int pos1, pos2, errors;
    avector<int> gaps, lengths;
    a.Unpack( pos1, pos2, errors, gaps, lengths );
    ForceAssert(lengths.length==gaps.length);
    const size_t nL=lengths.length;
    ForceAssert(gaps.x[0]==0);

    size_t longest=0;
    for(size_t ll=0;ll<nL;++ll){
        longest+=lengths.x[ll];
        longest+=abs(gaps.x[ll]);
    }
    map.resize(longest);

    size_t next1=0;
    size_t next2=pos2;

    for(;next1<size_t(pos1);++next1){
        map[next1]=next2;
    }
    for(size_t ll=0;ll<nL;++ll){
        if(gaps.x[ll] >0 ){
            next2 += gaps.x[ll];
        }
        else if(gaps.x[ll]<0){
            for(size_t aa=0,bb=-gaps.x[ll];aa<bb;++aa){
                map[next1++]=next2;
            }
        }
        for(size_t aa=0,bb=lengths.x[ll];aa<bb;++aa){
            map[next1++]=next2++;
        }
    }
    map.resize(next1);
}
void alternatives_t::pileup(const vecbasevector& reads,const vecqualvector& quals, const vec<vec<size_t>>& pos_maps
                           ,const typename qualvector::value_type threshold, unsigned int flank
                           ,int64_t pos_dev){
    pos_dev=abs(pos_dev);
    std::cout << Date() << ": Using kmerized pile-up strategy" << std::endl;
    ForceAssert(reads.size()==quals.size());
    ForceAssert(reads.size()==pos_maps.size());
    std::unordered_map<size_t,kp_analysis_t> k_kpa;
    for(const auto& elem_t: mElements){
        size_t k=elem_t.qseq.size();
        auto itr = k_kpa.find(k);
        if(itr==k_kpa.end()) k_kpa.insert(make_pair(k,kp_analysis_t(k)));

        k=elem_t.rseq.size();
        itr = k_kpa.find(k);
        if(itr==k_kpa.end()) k_kpa.insert(make_pair(k,kp_analysis_t(k)));
    }
    std::cout << Date() << ": kmerizing " << k_kpa.size() << " different k values"<< std::endl;
    #pragma omp parallel for schedule(dynamic,1)
    for(size_t kk=0;kk<k_kpa.size();++kk) {
        auto itr = k_kpa.begin();
        for(size_t kkk=0;kkk<kk;++kkk,++itr){};
        #pragma omp critical
        {
            std::cout << Date() << ": kmerizing K=" << (*itr).second.K() << std::endl;
        }
        (*itr).second.processReads(reads,quals,pos_maps,threshold,flank);
    }

    auto getCounts=[](const basevector&seq,const int64_t front, const int64_t back, const kp_analysis_t& kpa){
        ForceAssert(seq.size()==kpa.K());
        const auto& kpc = kpa.getKPC();
        size_t out=0;
        auto itr = kpc.find(KmerVal(seq,seq.size()));
        if( itr==kpc.end()) return out;
        for(const auto& p_c: (*itr).second){
            if( front <= p_c.first && p_c.first <= back){
                out+=p_c.second;
            }
        }
        return out;
    };

    std::cout << Date() << ": tabulating votes" << std::endl;
    for(auto& elem_t: mElements){
        auto itr =k_kpa.find(elem_t.qseq.size());
        ForceAssert( itr != k_kpa.end());
        elem_t.qcount += getCounts(elem_t.qseq,elem_t.rfront-pos_dev,elem_t.rback+pos_dev,(*itr).second);
        itr =k_kpa.find(elem_t.rseq.size());
        ForceAssert( itr != k_kpa.end());
        elem_t.rcount += getCounts(elem_t.rseq,elem_t.rfront-pos_dev,elem_t.rback+pos_dev,(*itr).second);
    }
    std::cout << Date() << ": Done" << std::endl;
}

}
