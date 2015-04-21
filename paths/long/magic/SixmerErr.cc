// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
//
// throw-away tool by neilw to count how many "good" and "bad" sixmers
// in a set of reads.  Good means that it falls fully within a non-gap
// portion of an alignment and all bases match.

#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "kmers/KmerRecord.h"



int main(int argc, char* argv[])
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String(READS_FASTB);
    CommandArgument_String(REF_FASTB);
    CommandArgument_String(QLTOUT);
    CommandArgument_Int(KMER_THRESH);
    EndCommandArguments;

    vecbasevector bases(READS_FASTB);
    vecbasevector ref(REF_FASTB);
    basevector const& ref0 = ref[0];
    vec<look_align> aligns;
    vec<vec<int>> aligns_index;
    LoadLookAligns(QLTOUT, aligns, aligns_index, bases.size() );

    ForceAssertEq( aligns.size(), bases.size() );



    // for each read, walk the read and the reference via alignment
    // each increment on the reference side pushes a new six-mer.
    // if it's not in a LEN run, then it's an error six-mer.
    const int K = 6;
    typedef kmer<K> KmerType;
    typedef triple<KmerType,int,int> KmerStat;  // kmer, bads, goods
    vec<KmerStat> pile1;

    for ( size_t id = 0; id < aligns.size(); ++id ) {
        ForceAssertEq(aligns[id].target_id, 0);
        align const& x = aligns[id].a;
        basevector const& read = bases[id];
        int p1 = x.pos1( ), p2 = x.pos2( );

        for ( int j = 0; j < x.Nblocks(); j++ ) {
            KmerType ktmp;

            if ( x.Gaps(j) > 0 ) {
                for ( int l = 0; l < x.Gaps(j); ++l ) {
                    ktmp.SetToSubOf( ref0, p2+l );
                    pile1.push( ktmp, 1, 0 );   // bad
                    p2++;
                }
            } else if ( x.Gaps(j) < 0 ) {
                p1 -= x.Gaps(j);
            }

            for (int l = 0; l < x.Lengths(j); ++l ) {
                ktmp.SetToSubOf( ref0, p2+l );
                if ( l < x.Lengths(j)-K ) {
                    // inefficient
                    bool good = true;
                    for ( int ll = 0; good && ll < K; ++ll )
                        if ( read[p1+l+ll] != ref0[p2+l+ll] ) good = false;
                    if ( good ) pile1.push(ktmp, 0, 1 );    // good
                    else pile1.push(ktmp, 1, 0 ); // bad
                } else
                    pile1.push(ktmp, 1, 0);     // bad for the remaining <K bases
            }

            p1 += x.Lengths(j);
            p2 += x.Lengths(j);
        }
    }

    cout << Date() << ": pile1 has " << pile1.size() << " kmers." << endl;
    ReverseSort(pile1);
    cout << Date() << ": done sorting" << endl;

    // compress counts into one record per kmer
    vec<KmerStat> pile2;
    for ( size_t i = 0; i < pile1.size(); ++i ) {
        size_t j = i+1;
        while ( j < pile1.size() && pile1[i].first == pile1[j].first ) j++;

        int goods = 0, bads = 0;
        for ( size_t l = i; l < j; ++l )  {
            bads += pile1[l].second;
            goods += pile1[l].third;
        }

        if ( bads + goods > KMER_THRESH )
            pile2.push( pile1[i].first, bads, goods );

        i = j;
    }

    // sort pile2 by bad frequency
    Sort( pile2, []( KmerStat const&a, KmerStat const& b) { 
            return ( static_cast<double>(a.second)/(a.second+a.third) ) > 
            ( static_cast<double>(b.second) / (b.second+b.third) );
            } );

    for ( size_t i = 0; i < pile2.size(); ++i )
//        cout << pile2[i].first.ToString() << " " << static_cast<double>(pile2[i].second) / (pile2[i].second+pile2[i].third) << endl;
        cout << pile2[i].first.ToString() << " " << pile2[i].second << " " << pile2[i].third << endl;

    return 0;
}
