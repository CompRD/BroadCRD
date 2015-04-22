///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * Zinger.cc - a graph singer
 *
 *  Created on: Nov 26, 2013
 *      Author: neilw
 */
#include "Basevector.h"
#include "FastaFileset.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/HyperBasevector.h"
#include "paths/UnibaseUtils.h"
#include <vector>




class SingeEdges {
public:
    explicit SingeEdges( HyperBasevector& hbv  ) : _hbv(hbv) {};

    HyperBasevector const& hbv() const { return _hbv; }

    void Singe(unsigned max_len, unsigned max_comp_len, bool renumber,
            bool remove_rc, bool debug = false) {
        vec<vec< int > > comps;
        cout << "number of edges " << _hbv.EdgeObjectCount() << endl;
        cout << "computing components..." << endl;
        _hbv.ComponentsE( comps );
        cout << "number of components " << comps.size() << endl;
        vec<int> to_left, to_right;
        _hbv.ToLeft(to_left);
        _hbv.ToRight(to_right);
        vec<Bool> to_delete_p( _hbv.Edges().size(), false );
        for ( size_t i = 0; i < comps.size(); ++i ) {

            if ( debug ) cout << i << ":" << endl;

            bool done;
            do {
                done = true;
                for ( size_t j = 0; j < comps[i].size(); ++j ) {

                    size_t edge_no = comps[i][j];
                    BaseVec const& bv = _hbv.Edges()[ edge_no ];
                    int vl = to_left[ edge_no ];
                    int vr = to_right[ edge_no ];

                    if ( !to_delete_p[edge_no] && bv.size() <= max_len ) {

                        bool saved = false;
                        for ( auto const& f : _hbv.FromEdgeObj( vr ) ) {
                            if ( !to_delete_p[f] ) {
                                saved = true;
                                break;
                            }
                        }

                        if ( saved ) {
                            saved = false;
                            for ( auto const& f : _hbv.ToEdgeObj( vl ) ) {
                                if ( !to_delete_p[f] ) {
                                    saved = true;
                                    break;
                                }
                            }
                        }

                        if ( !saved ) {
                            to_delete_p[edge_no] = true;
                            done = false;
                        }

                    }
                }
            } while ( !done );

            if ( debug ) {
                for ( size_t j = 0; j < comps[i].size(); ++j ) {
                    size_t edge_no = comps[i][j];
                    size_t bv_size = _hbv.Edges()[ edge_no ].size();
                    cout << "\t" << edge_no << "(" << bv_size << ")";
                    if ( to_delete_p[edge_no] ) cout << "*";
                    cout << endl;
                }
            }

        }

        vec<int> to_delete;
        for ( size_t i = 0; i < to_delete_p.size(); ++i )
            if ( to_delete_p[i] ) to_delete.push_back(i);

        _hbv.DeleteEdges( to_delete );
        if ( renumber ) _hbv.RemoveDeadEdgeObjects();
        _hbv.RemoveEdgelessVertices();

        _hbv.TestValid();

        if (debug) {
            Ofstream( dot_out, "beforemerge.dot" );
            _hbv.PrintSummaryDOT0w( dot_out, true, true, true );
            _hbv.DumpFasta( "beforemerge.fasta", false );
        }

        // merge contiguous edges
        _hbv.ToLeft(to_left);
        _hbv.ToRight(to_right);
        for ( int vi = 0; vi < _hbv.N() ; vi++ ) {
            if ( _hbv.FromSize(vi) == 1 &&
                    _hbv.ToSize(vi) == 1 &&
                    _hbv.From(vi)[0] != vi      // excl loops
                        ) {
                // in-degree = out-degree = 1
                // so merge
                if (debug) cout << "MERGING: " << endl;
                int edge_l = _hbv.EdgeObjectIndexByIndexTo( vi, 0 );
                int edge_r = _hbv.EdgeObjectIndexByIndexFrom( vi, 0 );
                int vert_l = _hbv.To(vi)[0];
                int vert_r = _hbv.From(vi)[0];
                if (debug) PRINT5( vi, edge_l, edge_r, vert_l, vert_r );
                BaseVec const& bv_l = _hbv.EdgeObject(edge_l);
                BaseVec const& bv_r = _hbv.EdgeObject(edge_r);
                _hbv.AddEdge( vert_l, vert_r, TrimCat( _hbv.K(), bv_l, bv_r ) );
                _hbv.DeleteEdgeTo( vi, 0 );
                _hbv.DeleteEdgeFrom( vi, 0 );
            }
        }

        _hbv.RemoveEdgelessVertices();

        cout << "computing components again..." << endl;
        _hbv.ComponentsE( comps );
        cout << "number of components " << comps.size() << endl;

        if ( max_comp_len ) {
            cout << "removing components equal or smaller than " << max_comp_len << endl;
            vec<int> to_delete;
            size_t count = 0;
            for ( auto const& comp : comps ) {
                size_t total = 0;
                std::for_each( comp.begin(), comp.end(),
                        [&total,this]( int bi ) {
                        total+=_hbv.Edges()[bi].size();
                        } );
                if ( total <= max_comp_len )  {
                    for ( auto id : comp ) to_delete.push_back(id);
                    count++;
                }
            }

            if ( to_delete.size() > 0 ) {
                UniqueSort(to_delete);
                _hbv.DeleteEdges( to_delete );
            }
            cout << "we removed " << count << " small components" << endl;
        }

        if ( renumber ) _hbv.RemoveDeadEdgeObjects();
        _hbv.RemoveEdgelessVertices();

        if ( remove_rc ) _hbv.DeleteReverseComplementComponents();
    }

private:
    HyperBasevector& _hbv;
};




#include "MainTools.h"

int main( int argc, char** argv )
{
    String empty;
    RunTime();
    BeginCommandArguments;
    CommandArgument_Bool_OrDefault_Doc( RENUMBER, True,
        "whether to remove dead edges so that the edges are renumbered contiguously" );
    CommandArgument_Bool_OrDefault_Doc( REMOVE_RC, True,
        "remove RC components" );
    CommandArgument_String_Doc(IN_HBV,
        "input HyperBasevector");
    CommandArgument_String_Doc(OUT_HEAD,
        "writes something");
    CommandArgument_UnsignedInt_OrDefault(NUM_THREADS,0u);
    CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB,0);
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_LEN,100,"largest edge to remove");
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_COMP_LEN,100,"largest component to remove");
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);
    SetMaxMemory(MAX_MEMORY_GB<<30);

    HyperBasevector hbv(IN_HBV);
    SingeEdges( hbv ).Singe( MAX_LEN, MAX_COMP_LEN, RENUMBER, REMOVE_RC );

    BinaryWriter writer( OUT_HEAD+".hbv" );
    hbv.writeBinary( writer );
    Ofstream( dot_out, OUT_HEAD+".dot" );
    hbv.PrintSummaryDOT0w( dot_out, true, false, true );
    hbv.DumpFasta( OUT_HEAD+".fasta", false );
}

