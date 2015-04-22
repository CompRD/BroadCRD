///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file DigraphTests.h
 * \author tsharpe
 * \date Apr 7, 2010
 *
 * \brief
 */
#ifndef DIGRAPHTESTS_H_
#define DIGRAPHTESTS_H_

#include <cxxtest/TestSuite.h>
#include "graph/DigraphTemplate.h"
#include "feudal/BinaryStream.h"
#include <unistd.h>

class DigraphTests : public CxxTest::TestSuite
{
public:
    void testEquals()
    {
        digraphE<double> g1;
        g1.AddVertices(4);
        g1.AddEdge(0,1,10.);
        g1.AddEdge(1,2,8.);
        g1.AddEdge(2,1,6.);
        g1.AddEdge(2,3,4.);

        digraphE<double> g2;
        g2.AddVertices(4);
        g2.AddEdge(2,3,4.);
        g2.AddEdge(1,2,8.);
        g2.AddEdge(2,1,6.);
        g2.AddEdge(0,1,10.);

        TS_ASSERT(EqualExceptEdgeObjectOrder(g1,g2));
        g1.Reverse();
        g2.Reverse();
        TS_ASSERT(EqualExceptEdgeObjectOrder(g1,g2));
    }

    void testBinaryStreaming()
    {
        digraphE<double> g2;
        g2.AddVertices(4);
        g2.AddEdge(2,3,4.);
        g2.AddEdge(1,2,8.);
        g2.AddEdge(2,1,6.);
        g2.AddEdge(0,1,10.);
        digraphE<double> g3;

        if ( true )
        {
            BinaryWriter w("tmp");
            w.write(g2);
            w.close();

            BinaryReader r("tmp");
            r.read(&g3);

            unlink("tmp");
        }

        TS_ASSERT_EQUALS(g2,g3);
        g2.Reverse();
        g3.Reverse();
        TS_ASSERT_EQUALS(g2,g3);
    }
};

#endif /* DIGRAPHTESTS_H_ */
