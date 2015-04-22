#ifndef __BASEVEC_TESTS_H
#define __BASEVEC_TESTS_H

#include <cxxtest/TestSuite.h>
#include "dna/Bases.h"
#include <iostream>
#include <fstream>
#include "feudal/BaseVec.h"
#include "Fastavector.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "paths/Alignlet.h"
#include <stdlib.h>
#include <time.h>
#include "paths/ContigsManager.h"
#include "paths/Alignlet.h"

// build a set of related fasta vectors
// align the vectors, build alignment class

class ContigsTests : public CxxTest::TestSuite
{
public:
    void setUp()
    { 
        srandom(time(NULL));
    }

    void tearDown()
    { }

    void randomize_bases(fastavector::iterator pos, fastavector::iterator end)
    {
        while(pos != end)
        {
            const GeneralizedBase &gb = GeneralizedBase::fromChar('N');
            *pos = as_base(gb.random());
            ++pos;
        }

    }

    void generate_fasta(fastavector &fv)
    {
        randomize_bases(fv.begin(), fv.end());
    }

    void generate_alignable_fasta(fastavector &source, fastavector &target, 
                                    size_t overlap) 
    {
        fastavector::iterator target_itr = target.begin();
        fastavector::iterator source_itr = source.begin(); 

        randomize_bases(target_itr, target_itr + overlap);
        target_itr += overlap;

        for( ; target_itr != target.end(); ++target_itr)
        {
            *target_itr = *source_itr;
            ++source_itr;
        }
    }

    void make_reads_from_fasta(vec<fastavector> &contigs, size_t contig_id,
                                vec<fastavector> &reads,
                                vec<alignlet> &aligns,
                                vec<int> &aligns_index, 
                                size_t readlen, size_t stepsize) 
    {
        fastavector &contig = contigs[contig_id];
        size_t steps = (contig.size() - readlen)  / stepsize;
        size_t offset = aligns.size();

        aligns_index.resize(steps + offset);
        aligns.resize(steps + offset);
        reads.resize(steps + offset);

        for(size_t step = 0; step < steps; ++step)
        {
            size_t read_id = step + offset;
            
            // make the read
            reads[read_id] = fastavector(contig.begin() + (step * stepsize),
                            contig.begin() + (step * stepsize) + readlen);
            // make the alignlet
            alignlet al = alignlet(step * stepsize,             
                                    step * stepsize + readlen, 
                                    contig_id, 
                                    readlen, true);
            // put the alignment into the aligns list
            aligns[read_id] = al;
            // relate the alignment with the associated read_id
            aligns_index[read_id] = read_id;
        }
    }


    ContigsManager prepare_dataset(vec<fastavector> &contigs,
                                    vec<fastavector> &reads,
                                    vec<alignlet> &aligns,
                                    vec<int> &aligns_index)
    {
        const size_t contig_len = 1000;
        size_t contig_overlap = 200;
        size_t contig_count = 100;
        size_t read_step = 10;
        size_t read_len = 100;

        contigs.resize(contig_count);
        fastavector fv;
        fv.resize(contig_len);
        generate_fasta(fv);
        contigs[0] = fv;

        size_t contig_id;
        for(contig_id = 0; contig_id < (contig_count - 1); ++contig_id)
        {
            make_reads_from_fasta(contigs, contig_id, reads, aligns, 
                                    aligns_index, read_len, read_step);
            fastavector &fv1 = contigs[contig_id];
            fastavector fv2;
            fv2.resize(contig_len);
            generate_alignable_fasta(fv1, fv2, contig_overlap);
            contigs[contig_id + 1] = fv2;
        }
        make_reads_from_fasta(contigs, contig_id, reads, aligns, 
                                aligns_index, read_len, read_step);
	
	vec<alignlet> dummy(0);
	vec< vec<int> > dummy_index(0);
	return ContigsManager(contigs, aligns, aligns_index, dummy, dummy_index );
    }

    void print_fasta(fastavector &fv)
    {
        printf("\n\n%p:\n", static_cast<void*>(&fv));
        for(fastavector::iterator itr = fv.begin(); itr != fv.end(); ++itr)
        {
            printf("%c", *itr);
        }
        printf("\n");
    }

    void verify_fasta(fastavector &key, fastavector &query)
    {
        TS_ASSERT_EQUALS(key.size(), query.size());
        for(size_t idx = 0; idx < key.size(); ++idx)
        {
            TS_ASSERT_EQUALS(key[idx], query[idx]);
        }

    }

    void test_creation()
    {
        vec<fastavector> contigs;
        vec<fastavector> reads;
        vec<alignlet> aligns;
        vec<int> aligns_index;
        ContigsManager cm = prepare_dataset(contigs, reads, aligns, 
                                                        aligns_index);
        vec<fastavector> fv;
        vec<size_t> res = cm.SplitContig(0, 200);
        
        print_fasta(contigs[res[0]]);
        print_fasta(contigs[res[1]]);
    }
};

#endif
