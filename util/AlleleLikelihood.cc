/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// AlleleLikelihood.  
/// read vectors of dumbcalls and 
/// generate ranked hypotheses of allele calls and confidence ratios.

#include "Basevector.h"
#include "Bitvector.h"
#include "MainTools.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "polymorphism/DumbCall.h"
#include "solexa/SolexaPipeline.h"
#include "util/BaitMap.h"
#include "util/AlleleLikelihood.h"
#include <cassert>

double AlleleCall(const dumbcall& coverage, String allele, double epsilon)
{
    return AlleleCallWithPrior(coverage, allele, NULL, "NULL", 1, epsilon);
}

double AlleleCallWithPrior(const dumbcall& coverage, String allele, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior_table, String position, char reference_base, double epsilon)
{
    QualRecord quals;
    return AlleleCallWithPriorAndQuals(coverage, allele, prior_table, position, reference_base, epsilon, quals);
}

/// Returns the log-likelihood (without the choose(N,n) factor) of a single genotype <allele> (e.g. "AC") 
/// at the given <position> (e.g. "Chr1:123456"), given the observed base counts <coverage>, prior
/// probabilities of the alleles <prior_table>, base qualities <quals>, reference base, and probability
/// of sequencing error <epsilon>.  Error probability <epsilon> is used only if <quals> are empty, in
/// wich case it gives the probability of incorrect basecall for any given base. <position> is used to compute
/// prior and only
/// if <prior_table> is not empty; if <position>=="NULL", prior table is ignored even if it is
/// non-empty. Finally, <reference_base> is used only to compute default prior and only if prior table is empty or
/// if it contains no record for the specified position.
double AlleleCallWithPriorAndQuals(const dumbcall& coverage, String allele, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior_table, String position, char reference_base, double epsilon, QualRecord& quals)
{
    char h1 = allele.c_str()[0];  
    char h2 = allele.c_str()[1];  
    
    //double epsilon = 0.03; // 3*P(error); P(error) = Q20.
    //double epsilon = 0.000001; // for comparison to previous result.

    // total number of all observations (i.e. 35 A's, 28 C's, 3 G's and 0 T's -> 35+28+3) 
    int total_coverage = coverage.A() 
                       + coverage.C()
                       + coverage.G()
                       + coverage.T();

    int c_h1 = 0; // observed count of the first base of the genotype
    int c_h2 = 0; // observed count of the second base of the genotype

    // retrieve counts for the first and second base of the genotype (i.e. genotype "AG" -> 
    // retrieve how many A's and how many G's were seen)
    switch(h1)
    {
        case 'A' : c_h1 = coverage.A(); break;
        case 'C' : c_h1 = coverage.C(); break;
        case 'G' : c_h1 = coverage.G(); break;
        case 'T' : c_h1 = coverage.T(); break;
        default  : assert(0);
    }

    switch(h2)
    {
        case 'A' : c_h2 = coverage.A(); break;
        case 'C' : c_h2 = coverage.C(); break;
        case 'G' : c_h2 = coverage.G(); break;
        case 'T' : c_h2 = coverage.T(); break;
        default  : assert(0);
    }

    // Compute the prior.
    double p_prior;
    if (position != "NULL") { 
      p_prior = log10(Prior(position, allele, prior_table, reference_base)); 
    }
    else { p_prior = 0.0; }

    //printf("%s log10(P(%s)) = %f\n", position.c_str(), allele.c_str(), p_prior); fflush(stdout);
    
    // Compute the likelihood.
    double p_hyp;
    if (h1 == h2)
    {
        // Homozygous hypothesis.
        
        if (quals.size() == 0)
        {
            p_hyp = log10(1-epsilon)*c_h1 + log10(epsilon)*(total_coverage - c_h1);
        }
        else
        {
            p_hyp = 0.0;
            //for (unsigned int i = 0; i < quals.As.size(); i++) 
            for (int i = 0; i < coverage.A(); i++) 
            { 
                if (h1 == 'A' || h2 == 'A') { p_hyp += log10(1.0-pow(10,((float)quals.As[i]/-10.0))); }
                else                        { p_hyp += log10(pow(10,((float)quals.As[i]/-10.0))); }
            }
            //for (unsigned int i = 0; i < quals.Cs.size(); i++) 
            for (int i = 0; i < coverage.C(); i++)
            { 
                if (h1 == 'C' || h2 == 'C') { p_hyp += log10(1.0-pow(10,((float)quals.Cs[i]/-10.0))); }
                else                        { p_hyp += log10(pow(10,((float)quals.Cs[i]/-10.0))); }
            }
            //for (unsigned int i = 0; i < quals.Gs.size(); i++) 
            for (int i = 0; i < coverage.G(); i++)
            { 
                if (h1 == 'G' || h2 == 'G') { p_hyp += log10(1.0-pow(10,((float)quals.Gs[i]/-10.0))); }
                else                        { p_hyp += log10(pow(10,((float)quals.Gs[i]/-10.0))); }
            }
            //for (unsigned int i = 0; i < quals.Ts.size(); i++) 
            for (int i = 0; i < coverage.T(); i++)
            { 
                if (h1 == 'T' || h2 == 'T') { p_hyp += log10(1.0-pow(10,((float)quals.Ts[i]/-10.0))); }
                else                        { p_hyp += log10(pow(10,((float)quals.Ts[i]/-10.0))); }
            }
        }
    }
    else
    {
        // Heterozygous hypothesis.
        if (quals.size() == 0)
        {
	    p_hyp = log10(0.5-(epsilon/2))*(c_h1 + c_h2) +
                    log10(epsilon)*(total_coverage-c_h1-c_h2);
        }
        else
        {
            p_hyp = 0.0;
            //for (unsigned int i = 0; i < quals.As.size(); i++) 
            for (int i = 0; i < coverage.A(); i++)
            { 
                if (h1 == 'A' || h2 == 'A') { p_hyp += log10(0.5-pow(10,((float)quals.As[i]/-10.0))/2.0); }
                else                        { p_hyp += log10(pow(10,((float)quals.As[i]/-10.0)));         }
            }
            //for (unsigned int i = 0; i < quals.Cs.size(); i++) 
            for (int i = 0; i < coverage.C(); i++)
            { 
                if (h1 == 'C' || h2 == 'C') { p_hyp += log10(0.5-pow(10,((float)quals.Cs[i]/-10.0))/2.0); }
                else                        { p_hyp += log10(pow(10,((float)quals.Cs[i]/-10.0)));         }
            }
            //for (unsigned int i = 0; i < quals.Gs.size(); i++) 
            for (int i = 0; i < coverage.G(); i++)
            { 
                if (h1 == 'G' || h2 == 'G') { p_hyp += log10(0.5-pow(10,((float)quals.Gs[i]/-10.0))/2.0); }
                else                        { p_hyp += log10(pow(10,((float)quals.Gs[i]/-10.0)));         }
            }
            //for (unsigned int i = 0; i < quals.Ts.size(); i++) 
            for (int i = 0; i < coverage.T(); i++)
            { 
                if (h1 == 'T' || h2 == 'T') { p_hyp += log10(0.5-pow(10,((float)quals.Ts[i]/-10.0))/2.0); }
                else                        { p_hyp += log10(pow(10,((float)quals.Ts[i]/-10.0)));         }
            }
        }

    }

    // Compute the posterior.
    double p_posterior = p_hyp + p_prior;

    //printf("** %f %f %f %c %c %d %d %d\n", p_hyp, p_prior, p_posterior, h1, h2, c_h1, c_h2, total_coverage); fflush(stdout);

    return p_posterior;
}

/// Takes the counts of observed bases (how many A's, C's, ect) <coverage> and fills <alleles> and <likelihoods>
/// with all possible genotypes ("AA", "AC", "AG",...) and their corresponding log-likelihoods, given the 
/// <coverage> data. Upon return, <alleles> and <likelihoods> are synchronized and sorted in descending order
/// (so that the highest-likelihood allele goes first). <epsilon> is the probability of an incorrect basecall
/// for any given base counted in <coverage>. This method does not use priors or individual base qualities.

int ComputeAlleleCalls(const dumbcall& coverage, vec<double>& likelihoods, vec<String>& alleles, double epsilon)
{
    alleles.resize(10);
    likelihoods.resize(10);

    vec<String> bases(4);
    bases[0] = "A";
    bases[1] = "C";
    bases[2] = "G";
    bases[3] = "T";

    int k = 0;
    for (int i = 0; i < bases.isize(); i++)
    {
        for (int j = i; j < bases.isize(); j++)
        {
            alleles[k] = bases[i] + bases[j];
            k++;
        }
    }
    // alleles now filled with all possible genotypes: "AA", "AC", "AG", "AT", "CC", "CG", ...

    for (int k = 0; k < alleles.isize(); k++)
    {
        likelihoods[k] = AlleleCall(coverage, alleles[k], epsilon);
    }

    // likelihoods is now filled with log-likelihoods for each genotype from <alleles>

    // sort log-likelihoods and alleles in parallel, so that the largest likelihood (and its genotype) go first...
    ReverseSortSync(likelihoods, alleles);  
    return 1;
}

int ComputeAlleleCallsWithPrior(const dumbcall& coverage, vec<double>& likelihoods, vec<String>& alleles, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior, String position, char reference_base, double epsilon)
{
    alleles.resize(10);
    likelihoods.resize(10);

    vec<String> bases(4);
    bases[0] = "A";
    bases[1] = "C";
    bases[2] = "G";
    bases[3] = "T";

    int k = 0;
    for (int i = 0; i < bases.isize(); i++)
    {
        for (int j = i; j < bases.isize(); j++)
        {
            alleles[k] = bases[i] + bases[j];
            k++;
        }
    }

    for (int k = 0; k < alleles.isize(); k++)
    {
        likelihoods[k] = AlleleCallWithPrior(coverage, alleles[k], prior, position, reference_base, epsilon);
    }

    ReverseSortSync(likelihoods, alleles);  
    return 1;
}

int ComputeAlleleCallsWithPriorAndQuals(const dumbcall& coverage, vec<double>& likelihoods, vec<String>& alleles, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior, String position, char reference_base, double epsilon, QualRecord& quals)
{
    alleles.resize(10);
    likelihoods.resize(10);

    vec<String> bases(4);
    bases[0] = "A";
    bases[1] = "C";
    bases[2] = "G";
    bases[3] = "T";

    int k = 0;
    for (int i = 0; i < bases.isize(); i++)
    {
        for (int j = i; j < bases.isize(); j++)
        {
            alleles[k] = bases[i] + bases[j];
            k++;
        }
    }

    for (int k = 0; k < alleles.isize(); k++)
    {
        likelihoods[k] = AlleleCallWithPriorAndQuals(coverage, alleles[k], prior, position, reference_base, epsilon, quals);
    }

    ReverseSortSync(likelihoods, alleles);  
    return 1;
}

void LoadPriorTable(String file_name, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior)
{
    std::ifstream input(file_name.c_str());
    ForceAssert(input);

    float epsilon = 0.01; // ??

    String chromosome;
    int offset;
    String name;
    String strand;

    String geno1;
    float p_1; 
    int   c_1; 

    String geno2;
    float p_2; 
    int   c_2; 

    String geno3;
    float p_3; 
    int   c_3; 

    int total;

    float counts[128];
    float probs[128];

    prior = new __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>;

    while( input >> chromosome >> offset >> name >> strand
                  >> geno1 >> p_1 >> c_1
                  >> geno2 >> p_2 >> c_2
                  >> geno3 >> p_3 >> c_3 >> total )
    {
        ForceAssert(geno1.size()==3 && geno1[1]=='/');
        ForceAssert(geno2.size()==3 && geno2[1]=='/');
        ForceAssert(geno3.size()==3 && geno3[1]=='/');

        char bases[] = {'A', 'C', 'G', 'T'};
        for (int i = 0; i < 4; i++)
        {   
            counts[(int)bases[i]] = 0;
            probs[(int)bases[i]]  = 0;
        }

        counts[(int)geno1[0]] += c_1;
        counts[(int)geno1[2]] += c_1;
        counts[(int)geno2[0]] += c_2;
        counts[(int)geno2[2]] += c_2;
        counts[(int)geno3[0]] += c_3;
        counts[(int)geno3[2]] += c_3;

        probs[(int)geno1[0]] = ((float)counts[(int)geno1[0]] / (float)total) - epsilon;
        probs[(int)geno1[2]] = ((float)counts[(int)geno1[2]] / (float)total) - epsilon;
        probs[(int)geno2[0]] = ((float)counts[(int)geno2[0]] / (float)total) - epsilon;
        probs[(int)geno2[2]] = ((float)counts[(int)geno2[2]] / (float)total) - epsilon;
        probs[(int)geno3[0]] = ((float)counts[(int)geno3[0]] / (float)total) - epsilon;
        probs[(int)geno3[2]] = ((float)counts[(int)geno3[2]] / (float)total) - epsilon;

        int num_bases = 0;
        for (int i = 0; i < 4; i++)
        {   
            if ((int)counts[(int)bases[i]] != 0)
            {
                num_bases += 1;
            }
        }
        double p_n = (4.0 - num_bases)*epsilon;
        for (int i = 0; i < 4; i++)
        {
            if (counts[(int)bases[i]] == 0)
            {
                probs[(int)bases[i]] = num_bases * epsilon;
            }
        }

        p_genotype p;
        for (int i = 0; i < 4; i++)
        {
            for (int j = i; j < 4; j++)
            {
                char g[3];
                sprintf(g, "%c%c", bases[i], bases[j]);

                float p_g = probs[(int)g[0]] * probs[(int)g[1]];
                if (i != j) 
                {
                    p_g *= 2;
                }
    
                String g_str(g);
                p.genotypes.push_back(g_str);
                p.probabilities.push_back(p_g);
            }
        }
        p.count = (float)total;
        String position(chromosome + ':' + ToString(offset-1));
        (*prior)[position] = p;
//        printf("PRIOR: (%s) ", buf); for (int i = 0; i < p.probabilities.size(); i++) { printf("%s:%f ", p.genotypes[i].c_str(), p.probabilities[i]); } printf("\n"); fflush(stdout);
    }
    ForceAssert(input.eof());
}


double Prior(String genotype, char reference_base, bool using_other_prior)
{
    //printf(">> %d\n", prior.size()); fflush(stdout);

        // Fall back to the null-record.
        char h1 = genotype[0];   
        char h2 = genotype[1];   

        if (!using_other_prior)
        {
	        if      ((h1 == reference_base) && (h2 == reference_base)) { return 0.999; }
	        else if ((h1 != reference_base) && (h2 != reference_base)) { return  1e-5; }
	        else                                                       { return  1e-3; }
        }
        else
        {
	        if      ((h1 == reference_base) && (h2 == reference_base)) { return 0.9999; }
	        else if ((h1 != reference_base) && (h2 != reference_base)) { return   1e-5; }
	        else                                                       { return   1e-4; }
        }

    assert(0);
    return -1;
}


double Prior(int genotype, char reference_base, bool using_other_prior)
{
        char genotype_lookup[10][3] = {"AA",
                                       "AC",
                                       "AG",
                                       "AT",
                                       "CC",
                                       "CG",
                                       "CT",
                                       "GG",
                                       "GT",
                                       "TT"};


    //printf(">> %d\n", prior.size()); fflush(stdout);

        // Fall back to the null-record.
        char h1 = genotype_lookup[genotype][0];   
        char h2 = genotype_lookup[genotype][1];   

        if (!using_other_prior)
        {
	        if      ((h1 == reference_base) && (h2 == reference_base)) { return 0.999; }
	        else if ((h1 != reference_base) && (h2 != reference_base)) { return  1e-5; }
	        else                                                       { return  1e-3; }
        }
        else
        {
	        if      ((h1 == reference_base) && (h2 == reference_base)) { return 0.9999; }
	        else if ((h1 != reference_base) && (h2 != reference_base)) { return   1e-5; }
	        else                                                       { return   1e-4; }
        }

    assert(0);
    return -1;
}


double Prior(String position, String genotype, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior, char reference_base)
{
    //printf(">> %d\n", prior.size()); fflush(stdout);

    if ((prior != NULL) && (prior->count(position) != 0))
    {
        // We have a record for this position.
        p_genotype p = (*prior)[position];            
        for (unsigned int i = 0; i < p.genotypes.size(); i++)
        {
            if (genotype == p.genotypes[i]) { return p.probabilities[i]; }
        }
    }
    else
    {
        return Prior(genotype, reference_base);
#if 0    
        // Fall back to the null-record.
        char h1 = genotype[0];   
        char h2 = genotype[1];   

        if ((prior == NULL) || (prior->size() == 0))
        {
	        if      ((h1 == reference_base) && (h2 == reference_base)) { return 0.999; }
	        else if ((h1 != reference_base) && (h2 != reference_base)) { return  1e-5; }
	        else                                                       { return  1e-3; }
        }
        else
        {
	        if      ((h1 == reference_base) && (h2 == reference_base)) { return 0.9999; }
	        else if ((h1 != reference_base) && (h2 != reference_base)) { return   1e-5; }
	        else                                                       { return   1e-4; }
        }
#endif
    }

    assert(0);
    return -1;
}


