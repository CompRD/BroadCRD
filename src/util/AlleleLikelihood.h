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

#ifndef ALLELE_LIKELIHOOD_H
#define ALLELE_LIKELIHOOD_H

#include "Basevector.h"
#include "Bitvector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "polymorphism/DumbCall.h"
#include "solexa/SolexaPipeline.h"
#include "util/BaitMap.h"
#include <ext/hash_map>


struct StringHashFunction
{
  size_t operator()(const String& s) const
  {
    __gnu_cxx::hash<const char*> h;
    return h(s.c_str());
  }
};

struct p_genotype
{
    vec<float>  probabilities;
    vec<String> genotypes;
    float       count;
};

struct QualRecord
{
    vec<qual_t> As;
    vec<qual_t> Cs;
    vec<qual_t> Gs;
    vec<qual_t> Ts;

    int size() { return As.size() + Cs.size() + Gs.size() + Ts.size(); }
};

double AlleleCall(const dumbcall& coverage, String allele, double epsilon=0.03);
double AlleleCallWithPrior(const dumbcall& coverage, String allele, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior_table, String position, char reference_base, double epsilon=0.03);
double AlleleCallWithPriorAndQuals(const dumbcall& coverage, String allele, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior_table, String position, char reference_base, double epsilon, QualRecord& quals);
int ComputeAlleleCalls(const dumbcall& coverage, vec<double>& likelihoods, vec<String>& alleles, double epsilon=0.03);
int ComputeAlleleCallsWithPrior(const dumbcall& coverage, vec<double>& likelihoods, vec<String>& alleles, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior, String position, char reference_base, double epsilon=0.03);
int ComputeAlleleCallsWithPriorAndQuals(const dumbcall& coverage, vec<double>& likelihoods, vec<String>& alleles, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior, String position, char reference_base, double epsilon, QualRecord& quals);

// Note, new values added to prior will overwrite old ones.
void LoadPriorTable(String file_name, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior);
double Prior(String genotype, char reference_base, bool using_other_prior=false);
double Prior(int genotype, char reference_base, bool using_other_prior=false);
double Prior(String position, String genotype, __gnu_cxx::hash_map<String,p_genotype,StringHashFunction>* prior, char reference_base);


/////////////////////////////////////////
// Code for calling with quals.


#define GENO_AA 0x00
#define GENO_AC 0x01
#define GENO_AG 0x02
#define GENO_AT 0x03
#define GENO_CC 0x04
#define GENO_CG 0x05
#define GENO_CT 0x06
#define GENO_GG 0x07
#define GENO_GT 0x08
#define GENO_TT 0x09


struct genotype_likelihoods
{
    float likelihoods[10];


    genotype_likelihoods()
    {
        memset(likelihoods, 0, sizeof(likelihoods));
    }

    genotype_likelihoods(char ref, char base, unsigned char qual)
    { init(ref,base,qual); }

    // Computing P(T|G), for all genotypes G. T is the base contributed by the current alignment A.
    // P(T|G) = P(T|G,A) + P(T|G,~A)
    // P(T|G,A) = P(A)L(T|G)
    // P(T|G,~A) = P(T|~A) = (1 - P(A)) * sum_over_G( L(T|G)P(G) ) , where P(G) = 0.1

    genotype_likelihoods(char ref, char base, unsigned char qual, const vec<double>& probs) 
    { 
        init(ref, base, qual);

        double p_genotype_uniform = 0.1; 

        double p_base = 0.0;
        // sum over G( L(T|G) * P(G) )
        for ( unsigned int i = 0; i < 10; i++)
        {
	    p_base = pow(10.0, this->likelihoods[i]) * p_genotype_uniform;
        }
        // P(A)
        const double p_best_alignment  = probs[0];
        // 1 - P(A)
        // Could assert that probs sums to 0 here 
        const double p_alt_alignment = 1 - p_best_alignment; 
        // P(T|G),~A) = P(T|~A)
        const double p_base_alt = p_alt_alignment * p_base;
        // P(T|G) = P(T|G,A) + P(T|G,~A), for all G
        // Note that that P(T|G,~A) ==  P(T|~A) makes this term loop-invariant below
        for ( unsigned int i = 0; i < 10; i++)
        {
	    this->likelihoods[i] = p_best_alignment * this->likelihoods[i] + p_base_alt; 
        }
    }

    void init(char ref, char base, unsigned char qual)
    {
        dumbcall C;

        int base_int;
        switch(base)
        {
                case 'A' : base_int = 0; break;
                case 'C' : base_int = 1; break;
                case 'G' : base_int = 2; break;
                case 'T' : base_int = 3; break;
                case 'D' : base_int = 4; break;
                case 'I' : base_int = 5; break;
                default: FatalErr("Unrecognized base: " << base);
        }

        C.base[base_int] += 1;

        if (qual == 0) { qual = 1; } // map Q0 to something that won't make log(0) happen.

        float epsilon = pow(10.0, (double)qual / -10);

        //printf("EPSILON: %f %f\n", (double)qual, epsilon); fflush(stdout);

        vec<double> local_likelihoods;
        vec<String> alleles;
        ComputeAlleleCalls(C, local_likelihoods, alleles, epsilon);
        //AtomicAlleleCalls(base, local_likelihoods, alleles, epsilon);

        SortSync(alleles, local_likelihoods);
        for ( unsigned int i = 0; i < local_likelihoods.size(); i++)
        {
            this->likelihoods[i] = local_likelihoods[i];
        }

        //printf("\nFNORD %c | %s | %s\n", ref, C.TabularString(ref).c_str(), this->TabularString(ref).c_str()); fflush(stdout);
    }

    static vec<String> atomic_alleles_A;
    static vec<String> atomic_alleles_C;
    static vec<String> atomic_alleles_G;
    static vec<String> atomic_alleles_T;
    static vec< vec<double> > atomic_likelihoods;


    static void InitLookupTables()
    {  

        atomic_likelihoods.resize(10000);

	    for (int i = 0; i < 10000; i++)
	    {
            atomic_likelihoods[i].resize(10);

	        double eps = (1.0 / (double)i);
	        double eps_over_2 = eps / 2.0;
	
	        atomic_likelihoods[i][0] = log10(1-eps);
	        atomic_likelihoods[i][1] = log10(0.5-eps_over_2);
	        atomic_likelihoods[i][2] = log10(0.5-eps_over_2);
	        atomic_likelihoods[i][3] = log10(0.5-eps_over_2);
	        atomic_likelihoods[i][4] = log10(eps);
	        atomic_likelihoods[i][5] = log10(eps);
	        atomic_likelihoods[i][6] = log10(eps);
	        atomic_likelihoods[i][7] = log10(eps);
	        atomic_likelihoods[i][8] = log10(eps);
	        atomic_likelihoods[i][9] = log10(eps);
	    }

        atomic_alleles_A.resize(10);

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
                atomic_alleles_A[k] = bases[i] + bases[j];
                k++;
            }
        }

        atomic_alleles_C = atomic_alleles_A; 
        atomic_alleles_G = atomic_alleles_A; 
        atomic_alleles_T = atomic_alleles_A; 

        Sort(atomic_alleles_A, atomic_allele_sorter('A'));
        Sort(atomic_alleles_C, atomic_allele_sorter('C'));
        Sort(atomic_alleles_G, atomic_allele_sorter('G'));
        Sort(atomic_alleles_T, atomic_allele_sorter('T'));
    }

    static void TestSpeedUp()
    {
        dumbcall C; C.base[0] = 1;
        double epsilon = 0.001;
    
        vec<String> A;
        vec<double> L;

        time_t start; 
        time_t end; 

        start = clock();
        for ( int i = 0; i < 1e6; i++)
        {
            ComputeAlleleCalls(C,  L, A, epsilon);
        }
        end = clock();
        printf("Old Way: %0.04f seconds\n", (((double)end - (double)start)/(double)CLOCKS_PER_SEC)); fflush(stdout);

        start = clock();
        for ( int i = 0; i < 1e6; i++)
        {
            genotype_likelihoods g('A', 'A', 30);
        }
        end = clock();
        printf("New Way: %0.04f seconds\n", (((double)end - (double)start)/(double)CLOCKS_PER_SEC)); fflush(stdout);
    }

    class atomic_allele_sorter
    {
        public:

        char base;
        atomic_allele_sorter(char base) { this->base = base; }
        bool operator()(const String& a, const String& b) const
        {
            int c_a = 0; 
            if (a[0] == base) { c_a += 1; } 
            if (a[1] == base) { c_a += 1; } 

            int c_b = 0; 
            if (b[0] == base) { c_b += 1; } 
            if (b[1] == base) { c_b += 1; } 

            return (c_a > c_b);
        }
    };  

    void AtomicAlleleCalls(char base, vec<double>& likelihoods, vec<String>& alleles, double epsilon)
    {
        switch (base)
        {
            case 'A': alleles = genotype_likelihoods::atomic_alleles_A; break;
            case 'C': alleles = genotype_likelihoods::atomic_alleles_C; break;
            case 'G': alleles = genotype_likelihoods::atomic_alleles_G; break;
            case 'T': alleles = genotype_likelihoods::atomic_alleles_T; break;
        }

        int eps = (int)(1.0 / epsilon);
    
//        printf(">> %f %d\n", epsilon, eps);

        likelihoods = genotype_likelihoods::atomic_likelihoods[eps];

//        for (int i = 0; i < 10; i++) { printf("%s:%f ", genotype_likelihoods::atomic_alleles[i].c_str(), genotype_likelihoods::atomic_likelihoods[eps][i]); } printf("\n"); fflush(stdout);
//        for (int i = 0; i < 10; i++) { printf("%s:%f ", alleles[i].c_str(), likelihoods[i]); } printf("\n"); fflush(stdout);

//        SortSync(alleles, likelihoods);
//        for (int i = 0; i < 10; i++) { printf("%s:%f ", alleles[i].c_str(), likelihoods[i]); } printf("\n"); fflush(stdout);
    }

    genotype_likelihoods& operator+=(genotype_likelihoods& G)
    {
        for (int i = 0; i < 10; i++)
        {
            this->likelihoods[i] += G.likelihoods[i];
        }
        return *this;
    } 

    void Accumulate(genotype_likelihoods& G)
    {
        for (int i = 0; i < 10; i++)
        {
            this->likelihoods[i] += G.likelihoods[i];
        }
    } 

    genotype_likelihoods Times(float x)
    {
        genotype_likelihoods G;
        for (int i = 0; i < 10; i++)
        {
            G.likelihoods[i] = this->likelihoods[i] + x;
        }
        return G;
    }

    genotype_likelihoods Plus(float x)
    {
        genotype_likelihoods G;
        for (int i = 0; i < 10; i++)
        {
            G.likelihoods[i] = log10(pow(10,this->likelihoods[i]) + pow(10,x));
        }
        return G;
    }

    double Sum()
    {
        double ans = 0;
        for (int i = 0; i < 10; i++) { ans += pow(10,this->likelihoods[i]); }
        return log10(ans);
    }

    String TabularString(char ref) const
    {
        char buf[8192];
        memset(buf, 0x00, 8192);

        vec<String> genotypes;
        genotypes.push_back("AA");
        genotypes.push_back("AC");
        genotypes.push_back("AG");
        genotypes.push_back("AT");
        genotypes.push_back("CC");
        genotypes.push_back("CG");
        genotypes.push_back("CT");
        genotypes.push_back("GG");
        genotypes.push_back("GT");
        genotypes.push_back("TT");
        
        sprintf(buf, "%c ", ref);
        for (int i = 0; i < 10; i++)
        {
            sprintf(buf, "%s %s:%0.5f", buf, genotypes[i].c_str(), likelihoods[i]);
        }

        return String(buf);
    }

    friend ostream& operator<<(ostream& os,const genotype_likelihoods& g) 
    {
        os << g.TabularString('N');
        return os;
    }

    void GetSorted(vec<String>& alleles, vec<double>& likelihoods_out)
    {
#if 1
        vec<int> allele_tokens(10);
        allele_tokens[0] = GENO_AA;
        allele_tokens[1] = GENO_AC;
        allele_tokens[2] = GENO_AG;
        allele_tokens[3] = GENO_AT;
        allele_tokens[4] = GENO_CC;
        allele_tokens[5] = GENO_CG;
        allele_tokens[6] = GENO_CT;
        allele_tokens[7] = GENO_GG;
        allele_tokens[8] = GENO_GT;
        allele_tokens[9] = GENO_TT;


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

#else
        alleles.resize(10);
        alleles[0] = "AA";
        alleles[1] = "AC";
        alleles[2] = "AG";
        alleles[3] = "AT";
        alleles[4] = "CC";
        alleles[5] = "CG";
        alleles[6] = "CT";
        alleles[7] = "GG";
        alleles[8] = "GT";
        alleles[9] = "TT";
#endif


        likelihoods_out.resize(10);
        for (int i = 0; i < 10; i++)
        {
            likelihoods_out[i] = likelihoods[i];
        }

        //ReverseSortSync(likelihoods_out, alleles);
        ReverseSortSync(likelihoods_out, allele_tokens);

#if 1
        alleles.resize(10);
        for (int i = 0; i < 10; i++)
        {
            alleles[i] = genotype_lookup[allele_tokens[i]];
        }
#endif
    }

    void GetSorted(vec<int>& allele_tokens, vec<double>& likelihoods_out)
    {
        allele_tokens.resize(10);
        allele_tokens[0] = GENO_AA;
        allele_tokens[1] = GENO_AC;
        allele_tokens[2] = GENO_AG;
        allele_tokens[3] = GENO_AT;
        allele_tokens[4] = GENO_CC;
        allele_tokens[5] = GENO_CG;
        allele_tokens[6] = GENO_CT;
        allele_tokens[7] = GENO_GG;
        allele_tokens[8] = GENO_GT;
        allele_tokens[9] = GENO_TT;


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


        likelihoods_out.resize(10);
        for (int i = 0; i < 10; i++)
        {
            likelihoods_out[i] = likelihoods[i];
        }

        //ReverseSortSync(likelihoods_out, alleles);
        ReverseSortSync(likelihoods_out, allele_tokens);
    }

    void ApplyPrior(char r)
    {
#if 0
        likelihoods[0] += log10(Prior("AA", r, false)); 
        likelihoods[1] += log10(Prior("AC", r, false)); 
        likelihoods[2] += log10(Prior("AG", r, false)); 
        likelihoods[3] += log10(Prior("AT", r, false)); 
        likelihoods[4] += log10(Prior("CC", r, false)); 
        likelihoods[5] += log10(Prior("CG", r, false)); 
        likelihoods[6] += log10(Prior("CT", r, false)); 
        likelihoods[7] += log10(Prior("GG", r, false)); 
        likelihoods[8] += log10(Prior("GT", r, false)); 
        likelihoods[9] += log10(Prior("TT", r, false)); 
#else
        likelihoods[0] += log10(Prior(GENO_AA, r, false)); 
        likelihoods[1] += log10(Prior(GENO_AC, r, false)); 
        likelihoods[2] += log10(Prior(GENO_AG, r, false)); 
        likelihoods[3] += log10(Prior(GENO_AT, r, false)); 
        likelihoods[4] += log10(Prior(GENO_CC, r, false)); 
        likelihoods[5] += log10(Prior(GENO_CG, r, false)); 
        likelihoods[6] += log10(Prior(GENO_CT, r, false)); 
        likelihoods[7] += log10(Prior(GENO_GG, r, false)); 
        likelihoods[8] += log10(Prior(GENO_GT, r, false)); 
        likelihoods[9] += log10(Prior(GENO_TT, r, false)); 
#endif
    }

    void ApplyPrior(String genotype, double penalty)
    {
        double match = log10(1.0 - penalty);
        double mismatch = log10(penalty);

        likelihoods[0] += (genotype == "AA" ? match : mismatch);
        likelihoods[1] += (genotype == "AC" ? match : mismatch);
        likelihoods[2] += (genotype == "AG" ? match : mismatch);
        likelihoods[3] += (genotype == "AT" ? match : mismatch);
        likelihoods[4] += (genotype == "CC" ? match : mismatch);
        likelihoods[5] += (genotype == "CG" ? match : mismatch);
        likelihoods[6] += (genotype == "CT" ? match : mismatch);
        likelihoods[7] += (genotype == "GG" ? match : mismatch);
        likelihoods[8] += (genotype == "GT" ? match : mismatch);
        likelihoods[9] += (genotype == "TT" ? match : mismatch);
    }

    double GetLikelihood(String genotype)
    {
        if (genotype == "AA") { return likelihoods[0]; }
        if (genotype == "AC") { return likelihoods[1]; }
        if (genotype == "AG") { return likelihoods[2]; }
        if (genotype == "AT") { return likelihoods[3]; }
        if (genotype == "CC") { return likelihoods[4]; }
        if (genotype == "CG") { return likelihoods[5]; }
        if (genotype == "CT") { return likelihoods[6]; }
        if (genotype == "GG") { return likelihoods[7]; }
        if (genotype == "GT") { return likelihoods[8]; }
        if (genotype == "TT") { return likelihoods[9]; }
   
        ForceAssert(false);
        return 0; 
    }

    double GetLikelihood(int genotype)
    {
        return likelihoods[genotype];
    }

    void Print(ostream& out)
    {
        vec<String> alleles;
        vec<double> likelihoods;
        this->GetSorted(alleles, likelihoods);
        SortSync(alleles, likelihoods);
        vec< vec<String> > rows(2); 
        rows[0] = alleles;
        for (int i = 0; i < 10; i++) { rows[1].push_back(ToString(likelihoods[i])); }
        PrintTabular(out, rows, 2, "r");
    }

};



#endif





















