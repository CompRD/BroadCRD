%module scaffold 

%{
#include <fcntl.h>

#include "Superb.h"
#include "PairsManager.h"
#include "paths/Alignlet.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/reporting/CLinkBundle.h"
#include "paths/BuildScaffoldGraph.h"
#include "paths/Sepdev.h"

#include "Vec.h"
#include "String.h"
%}

%include "std_vector.i"
%include "std_pair.i"
%include "std_string.i"

/* superb */
class superb 
{
public:
    superb();
    int Ntigs() const;
    void SetNtigs(int n);
    int TrueLength() const;
};

/* CLinkBundle */
class CLinkBundle
{
public:
    CLinkBundle();
    int Sep() const;
    int Dev() const;
    int Weight() const;
    float Score() const;
    std::pair<double,double> Win1() const;
    std::pair<double,double> Win2() const;
};

/* digraph */
class digraph
{
public:
    digraph();
    explicit digraph(int n);
};

/* vec */
template<class T> class vec: public std::vector<T>
{
public:
    vec();
    size_t size();
};

/* digraphE */
template<class E> class digraphE: public digraph
{
public:
    digraphE();
    int N();
    const vec<int> &To(int v);
    const vec<int> &From(int v);
    const vec<int> &FromEdgeObj(int v);
    const vec<int> &ToEdgeObj(int v);
    const E& EdgeObject(int i) const;
};


void BuildScaffold(const std::string &root, digraphE<CLinkBundle> &graph, vec<superb> &supers);

/* file IO */
void WriteLinkGraph(char *fn, const digraphE<CLinkBundle> &g);
void ReadLinkGraph(char *fn, digraphE<CLinkBundle> &g);
void ReadSuperbs(const std::string &fn, vec<superb> &supers);
void WriteSuperbs(const std::string &fn, const vec<superb> &supers);

%{

/* wrapper functions */
void WriteLinkGraph(char *fn, const digraphE<CLinkBundle> &g)
{
    int fd = creat(fn, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);

    BinaryWrite(fd, g);
}

void ReadLinkGraph(char *fn, digraphE<CLinkBundle> &g)
{
    int fd = open(fn,  O_RDONLY);
    BinaryRead(fd, g);
}

void BuildScaffold(const std::string &root, digraphE<CLinkBundle> &graph, vec<superb> &supers)
{
    String pairs_file = root + "/../../" + "scaffold_reads.pairs";
    String aligns_file = root + "/" + "scaffold_reads.qltoutlet";
    String index_file =  root + "/" + "scaffold_reads.qltoutlet.index";
    String supers_file =  root + "/" + "initial_scaffolds.superb";
    vec<alignlet> aligns;
    vec<int> index;
    digraphE<sepdev> unused;

    PairsManager pairs(pairs_file);
    ReadSuperbs(supers_file, supers);
    BinaryRead3(aligns_file, aligns);
    BinaryRead3(index_file, index);

    BuildScaffoldGraph(pairs, supers, aligns, index, unused, &graph);
}


%}

/* template instantiations */
%template(std_vector_int) std::vector<int>;
%template(std_vector_superb) std::vector<superb>;
%template(std_pair_double_double) std::pair<double, double>;
%template(int_vector) vec<int>;
%template(superb_vector) vec<superb>;
%template(LinkGraph) digraphE<CLinkBundle>;
