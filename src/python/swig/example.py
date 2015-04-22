#!/usr/bin/env python
import scaffold

g = scaffold.LinkGraph()

scaffold.read_link_graph("/wga/dev/WGAdata/projects/ALLPATHS/bacteria/Escherichia_coli_29193/3031H/ecoli_dexter_9227_05_06_10/ASSEMBLIES/test/initial_scaffolds.graph", g)

supers = scaffold.superb_vector()

scaffold.ReadSuperbs("/wga/dev/WGAdata/projects/ALLPATHS/bacteria/Escherichia_coli_29193/3031H/ecoli_dexter_9227_05_06_10/ASSEMBLIES/test/initial_scaffolds.superb", supers)

def print_superbs():
    print len(supers)
    print g.N()
    for super in supers:
        print super.TrueLength()

def print_graph():
    for vertex in range(g.N()):
        to_verticies =  g.To(vertex)
        to_edges = g.ToEdgeObj(vertex)
        from_list =  g.From(vertex)
        from_edges = g.FromEdgeObj(vertex)
        nodes = str.join(', ', map(str, to_verticies))
        print "  To %d: %s" % (vertex, nodes)
        # get edge objects
        for edge in to_edges:
            link = g.EdgeObject(edge)
            print link.Sep(), link.Dev(), link.Weight(), link.Score()
            print link.Win1()
        nodes = str.join(', ', map(str, from_list))
        print "From %d: %s" % (vertex, nodes)

print_graph()
print_superbs()
