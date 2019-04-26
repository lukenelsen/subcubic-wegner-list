# To-Do:
    # Revise graph_power and core_square_graph, copy to wegner.py
    # After paper is finished, rename lemma_9 appropriately.


#------------------------------------------------------
# Contents
#------------------------------------------------------

#    core_square_graph
#    length_equals_8
#    length_equals_9
#    (main command)



#------------------------------------------------------
# Libraries, Comments, etc.
#------------------------------------------------------

#import time
## For displaying runtime in run_lemma_28.

#def timestring(time):
    ## Displays time length:  ?m ?.?s
    
    #t = int(time)
    #d = str(int((10*(time-t))%10))
    #s = str(t%60)
    #m = int(t/60)
    #string = str(m)+"m "+s+"."+d+"s"
    #return string



#------------------------------------------------------
# The Routines
#------------------------------------------------------

import choosability as ch
# We use fChoosable as a black box.



import copy
# We use deepcopy in graph_power.



#Takes a graph G and integer power as input, returns a revised graph with edges between pairs of vertices whose original distance was at most the power.
def graph_power(G,power):
    if power == 0:
        return Graph(matrix.zero(G.order()))
    elif power == 1:
        return G
    elif power >= 2:
        E = copy.deepcopy(G.edges())
        M = G.distance_matrix()
        for j in range(1,G.order()):
            for i in range(j):
                if 2 <= M[i][j] <= power:
                    E.append((i,j,None))
        return Graph(E)
    else:
        print 'Error:  Given power is invalid.'







#identifications is a list of lists:  
#   [stem-to-root 2-identifications, stem-to-stem 2-identifications, stem-to-stem 3-identifications].
#
def core_square_graph(underlying_graph,identifications):
    #Start with f
    f = []
    ids1 = {x for y in identifications[0] for x in y}
    ids2 = {x for y in identifications[1] for x in y}
    ids3 = {x for y in identifications[2] for x in y}
    #2-vertices in the graph which are stem-to-root 2-identified become 3-vertices.
    twos = {x for x in underlying_graph.vertices() if underlying_graph.degree(x)==2 if not x in ids1}
    for v in underlying_graph.vertices():
        #two_nbrs is closed neighborhood of v in underlying_graph.
        two_nbrs = set([])
        if v in twos:
            two_nbrs |= {v,}
        two_nbrs |= {x for x in underlying_graph.neighbors(v) if x in twos}
        #We count how many stem-to-stem identifications the elements of two_nbrs are part of.
        #Each contributes 1 to v's restrictions.
        count = 0
        for x in identifications[1]:
            if two_nbrs & set(x):
                count += 1
                two_nbrs -= set(x)
        for x in identifications[2]:
            if two_nbrs & set(x):
                count += 1
                two_nbrs -= set(x)
        #Anything remaining in two-nbrs is a 2-vtx without identification.  These too contribute 1 restriction.
        count += len(two_nbrs)
        #If v is a 3-vtx or part of a stem-to-root 2-ident. or a stem-to-stem 3-ident., that's all the restrictions.
        #But if v is part of a stem-to-stem 2-ident., then there is one more restriction.
        #Or if v is a 2-vtx without any stem identification, then there are two more restrictions.
        if v in twos:
            if v in ids2:
                count += 1
            elif not v in ids1|ids3:
                count += 2
        #The f-value for v is 7 minus the number of restrictions from outside the graph.
        f.append(7-count)
    
    #Now for the graph.  For a stem-to-root 2-identification, the vertices are adjacent and affect the square.
    #If two vertices have been stem-to-root identified, then there is no edge between them in underlying_graph.
    G = Graph(underlying_graph)
    for e in identifications[0]:
        G.add_edge(e)
    #Now square the resulting graph.
    new_G = graph_power(G,2)
    #The only changes left to be made are edges which might show up from stem-to-stem identifications.
    for e in identifications[1]:
        #We assume that stem-to-stem 2-identifications are done with roots at distance at least 3 from each other.
        new_G.add_edge(e)
    for t in identifications[2]:
        #Stem-to-stem 3-identifications might happen among vertices of which some of which are at distance 2.
        for e in [t[:2],[t[0],t[2]],t[1:]]:
            if not new_G.has_edge(e):
                new_G.add_edge(e)

    return new_G,f





def length_equals_8():
    # The case ell(F) = 8.

    # The graph is a three cycle (on 0,1,2) connected to another 3-cycle (on 3,4,5) by an edge (2,3).
    g = Graph([(0,1),(0,2),(1,2),(2,3),(3,4),(3,5),(4,5)])

    # Recall that in this lemma, new edges and stems must occur within the two regions separated by F.
    # There are only three neighborhood structures, up to symmetry.
    # Each neighborhood structure is a list of three lists:  pairs of vertices that form a new edge, pairs of vertices that have a common stem with only each other, and triples of vertices that have a common stem.
    neighborhood_structures = [ 
        [[],[],[]],
        [[],[[0,1]],[]],
        [[],[[0,1],[4,5]],[]]
        ]

    print "Checking the case for length(F) = 8."
    for ns in neighborhood_structures:
        core_square,f = core_square_graph(g,ns)
        print "Checking the neighborhood structure",ns
        ch.fChoosable(core_square,f)
        print
    print "\n"*4





def length_equals_9():
    # The case ell(F) = 9.

    # The graph is a three cycle (on 0,1,2) connected to a 4-cycle (on 3,4,5,6) by an edge (2,3).
    g = Graph([(0,1),(0,2),(1,2),(2,3),(3,4),(4,5),(5,6),(3,6)])

    # Recall that in this lemma, new edges and stems must occur within the two regions separated by F.
    # There are only ten neighborhood structures, up to symmetry (two within the 3-cycle, times five within the 4-cycle).
    # Each neighborhood structure is a list of three lists:  pairs of vertices that form a new edge, pairs of vertices that have a common stem with only each other, and triples of vertices that have a common stem.
    neighborhood_structures = [ 
        [[],[],[]],
        [[],[[0,1]],[]],
        [[[4,6]],[],[]],
        [[[4,6]],[[0,1]],[]],
        [[],[[4,6]],[]],
        [[],[[4,6],[0,1]],[]],
        [[],[[4,5]],[]],
        [[],[[4,5],[0,1]],[]],
        [[],[],[[4,5,6]]],
        [[],[[0,1]],[[4,5,6]]]
        ]
    
    print "Checking the case for length(F) = 9."
    for ns in neighborhood_structures:
        core_square,f = core_square_graph(g,ns)
        print "Checking the neighborhood structure",ns
        ch.fChoosable(core_square,f)
        print
    print "\n"*4

    
    





if __name__ == "__main__":
    length_equals_8()
    length_equals_9()





