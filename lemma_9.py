# To-Do:
    # Possible change:  Collapse all cases down to noninduced cases, then prove with CNS?
    # After paper is finished, check references to Lemmas/Observations/Definitions.


#------------------------------------------------------
# Contents
#------------------------------------------------------

#    length_equals_8
#    length_equals_9
#    (main command)



#------------------------------------------------------
# Libraries, Comments, etc.
#------------------------------------------------------

# These scripts were written with the intent of using Python 2.7 before the switch to Python 3.  There should not be many differences to implement our code in Python 3 once this change is made.

from checking_realizations import *
# We use checking_realizations.py as a module which contains the common code for Lemmas ??, ??, and ??.

import choosability as ch
# We use fChoosableNoPrint as a black box.



#------------------------------------------------------
# The Routines
#------------------------------------------------------



def length_equals_8():
    # The case ell(F) = 8.
    
    # The graph is a three cycle (on 0,1,2) connected to another 3-cycle (on 3,4,5) by an edge (2,3).
    edges = [(0,1),(0,2),(1,2),(2,3),(3,4),(3,5),(4,5)]
    
    # Recall that in this lemma, new edges and stems must occur within the two regions separated by F.
    # There are only three stem structures, up to symmetry.
    # Each neighborhood structure is a list of three lists:  pairs of vertices that form a new edge, pairs of vertices that have a common stem with only each other, and triples of vertices that have a common stem.
    stem_structures = [ 
        [[],[],[]],
        [[],[{0,1}],[]],
        [[],[{0,1},{4,5}],[]]
        ]
    
    print "*** Checking the case for length(F) = 8. ***"
    print "The plane graph in question is given by the edges "+str(edges)+".\n"
    for SS in stem_structures:
        core_square,f = core_square_graph(Graph(edges),SS)
        print "   >>> Checking the neighborhood structure "+str(SS)+"."
        print "      >>> The induced subgraph of the square on the core vertices has the following edge set:"
        print "          "+str([e[:2] for e in core_square.edges()])
        print "      >>> The list size function is given by the following list:"
        print "          "+str(f)
        x = ch.fChoosableNoPrint(core_square,f)
        if x[0]:
            print "      >>> Verified ("+x[1]+") -- the induced subgraph of the square is choosable for the given list size function."
        else:
            print "      >>> !PROBLEM! -- the induced subgraph of the square is NOT choosable for the given list size function."
        print
    print "\n"*4





def length_equals_9():
    # The case ell(F) = 9.
    
    # The graph is a three cycle (on 0,1,2) connected to a 4-cycle (on 3,4,5,6) by an edge (2,3).
    edges = [(0,1),(0,2),(1,2),(2,3),(3,4),(4,5),(5,6),(3,6)]
    
    # Recall that in this lemma, new edges and stems must occur within the two regions separated by F.
    # There are only ten neighborhood structures, up to symmetry (two within the 3-cycle, times five within the 4-cycle).
    # Each neighborhood structure is a list of three lists:  pairs of vertices that form a new edge, pairs of vertices that have a common stem with only each other, and triples of vertices that have a common stem.
    stem_structures = [ 
        [[],[],[]],
        [[],[{0,1}],[]],
        [[{4,6}],[],[]],
        [[{4,6}],[{0,1}],[]],
        [[],[{4,6}],[]],
        [[],[{4,6},{0,1}],[]],
        [[],[{4,5}],[]],
        [[],[{4,5},{0,1}],[]],
        [[],[],[{4,5,6}]],
        [[],[{0,1}],[{4,5,6}]]
        ]
    
    print "*** Checking the case for length(F) = 9. ***"
    print "The plane graph in question is given by the edges "+str(edges)+".\n"
    for SS in stem_structures:
        core_square,f = core_square_graph(Graph(edges),SS)
        print "   >>> Checking the neighborhood structure "+str(SS)+"."
        print "      >>> The induced subgraph of the square on the core vertices has the following edge set:"
        print "          "+str([e[:2] for e in core_square.edges()])
        print "      >>> The list size function is given by the following list:"
        print "          "+str(f)
        x = ch.fChoosableNoPrint(core_square,f)
        if x[0]:
            print "      >>> Verified ("+x[1]+") -- the induced subgraph of the square is choosable for the given list size function."
        else:
            print "      >>> !PROBLEM! -- the induced subgraph of the square is NOT choosable for the given list size function."
        print
    print "\n"*4








if __name__ == "__main__":
    length_equals_8()
    length_equals_9()





