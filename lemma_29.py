# To-Do:
    # Separate expanded c7a4x5x5 stuff into its own sage script.
    # After paper is finished, check references to Lemmas/Observations/Definitions.


#------------------------------------------------------
# Contents
#------------------------------------------------------

##    check_all_realizations_from_expanded_c7a4x5x5_case
#    run_lemma_29
#    (main command)



#------------------------------------------------------
# Libraries, Comments, Global Variables
#------------------------------------------------------

# These scripts were written with the intent of using Python 2.7 before the switch to Python 3.  There should not be many differences to implement our code in Python 3 once this change is made.

from checking_realizations import *
# We use checking_realizations.py as a module which contains the common code for Lemmas ?? and ??.

#import time
# For displaying runtime.



#------------------------------------------------------
# The Routines
#------------------------------------------------------







###Fix to work with NaturalCoreSubgraph class.
def check_all_realizations_from_expanded_c7a4x5x5_case(case):
    print "Expanding scope around c7a4x5x5:"

    order,edges,spine,roots,faces = makeGraph('c7a4x5x5')#spine not modified here.
    open_spine = False
    roots.sort()#why sort?
    
    #switch partition_di to forbidden_dict, remove stem identification dictionaries that don't affect partition restrictions.
    
    order = 15  # Number of vertices in the natural core subgraph of 7:4*5*5.

    #Initial formulations of the restrictions:  no root identifications, and no stem identifications 
    #except non-edge [6,8] and triple [6,8,12].
    partition_di = {x:set(range(order)) for x in range(order)}
    stem_1 = {x:set(roots)-{x,} for x in roots}
    stem_2 = {x:set(roots)-{x,} for x in roots}
    stem_2[6].remove(8)
    stem_2[8].remove(6)
    stem_3 = {x:{frozenset([roots[y],z]) for y in range(len(roots)-1) if not roots[y]==x for z in roots[y+1:] if not z==x} for x in roots}
    stem_3[6].remove(frozenset([8,12]))
    stem_3[8].remove(frozenset([6,12]))
    stem_3[12].remove(frozenset([6,8]))

    #Even for non-roots, we want to have the stem_restrictions defined in check_all_realizations_from_initial_plane_graph.
    for v in [x for x in range(order) if not x in roots]:
        stem_1[v] = set()
        stem_2[v] = set()
        stem_3[v] = set()
    
    

    #Now, which case are we in?
    #A \ne 3 by identification case.
    #A \ne 4 by identification case.
    #A \ne 5 by 4:57.
    if case == 1:
        print "Case A = 6."
        le = 6#Length of new face.
        base_border = [7,0,6]#In appropriate rotation for new face.
    
    #B \ne 3 by identification case.
    #B \ne 4 by 44.
    elif case == 2:
        print "B = 5."
        le = 5#Length of new face.
        base_border = [8,7]#In appropriate rotation for new face.
    elif case == 3:
        print "B = 6."
        le = 6#Length of new face.
        base_border = [8,7]#In appropriate rotation for new face.
    
    #C/D \ne 3 by identification case.
    #C/D \ne 4 by identification case.
    #C/D \ne 5 by identification case.
    elif case == 4:
        print "C/D = 6."
        le = 6#Length of new face.
        base_border = [9,2,1,8]#In appropriate rotation for new face.
    
    #E \ne 3 by identification case.
    elif case == 5:
        print "E = 4."
        le = 4#Length of new face.
        base_border = [10,9]#In appropriate rotation for new face.
    elif case == 6:
        print "E = 5."
        le = 5#Length of new face.
        base_border = [10,9]#In appropriate rotation for new face.
    elif case == 7:
        print "E = 6."
        le = 6#Length of new face.
        base_border = [10,9]#In appropriate rotation for new face.
    
    #G/H \ne 3 by identification case.
    #G/H \ne 4 by identification case.
    #G/H \ne 5 by identification case.
    elif case == 8:
        print "G/H = 6."
        le = 6#Length of new face.
        base_border = [12,4,3,11]#In appropriate rotation for new face.
    
    #I \ne 3 by identification case.
    elif case == 9:
        print "I = 4."
        le = 4#Length of new face.
        base_border = [13,12]#In appropriate rotation for new face.
    elif case == 10:
        print "I = 5."
        le = 5#Length of new face.
        base_border = [13,12]#In appropriate rotation for new face.
    elif case == 11:
        print "I = 6."
        le = 6#Length of new face.
        base_border = [13,12]#In appropriate rotation for new face.
    
    #K \ne 3 by identification case.
    #K \ne 4 by identification case.
    elif case == 12:
        print "K = 5."
        le = 5#Length of new face.
        base_border = [6,5,14]#In appropriate rotation for new face.
    elif case == 13:
        print "K = 6."
        le = 6#Length of new face.
        base_border = [6,5,14]#In appropriate rotation for new face.

    else:
        print "Out of range!"
        
    print



    new_num = le - len(base_border)
    old_order = order

    i = 0
    edges.append((base_border[-1],order))
    i += 1
    while i < new_num:
        edges.append((order,order+1))
        i += 1
        order += 1
    edges.append((base_border[0],order))
    order += 1

    del roots[roots.index(base_border[0])]#switch roots to set?
    del roots[roots.index(base_border[-1])]
    roots.extend(range(old_order,order))

    new_face = copy.copy(base_border)###avoid copy?  Do a loop if necessary
    new_face.extend(range(old_order,order))
    faces.append(new_face)

#     Graph(edges).show()
    g = (order,edges,spine,roots,faces)
#     print g
    #Initialize dictionaries for the new vertices:
    for v in range(old_order,order):
        partition_di[v] = set(faces[-1])
        #We could further restrict things that new vertices can identify with, but we're not going to.
        stem_1[v] = set()
        stem_2[v] = set()
        stem_3[v] = set()
    #Actually, we'll go ahead and add some basic restrictions: new vertices adjacent to the base can't 
    #be identified with vertices from the original that the base neighbors couldn't have an edge to, 
    #because this will make such an edge.
    partition_di[old_order] |= stem_1[base_border[-1]]
    partition_di[order-1] |= stem_1[base_border[0]]
    #Also, new vertices adjacent to the base can't be make an edge with vertices from the original that
    #the base neighbors couldn't have a 2-stem with, because this will make such a 2-stem.
    stem_1[old_order] |= stem_2[base_border[-1]]
    stem_1[order-1] |= stem_2[base_border[0]]
    #Also, new vertices distance 2 from the base can't be identified with vertices from the original 
    #that the base neighbors (distance 2) couldn't have a 2-stem with, because this will make such a 2-stem.
    if old_order - order > 1:
        partition_di[old_order+1] |= stem_2[base_border[-1]]
        partition_di[order-2] |= stem_2[base_border[0]]
    #Again, we could do more restrictions.  But it would probably not be worth the savings for such small additional faces.



    check_all_realizations_from_initial_plane_graph(NCS,forbid_identifications=partition_di)











###adjust for lemma 29
def run_lemma_20():
    # Verifies that each configuration in the target set is reducible by checking that all realizations are core-choosable.
    
    begin = time.clock()
    
    for rc in TargetSet:
        check_all_realizations_of_configuration(rc)
        print "\n"*20
    
    print "Finished with all configurations in the target set!"
    print "Total time: "+timestring(time.clock()-begin)










###adjust for lemma 29
if __name__ == "__main__":
    run_lemma_20()




