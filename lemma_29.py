# To-Do:
    # After paper is finished, check references to Lemmas/Observations/Definitions.


#------------------------------------------------------
# Contents
#------------------------------------------------------

#    check_all_realizations_from_expanded_c7a4x5x5_case
#    run_lemma_29
#    (main command)



#------------------------------------------------------
# Libraries, Comments, Global Variables
#------------------------------------------------------

# These scripts were written with the intent of using Python 2.7 before the switch to Python 3.  There should not be many differences to implement our code in Python 3 once this change is made.

from checking_realizations import *
# We use checking_realizations.py as a module which contains the common code for Lemmas ?? and ??.

#------------------------------------------------------
# The Routines
#------------------------------------------------------








def check_all_realizations_from_expanded_c7a4x5x5_case(case):
    # This takes an int between 1 and 12 corresponding to a particular case from Lemma ??.  It then builds the corresponding plane graph and dictionary of forbidden vertex identifications, and then checks all resulting realizations for core-choosability.
    
    
    print "Expanding scope around c7a4x5x5:"
    
    # We begin with the natural core subgraph of c7a4x5x5, but we later add a face and adjust PlaneGraph accordingly.
    PlaneGraph = NaturalCoreSubgraph('c7a4x5x5')
    
    # None of the vertices in the original core subgraph will be permitted to be identified with each other.
    forbidden_dict = {x:set(PlaneGraph.graph.vertices()) for x in PlaneGraph.graph.vertices()}
    
    
    # Now we manually preprocess our adjustments to PlaneGraph according to each of the twelve remaining cases.
    # new_face_length will be the length of the new face, and base_border will be the clockwise ordering of the part of the boundary of the new face which is already in PlaneGraph.
    
    # We begin with the faces which are necessarily adjacent to the central 7-face:  A, C, G, K.
    if case == 1:
        print "Case A = 6.\n"
        new_face_length = 6
        base_border = [7,0,6]
        PlaneGraph.name += "CaseA6"
    
    elif case == 2:
        print "C/D = 6.\n"
        new_face_length = 6
        base_border = [9,2,1,8]
        PlaneGraph.name += "CaseC6"
    
    elif case == 3:
        print "G/H = 6.\n"
        new_face_length = 6
        base_border = [12,4,3,11]
        PlaneGraph.name += "CaseG6"
    
    elif case == 4:
        print "K = 6.\n"
        new_face_length = 6
        base_border = [6,5,14]
        PlaneGraph.name += "CaseK6"
    
    # Then we move on to the other face adjacent to the 4-face:  B.
    elif case == 5:
        print "B = 5.\n"
        new_face_length = 5
        base_border = [8,7]
        PlaneGraph.name += "CaseB5"
    elif case == 6:
        print "B = 6.\n"
        new_face_length = 6
        base_border = [8,7]
        PlaneGraph.name += "CaseB6"
    
    # Then we move on to another face adjacent to the first 5-face:  E.
    elif case == 7:
        print "E = 4.\n"
        new_face_length = 4
        base_border = [10,9]
        PlaneGraph.name += "CaseE4"
    elif case == 8:
        print "E = 5.\n"
        new_face_length = 5
        base_border = [10,9]
        PlaneGraph.name += "CaseE5"
    elif case == 9:
        print "E = 6.\n"
        new_face_length = 6
        base_border = [10,9]
        PlaneGraph.name += "CaseE6"
    
    # Then we move on to another face adjacent to the second 5-face:  I.
    elif case == 10:
        print "I = 4.\n"
        new_face_length = 4
        base_border = [13,12]
        PlaneGraph.name += "CaseI4"
    elif case == 11:
        print "I = 5.\n"
        new_face_length = 5
        base_border = [13,12]
        PlaneGraph.name += "CaseI5"
    elif case == 12:
        print "I = 6.\n"
        new_face_length = 6
        base_border = [13,12]
        PlaneGraph.name += "CaseI6"
    
    
    # Now that we have the beginning of our new face, we can build the rest of it and change PlaneGraph accordingly.
    PlaneGraph.faces.append(base_border)  # We will extend the new face from base_border.
    num_new_vtxs = new_face_length - len(base_border)
    
    # Build the first edge out from the original plane graph.
    PlaneGraph.edges.append((base_border[-1],PlaneGraph.order))  # base_border has not yet been changed, so base_border[-1] is our starting point.  PlaneGraph.order is the name of the new vertex.
    PlaneGraph.faces[-1].append(PlaneGraph.order)  # Add the new vertex to our new face.
    PlaneGraph.order += 1
    
    # Build all remaining edges of the new face except the last closing edge.
    for i in range(num_new_vtxs-1):  # We already built one new vertex.
        PlaneGraph.edges.append((PlaneGraph.order-1,PlaneGraph.order))  # PlaneGraph.order is the name of the new vertex.
        PlaneGraph.faces[-1].append(PlaneGraph.order)  # Add the new vertex to our new face.
        PlaneGraph.order += 1
    
    # Build the last closing edge of the new face.
    PlaneGraph.edges.append((PlaneGraph.faces[-1][0],PlaneGraph.faces[-1][-1]))
    
    
    # Complete forbidden_dict.
    for v in range(15,15+num_new_vtxs):  # Initialize forbidden_dict for the new vertices.
        forbidden_dict[v] = set()
    for v in PlaneGraph.faces[-1]:  # For all the vertices on the new face, add the face to the forbidden set.
        forbidden_dict[v] |= set(PlaneGraph.faces[-1])
    
    
    # There is no need to update anything else from PlaneGraph before feeding into the checker, since it will use only .name, .order, .is_open (unchanged), .faces, and .edges.
    check_all_realizations_from_initial_plane_graph(PlaneGraph,forbid_identifications=forbidden_dict)










###add normal run for c7a4x5x5 first -- need to put check configuration routine into common pool.
def run_lemma_29():
    # Verifies that each configuration in the target set is reducible by checking that all realizations are core-choosable.
    
    begin = time.clock()
    
    # First, check the original configuration to see the two bad realizations.
    check_all_realizations_of_configuration('c7a4x5x5')
    print "\n"*10
    
    # Now check the remaining 12 cases of faces around the 4- and 5-faces.
    for i in range(1,13):
        check_all_realizations_from_expanded_c7a4x5x5_case(i)
        print "\n"*10
    
    print "Finished with all analysis of c7a4x5x5!"
    print "Total time: "+timestring(time.clock()-begin)










if __name__ == "__main__":
    run_lemma_29()




