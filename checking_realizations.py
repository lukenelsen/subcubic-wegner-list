# To-Do:
    # Make vertex coordinates for tikz.
    # Adjust the print-outs for our wants.  (remove last comma in edge list)
    # After paper is finished, check references to Lemmas/Observations/Definitions.


#------------------------------------------------------
# Contents
#------------------------------------------------------

# Section 1:  Graph Essentials
#    NaturalCoreSubgraph
#    core_square_graph

# Section 2:  Core Subgraph Generation
#    restricted_partition_generator
#    core_subgraph_generator

# Section 3:  Stem Structure Generation
#    get_open_region_roots
#    stem_structure_generator

# Section 4:  Checking Realizations
#    check_all_stem_structures_for_given_core_subgraph
#    check_all_realizations_from_initial_plane_graph



#------------------------------------------------------
# Libraries, Comments, Global Variables
#------------------------------------------------------

# These scripts were written with the intent of using Python 2.7 before the switch to Python 3.  There should not be many differences to implement our code in Python 3 once this change is made.

from sage.graphs.graph import Graph
# We use the Graph class to work with our realizations.  We use it for neighborhoods (core_square_graph), the distance matrix (core_square_graph), and checking planarity (core_subgraph_generator).

from itertools import combinations
# We use combinations as a pairs generator in stem_structure_generator when we are iterating through triples containing the last root, and also in core_square_graph when we are adding edges resulting from these triples.

import choosability as ch
# We use fChoosableNoPrint as a black box to verify core-choosability for all our realizations.

import time
# For displaying runtime.

def timestring(time):
    # Displays time length:  ?m ?.?s
    
    t = int(time)
    d = str(int((10*(time-t))%10))
    s = str(t%60)
    m = int(t/60)
    string = str(m)+"m "+s+"."+d+"s"
    return string

ShortTargetSet =[
    # open chain:  start with 3 and end with 3
    'c3a3',
    'cxa3x3',
    'cxa3xx3',
    'cxa3xxx3',
    'cxa3xx5x3',

    # open chain:  start with 3 and end with 4
    'c3a4',
    'cxa3x4',
    'cxa35x4',
    'cxa3x54',
    'cxa3xx54',
    'cxa3xxx54',

    # open chain:  start with 3 and end with 5
    'cxa355',
    'cxa375',
    'cxa3x555',

    # open chain:  start with 3 and end with 6
    'c3a6',
    'cxa356',

    # open chain:  start with 4 and end with 4
    'c4a4',
    'cxa454',
    'cxa464',
    'cxa474',
    'cxa4x54',
    #'cxa4x55x4',
    #'cxa4x555x4',
    #'cxa4x4x4x4',

    # open chain:  start with 4 and end with 5
    'cxa4x45',
    'cxa4555',
    #'cxa4x5555',

    # open chain:  start with 4 and end with 6
    'cxa456',

    # open chain:  start with 5 and end with 5
    'cxa535',
    #'cxa555555',

    # open chain:  start with 5 and end with 6
    'cxa546',

    # central 3-face
    'c3a57',

    # central 4-face
    'c4a55',
    'c4a56',
    'c4a57',
    'c4a66',
    'c4a585',
    'c4a676',
    'c4a686',

    # central 5-face
    'c5a555',
    'c5a556',
    'c5a565',
    'c5a566',
    'c5a575',
    'c5a656',
    'c5a666',
    'c5a4x55',    
    'c5a5x65',
    #'c5a5585',
    #'c5a5x66',
    'c5a55x6',
    #'c5a56x6',
    #'c5a6x66',

    # central 7-face
    'c7a3xx4',
    'c7a3xx5',
    'c7a4xx4',
    'c7a4x55',
    'c7a4xx55',
    'c7a4x5xx5',
    #'c7a55x5',

    # central 8-face
    'c8a3xx4',
    'c8a3xxx4',
    'c8a35x5',
    'c8a35xx5',
    'c8a35xxx5',
    'c8a35xxxx5',
    #'c8a3x55x55',
    'c8a45xx4',
    #'c8a455x5',
    #'c8a4x5x4x4',

    # central 9-face
    'c9a3xx4',
    'c9a3xxx4',
    'c9a35x5'
    #'c9a455x4',
    #'c9a455x5',
    #'c9a545x5',
        ]



# We store the target set as a list of strings.  Each string follows a slightly different format than our standard description in the paper.  The configuration b:s_1...s_t is encoded as "cbas_1...s_t" where any "*" is replaced with "x".  (The "c" precedes the (c)entral face and the "a" precedes the facial length list (a)round the central face.)  For example, "9:3**4" is encoded as "c9a3xx4" and "4*555*4" is encoded as "cxa4x555x4".

TargetSet =[
    # open chain:  start with 3 and end with 3
    'c3a3',  # In the paper, we list this configuration as cxa33.
    'cxa3x3',
    'cxa3xx3',
    'cxa3xxx3',
    'cxa3xx5x3',

    # open chain:  start with 3 and end with 4
    'c3a4',  # In the paper, we list this configuration as cxa34.
    'cxa3x4',
    'cxa35x4',
    'cxa3x54',
    'cxa3xx54',
    'cxa3xxx54',

    # open chain:  start with 3 and end with 5
    'cxa355',
    'cxa375',
    'cxa3x555',

    # open chain:  start with 3 and end with 6
    'c3a6',  # In the paper, we list this configuration as cxa36.
    'cxa356',

    # open chain:  start with 4 and end with 4
    'c4a4',  # In the paper, we list this configuration as cxa44.
    'cxa454',
    'cxa464',
    'cxa474',
    'cxa4x54',
    'cxa4x55x4',
    'cxa4x555x4',
    'cxa4x4x4x4',

    # open chain:  start with 4 and end with 5
    'cxa4x45',
    'cxa4555',
    'cxa4x5555',

    # open chain:  start with 4 and end with 6
    'cxa456',

    # open chain:  start with 5 and end with 5
    'cxa535',
    'cxa555555',

    # open chain:  start with 5 and end with 6
    'cxa546',

    # central 3-face
    'c3a57',

    # central 4-face
    'c4a55',
    'c4a56',
    'c4a57',
    'c4a66',
    'c4a585',
    'c4a676',
    'c4a686',

    # central 5-face
    'c5a555',
    'c5a556',
    'c5a565',
    'c5a566',
    'c5a575',
    'c5a656',
    'c5a666',
    'c5a4x55',    
    'c5a5x65',
    'c5a5585',
    'c5a5x66',
    'c5a55x6',
    'c5a56x6',
    'c5a6x66',

    # central 7-face
    'c7a3xx4',
    'c7a3xx5',
    'c7a4xx4',
    'c7a4x55',
    'c7a4xx55',
    'c7a4x5xx5',
    'c7a55x5',

    # central 8-face
    'c8a3xx4',
    'c8a3xxx4',
    'c8a35x5',
    'c8a35xx5',
    'c8a35xxx5',
    'c8a35xxxx5',
    'c8a3x55x55',
    'c8a45xx4',
    'c8a455x5',
    'c8a4x5x4x4',

    # central 9-face
    'c9a3xx4',
    'c9a3xxx4',
    'c9a35x5',
    'c9a455x4',
    'c9a455x5',
    'c9a545x5',
        ]



#------------------------------------------------------
# The Routines
#------------------------------------------------------


#------------------------------------------
# Section 1:  Graph Essentials
#------------------------------------------





class NaturalCoreSubgraph:
    # This class is for holding all of the necessary information we need to know from an initial plane graph.  For Lemma ??, this means the natural core subgraph of a configuration.  For Lemma ??, we will start with the natural core subgraph of c7a4x5x5 and edit the plane graph for our purposes.
    # The constructor takes the string following the notation described before the definition of TargetSet above and then builds the graph.
    # Note:  We assume for any configuration that we feed into the constructor, that any specified entry preceding or succeeding a 3 in the facial length list has length at least 5.  (So cxa33 is not permissible, but c3a3 is.)  We also assume that the facial length list around a central face of specified length is shorter then the central face length.  (So c4a85x5 is not permissible, but c4a585 is.)
    # The following are the attributes of an object from this class:
        # .name -- the string originally given in c/a notation
        # .is_open -- a boolean, True if the configuration is open
        # .central_face_length -- the integer value of the specified central face length, or None otherwise
        # .facial_length_list -- the list of entries of the facial length list around the central face
        # .edges -- the list of edges (as tuples) of the plane graph
        # .faces -- the list of clockwise-oriented boundaries of the spine and/or specified faces.  The first item is always the known boundary of the central face.  (If the configuration is open, this is the spine.  Otherwise, it is all of the central face.)  The remaining items are the noncentral specified faces, as ordered by self.facial_length_list.  Each boundary is itself a list of vertices.
        # .graph -- the initial plane graph as a Graph object
        # .order -- the number of vertices in the initial plane graph
        # .spine -- the actual spine given as a list of vertices.  If the central face has specified length, then spine is not necessarily the entirety of the central face.  (So spine is not always just self.faces[0].)
    
    def __init__(self,config_str):
        
        self.name = config_str
        
        # First, look at the characters between 'c' and 'a' to determine if the configuration is open or closed.
        s = config_str[1:config_str.index('a')]
        if s == 'x':  # Open (unspecified central face length).
            self.is_open = True
            self.central_face_length = None
        else:  # Closed (specified central face length).
            self.is_open = False
            self.central_face_length = int(s)
        
        # Then, look at the facial length list around the central face.
        self.facial_length_list = []
        for s in config_str[config_str.index('a')+1:]:
            if s == 'x':
                self.facial_length_list.append('x')
            else:
                self.facial_length_list.append(int(s))
        FLL = self.facial_length_list  # 'FLL' will be easier to read in the rest of this constructor than 'self.facial_length_list'.
        
        # Now we build the natural core subgraph by building our list of edges,.  As we build the edges, we also keep track of the central face/spine and the specified faces.
        self.edges = []  # Each entry will be a tuple of two endpoints.
        self.faces = []  # Each entry will be a list of vertices.
        
        # We begin with the central face.  There are two possibilities:  open or closed.  In the open case, the spine will be our central face.  In our later generation, we will need to treat the central face specially to account for this.
        # In the open case, we build a path and append this to faces.
        if self.is_open:
            for j in range(len(FLL)):  # Build the path from 0 to len(FLL) (inclusive).
                self.edges.append((j,j+1))
            self.faces.append(range(len(FLL)+1))  # Add the path to faces.
            vtx = len(FLL)  # vtx is the current vertex we are building edges to.  (For future use.)
        # In the closed case, we build a cycle and append this to faces.
        else:
            for j in range(self.central_face_length-1):  # Build the path from 0 to self.central_face_length-1 (inclusive).
                self.edges.append((j,j+1))
            self.edges.append((0,self.central_face_length-1))  # Close the path to make it a cycle.
            self.faces.append(range(self.central_face_length))  # Add the cycle to faces.
            vtx = self.central_face_length-1  # vtx is the current vertex we are building edges to.  (For future use.)
        
        # Now that we have built the central face, we move on to the other specified faces.
        for i in range(len(FLL)):
            # i is the index of the entry in FLL, but it is also the counterclockwise-most vertex on the corresponding edge.  i+1 is the clockwise-most vertex of that edge.
            
            # If the ith entry of the facial length list is 'x', then we do nothing.
            if FLL[i]=='x':
                continue
            
            # Otherwise, we want to build a face by building a path from i to i+1 of length FLL[i]-1.
            # Simultaneously, we build the list of the face's boundary.
            new_face = [i+1,i]
            # First, we check the previous entry(ies) to see how many shared edges have already been built.  There are three cases:  (i) the previous entry is not a specified face, (ii) the previous entry is a sandwiched 3-face, or (iii) the previous entry is a specified face that is not a sandwiched 3-face.
            # Based on these cases, we define the number of shared edges already built and the starting point for the rest of the path to i+1.
            if i==0 or FLL[i-1]=='x':
                num_shared = 0
                start = i
            elif i>1 and FLL[i-1]==3 and FLL[i-2]!='x':  # This only occurs for cxa535.  We could have used the alternative encoding c3a55 but prefer the former.
                num_shared = 2
                start = vtx-1  # vtx-1 is the second-to-last vertex that was created for a face.  This is the right vertex because the face previous to the 3-face has length at least 5.
                new_face.append(vtx)
                new_face.append(vtx-1)
            else:  # Previous entry is a specified face but not a sandwiched 3-face.
                num_shared = 1
                start = vtx  # vtx is the last vertex that was created.
                new_face.append(vtx)
            
            # Now we build any remaining edges except for the very last one.
            # FLL[i]-2-num_shared is the length of the face we are building minus:  the edge from the central face, the last edge we will build after this loop, and the edges already built by previous faces.
            for j in range(FLL[i]-2-num_shared):
                vtx += 1
                self.edges.append((start,vtx))
                new_face.append(vtx)
                start = vtx
            
            # Now we build the last edge.
            self.edges.append((i+1,vtx))
            
            # Note:  The last vertex/edge of the last specified face won't overlap with the first specified face, because we assume that configurations don't have full wrap-around.
            
            # The face is now complete.
            self.faces.append(new_face)
            
        # Add a few other things.
        self.graph = Graph(self.edges)
        self.order = self.graph.order()
        self.spine = range(len(FLL)+1)











def core_square_graph(core_subgraph,stem_structure):
    # core_subgraph is a Sage Graph object of the core subgraph of the realization, and stem_structure encodes the stem structure of the realization.  stem_structure is a list of three lists:  first a list of root-to-root edges (as sets), second a list of root pairs connected to stems of degree 2 (as sets), and third a list of root triples connected to stems of degree 3 (as sets).
    # core_square_graph takes the realization given by core_subgraph and stem_structure and returns the graph and list size function described in Definition ?(core-choosable)?.
    
    
    # Start with making the graph (the induced subgraph of the square of the host graph on the set of core vertices).  First, we add the root-to-root edges from the stem structure.  Then, we square that resulting graph.  Finally, we add any additional edges which would live in the square of the host graph as a result of the stems of degree 2 or 3 in the stem structure.
    new_G = Graph(core_subgraph)  # We make a copy since core_subgraph will be needed for other stem_structures.
    
    for edge in stem_structure[0]:  # Root-to-root edges from the stem structure.
        new_G.add_edge(edge)  # If edge is already an edge of the graph, then this does nothing.
    
    # The current graph is H[A], where H is the realization and A is the set of core vertices.  We save this graph for later when we calculate the list size function.
    H_A = Graph(new_G)
    
    # Now square the resulting graph.
    M = new_G.distance_matrix()
    for i,j in combinations(range(new_G.order()),2):
        if M[i][j] == 2:
            new_G.add_edge((i,j))

    # Finally, add any additional edges which would live in the square of the host graph.  (The current graph is (H[A])^2, but we might need to add some edges to get H^2[A].)
    for edge in stem_structure[1]:  # Roots connected to stems of degree 2 from the stem structure.
        new_G.add_edge(edge)  # If edge is already an edge of the graph, then this does nothing.
    for triple in stem_structure[2]:  # Roots connected to stems of degree 3 from the stem structure.
        for edge in combinations(triple,2):
            new_G.add_edge(edge)  # If edge is already an edge of the graph, then this does nothing.
    
    
    # Now make f, the list size function.  For each core vertex, we calculate the restriction index defined in Definition ?(core-choosable)?.  For every stem, we first see if the vertex is adjacent to the stem or if it is distance 2 from the stem.  
    f = []
    
    # We need to make note of which roots are attached to stems which have degree 1 to the core.
    roots = {x for x in core_subgraph.vertices() if core_subgraph.degree(x)==2}
    roots_already_observed = {x for stem_type in stem_structure for stem_group in stem_type for x in stem_group}  # "for" statements in a set comprehension are processed from left to right.  This set contains all the vertices appearing anywhere in stem_structure.
    one_stem_roots = roots - roots_already_observed  # Vertices that have not appeared anywhere in stem_structure are roots which are adjacent to a stem which has degree 1 to the core.
    
    for v in H_A.vertices():  # This iterates through the vertices in order:  0,...,n-1.
        restriction_index = 0
        v_nbrs = set(H_A.neighbors(v))  # v_nbrs is the set of neighbors of v in H[A].
        
        for one_stem in one_stem_roots:
            if v == one_stem:  # Then gamma_1(v) == 1.
                restriction_index += 3
            elif one_stem in v_nbrs:  # Then Gamma(v) gains 1 through one_stem.
                restriction_index += 1
        
        for two_stem in stem_structure[1]:
            if v in two_stem:  # Then gamma_2(v) == 1.
                restriction_index += 2
            elif two_stem & v_nbrs:  # Then Gamma(v) gains 1 through two_stem.
                restriction_index += 1
        
        for three_stem in stem_structure[2]:
            if v in three_stem:  # Then gamma_3(v) == 1.
                restriction_index += 1
            elif three_stem & v_nbrs:  # Then Gamma(v) gains 1 through two_stem.
                restriction_index += 1
        
        f.append(7-restriction_index)

    return new_G,f
















#------------------------------------------
# Section 2:  Core Subgraph Generation
#------------------------------------------




def restricted_partition_generator(li,last_index,forbidden_dict):
    # Given a list li of distinct objects, we generate all partitions of li which do not put certain objects into the same part.  The restrictions on which objects can be put into the same part are given by forbidden_dict, a dictionary which maps each object foo to the set of objects with which foo cannot share a part.  Each generated partition is a list of sets.
    # We recursively generate the partitions by conditioning on which part the last element goes into.  For efficiency, we avoid copying/slicing li and the partition.  We use last_index to control how much of li is being viewed, and edit the partition in place immediately before and immediately after it is yielded.
    # Note:  last_index is the index of the last element being viewed.
    # Note:  forbidden_dict should be symmetric in the sense that if bar in forbidden_dict[foo], then also foo in forbidden_dict[bar].
    # Note:  It might be the case that foo will be in its own set of forbidden objects, even though it is unsensible to forbid an object to be in its own part.  However, this causes no problems since we only compare distinct entries of li.
    
    # The only partition of the empty set is the empty set.
    if last_index < 0:  # The relevant portion of li is empty if and only if last_index < 0.
        yield []
        return
    
    last_element = li[last_index]
    
    # Otherwise, we look at the portion of li prior to last_element and consider any partitions of this smaller list.  For each of these, we look at the possible parts we can put last_element into.
    # Note:  We pass li and forbidden_dict without modifications because they still contain all the information we need.
    for partition in restricted_partition_generator(li,last_index-1,forbidden_dict):
        
        # First, which parts of this partition are eligible for the last_element?
        for i in range(len(partition)):  # i indexes the parts.
            if not partition[i] & forbidden_dict[last_element]:  # Check the eligible parts.
                partition[i].add(last_element)  # Add last_element to ith part.
                yield partition
                partition[i].remove(last_element)  # Return partition to previous state without last_element.
        
        # Lastly, we can always put last_element in a part by itself.
        partition.append(set([last_element]))
        yield partition
        partition.pop()  # Return partition to previous state without last_element.
    return
















def core_subgraph_generator(NCS,forbid_identifications={}):    
    # NCS is an object from the NaturalCoreSubgraph class which holds all the information we need.
    # forbid_identifications is a dictionary.  Any foo:bar entries (foo is a vertex, bar is a set) will forbid foo from being identified with any element of bar.  When running check_all_realizations_of_configuration, forbid_identifications is empty.  When running check_all_realizations_from_expanded_c7a4x5x5_case, we will have added entries to forbid_identifications.
    # We use Observation ? to make a dictionary of all the forbidden identifications, in addition to what was given by forbid_identifications.  Then, we generate all identifications (partitions of the vertex set) which avoid these forbidden identifications.  For each identification, we check each of the three conditions from Observation ??.  For any identifications that make it through the three conditions, we yield the graph formed by the identification and a list of its faces.
    
    
    # forbidden_dict will be the dictionary of all the forbidden identifications.
    # Initialize each vertex to have any empty set of forbidden vertices.
    forbidden_dict = {x:set() for x in range(NCS.order)}
    
    # If we have manually entered any forbidden identifications, go ahead and add those first.
    for x in forbid_identifications.keys():
        forbidden_dict[x] |= set(forbid_identifications[x])
    
    # If we have an open central face, then we apply Observation ? to the vertices on the spine before moving on to the specified faces.
    if NCS.is_open:
        init_spec_face = 1  # This will tell a later loop to skip the central face.
        
        # Vertices on the spine cannot be identified which are too close.
        # Handle vertices close to the ends separately.
        # Note:  for any open configurations in the target set, the spine has at least four vertices.
        forbidden_dict[NCS.spine[0]] |= set(NCS.spine[:3])
        forbidden_dict[NCS.spine[1]] |= set(NCS.spine[:4])
        forbidden_dict[NCS.spine[-2]] |= set(NCS.spine[-4:])
        forbidden_dict[NCS.spine[-1]] |= set(NCS.spine[-3:])
        # Now handle all middle vertices on the spine.
        for i in range(2,len(NCS.spine)-2):
            forbidden_dict[NCS.spine[i]] |= set(NCS.spine[i-2:i+3])
        
        # Furthermore, spine vertices which are joining two specified outside faces ('sandwiched' vertices) cannot be identified with anything else on the spine.
        sandwiched_vertices = {x for x in range(1,len(NCS.facial_length_list)) if NCS.facial_length_list[x-1]!='x' and NCS.facial_length_list[x]!='x'}
        for v in NCS.spine:
            if v in sandwiched_vertices:
                forbidden_dict[v] |= set(NCS.spine)
            else:
                forbidden_dict[v] |= set(sandwiched_vertices)
    else:  # If we have a specified central face length...
        init_spec_face = 0  # ...then we tell the next loop to include the central face.
    
    # Now, for all specified faces, we apply Observation ??.  No two vertices on the same specified face can be identified.
    for face in NCS.faces[init_spec_face:]:
        for v in face:
            forbidden_dict[v] |= set(face)
    
    
    
    # Now that we have finished prepared forbidden_dict, we feed it to the partition generator to inspect the identifications.
    
    # Initialize the printing variables.
    count = 0  # Counts the number of partitions being checked.
    countCS = 0  # Counts how many made it to be yielded as core subgraphs.
    print "#Partitions      #CS     Time (Cum.)"
    print "..........0        0     begin!"
    begin = time.clock()
    
    # restricted_partition_generator generates all identifications of the vertices (partitions of the vertex set).  Its output is a partition, stored as a list of lists of vertices.
    for partition in restricted_partition_generator(range(NCS.order),NCS.order-1,forbidden_dict):
        count += 1
        if count % 10000 == 0:  # Print every so often -- more helpful for large configurations.
            print "."*max(11-len(str(count)),0)+str(count)+" "*max(9-len(str(countCS)),0)+str(countCS)+"     "+timestring(time.clock()-begin)
        
        # We will want to move between the current labelings of the initial plane subgraph and the new graph formed by the identification.  ident_dict will help with this by storing foo:bar where foo is a vertex from the original graph and bar is a vertex in the identified graph.  If vertices v1 and v2 are identified, then ident_dict[v1] == ident_dict[v2].  The labels of the vertices in the identified graph are the indices of the parts in partition.
        ident_dict = {y:x for x in range(len(partition)) for y in partition[x]}
        
        # To start checking the identified graph, we make the set of edges in the identified graph.
        new_edges = set()
        # For each edge in the old graph between v1 and v2, we add an edge between the parts y1 and y2 (where y1 is the part containing v1 and y2 is the part containing v2).  Since we are using a set, any duplicate edges get thrown out.
        for e in NCS.edges:
            a = ident_dict[e[0]]
            b = ident_dict[e[1]]
            if a > b:  # We enforce an ordering so there are no duplicates via reversal.  Given a partition, an edge (X,Y) in the identified graph is represented by an edge (x1,y1) where x1 in X and y1 in Y.  However, it might also be represented by an edge (y2,x2) where x2 in X and y2 in Y.  We don't want to represent an edge more than once in our set new_edges.
                a,b = b,a
            new_edges.add((a,b))
            # Note:  Our dictionary of forbidden identifications ensures that a != b.
        
        some_condition_failed = False  # If some_condition_failed switches to True, then we stop checking this partition and continue to the next partition.
        
        # Now that we have our edges, our first condition to check is if the identified graph is subcubic.  (Observation ?.1)
        # For each vertex in the identified graph, we count how many edges it is in.
        for x in range(len(partition)):
            degree_of_x = 0
            for e in new_edges:
                if x in e:
                    degree_of_x += 1
            if degree_of_x > 3:  # If we find a new vertex of degree more than 3, break immediately and continue to the next partition.
                some_condition_failed = True
                break
        
        if some_condition_failed:
            continue
        
        # To check the next conditions, we use the actual Graph object and get information about its specified faces.
        
        # First, get the graph.
        J = Graph(list(new_edges))
        # The graph we will check (J_prime) might have an additional edge to close the central face.
        J_prime = Graph(J)
        # We add the edge if the configuration is open and the spine ends were not identified.
        new_end_1,new_end_2 = ident_dict[NCS.faces[0][0]],ident_dict[NCS.faces[0][-1]]
        if NCS.is_open and new_end_1 != new_end_2:
            J_prime.add_edge([new_end_1,new_end_2])  # If the graph already has this edge, that is okay.  Then the edge is already serving the desired purpose.
        
        # Then, get the list of central/specified faces for this new graph.  We will need to check for possible identification of faces in the identification process, as we want the faces to be distinct when we check the next two conditions of Observation ? and get the open regions of the graph.
        # Note:  We are not checking that the faces have been identified with the same orientations of their boundaries (clockwise vs. counterclockwise).  Although this happens to be a necessary condition for a core subgraph, it is enforced by our planarity check later.
        new_faces = []  # Initialize new faces -- it will be the faces of J_prime.
        new_faces.append([ident_dict[x] for x in NCS.faces[0]])  # Start with central face.
        central_face = new_faces[0]  # It will help readability to refer to 'central_face' instead of 'new_faces[0]'.
        # If the central face length is unspecified and the spine ends were identified, then we need to remove one of the spine ends to avoid duplicate vertices in the face list.
        if NCS.is_open and new_end_1 == new_end_2:
            central_face.pop()
        # Now we deal with the rest of the faces.
        for face in NCS.faces[1:]:
            new_face = [ident_dict[x] for x in face]  # Convert to new vertex labels.
            for prev_face in new_faces[1:]:  # We start comparison at index 1 to avoid comparing with the central face.  No specified face can be identified with the central face because the resulting face would contain a cut edge (the shared edge between the two faces in the natural core subgraph).  Since the face must also have a length at most 9, this violates Lemma ?.
                if set(new_face) == set(prev_face):  # If the set of vertices of the current new face is the same set of vertices as a previous face, then these faces must have been identified (because the graph is not a cycle).
                    break
            else:
                new_faces.append(new_face)
        
        # Now we can check to see if the identified graph has any trapped 2-vertices.  (Observation ?.3 -- we check this before Observation ?.2 because it is easier)
        for v in J_prime.vertices():
            if J_prime.degree(v) != 2:
                continue
            # The 2-vertex v is "trapped" if it appears more than once on the boundary walks of the fixed faces of J_prime.
            boundary_count_for_v = 0
            for face in new_faces:
                boundary_count_for_v += face.count(v)
            if boundary_count_for_v > 1:
                some_condition_failed = True
                break
        
        if some_condition_failed: 
            continue
        
        # Now we begin checking Observation ?.2.  The first thing we check is that the spine does not wrap on itself -- Observation ?.2 says that the spinal path should be a boundary walk of a face in J_prime.  A boundary walk does not repeat an edge in the same direction.
        # (We need to check this only when the spine was open to begin with.)
        if NCS.is_open:
            for i,j in combinations(range(len(central_face)-1),2):
                if central_face[i] == central_face[j] and central_face[i+1] == central_face[j+1]:
                    some_condition_failed = True
                    break
            if some_condition_failed:
                continue
        
        # Now we finish checking Observation ?.2.  For each fixed face of J_prime, we add a new hub vertex to be a representative for the face.  Then for each edge on the boundary of the face, another vertex is added which is connected to both the endpoints of the edge and to the midpoint of the edge and to the hub vertex representing the face.  In other words, we glue a wheel-like structure into each face in a particular way along the face boundary.  Then, we will check the resulting graph for planarity.
        
        i = J_prime.order()  # i will index the label of each new vertex.
        
        # First, we subdivide each edge so that we can connect some of our new vertices to the midpoints of our original edges.
        edge_dict = {}  # edge_dict will keep track of the labels of the midpoints of our original edges.
        for edge in J_prime.edges():  # The Graph edges iterator does not add edges part-way through the for loop.  Even though the edges of J_prime are changing over the course of the for loop, the objects in the iterator remain only those originally in J_prime when the for loop first started.
            J_prime.subdivide_edge((edge[0],edge[1]),1)
            edge_dict[(edge[0],edge[1])] = i
            edge_dict[(edge[1],edge[0])] = i  # Keep track of both orientations of the edge -- we do not know which orientations we will need.
            i += 1
        
        for face in new_faces:
            face_hub_vertex = i  # face_hub_vertex is the hub of the wheel-like structure being embedded on the face.
            i += 1
            J_prime.add_vertex(face_hub_vertex)
            for j in range(len(face)):
                edge_hub_vertex = i  # For each edge on the face, edge_hub_vertex is the vertex which is adjacent to face_hub_vertex and also to both endpoints and midpoint of the edge.
                i += 1
                J_prime.add_vertex(edge_hub_vertex)
                J_prime.add_edge((edge_hub_vertex,face[j-1]))  # First endpoint.  Still works when j==0 due to Python's indexing.
                J_prime.add_edge((edge_hub_vertex,face[j]))  # Second endpoint.
                J_prime.add_edge((edge_hub_vertex,edge_dict[(face[j-1],face[j])]))  # Midpoint.
                J_prime.add_edge((edge_hub_vertex,face_hub_vertex))
        
        
        # Finally, check for planarity.
        if J_prime.is_planar():
            # If we made it this far, then all conditions were passed!
            countCS += 1
            # Print info and yield J (not J_prime).
            print
            print "      >>> "+"Core Subgraph #%d"%(countCS),"(Partition #%d)"%(count)
            print " "*10+"Partition: %s"%(str(partition))  
            print " "*10+"Edges: %s"%(str([e[:2] for e in J.edges()]))
            yield J,new_faces
    
    print "Total #partitions:   ",count











#------------------------------------------
#  Section 3:  Stem Structure Generation
#------------------------------------------
    






def get_open_region_roots(graph,given_faces):
    # graph is just the Graph object from the core subgraph generator.
    # given_faces is a list of some faces in graph (the central face/spine as well as any other specified faces).  Some faces might have been identified by core_subgraph_generator -- in this case, only one representative remains in given_faces.  The cyclic ordering for each face is in the same clockwise orientation.
    # given_faces[0] is the central face/spine.
    # We turn each edge of the graph into two oppositely directed edges and then remove any of the edges going in the direction of given_faces.  The remaining directed edges are then partitioned into disjoint directed cycles, each of which is an open region of the graph.
    # We return a list of open regions, each given as a list of roots in that region in clockwise order.
    
    
    # First, we consider all directed edges in the graph.
    all_directed_edges = {(e[0],e[1]) for e in graph.edges()} | {(e[1],e[0]) for e in graph.edges()}
    # (e[0],e[1]) is one direction, (e[1],e[0]) is the other direction.  (These directions are based on how edges are stored in graph, not on the orientation in given_faces.)
    
    # In addition to the edges in graph, we also want to close up the spine if the central face is open.  This is because the spine ends are connected by some closed curve, separating the central face from other open regions.
    # If the central face was specified to begin with or was closed by core_subgraph_generator, then this next line is redundant.
    all_directed_edges |= {(given_faces[0][0],given_faces[0][-1]),(given_faces[0][-1],given_faces[0][0])}
    
    # Second, we remove all directed edges from given_faces (including across the spine ends).
    facial_directed_edges = {(face[j-1],face[j]) for face in given_faces for j in range(len(face))}
    # pointer_dict is a dictionary which assigns each tail of a remaining directed edge to the head of the directed edge.
    pointer_dict = {e[0]:e[1] for e in all_directed_edges-facial_directed_edges}
    # The directed graph induced by all_directed_edges-facial_directed_edges is a collection of disjoint directed cycles because the remaining directed edges are precisely the boundaries of the faces not in given_faces and because each edge of the graph is on the boundary of at least one face from given_faces.  Hence pointer_dict assigns a value to each key precisely once.
    outer_regions = []  # outer_regions will be our list of open regions.
    
    while pointer_dict:
        v = pointer_dict.keys()[0]  # v is an arbitrary initial vertex.
        new_region = [v]  # new_region will eventually be a directed cycle.  Initialize it as the list containing just the element v.
        vtx = pointer_dict[v]  # vtx is the next vertex in the directed cycle.
        del pointer_dict[v]  # Now that we know vtx, delete the directed edge we just looked at.
        # At the top of the loop, we know what the next vertex is, but we have not yet added it to new_region.
        while vtx != v:  # Until our next vertex is our initial vertex.
            new_region.append(vtx)  # If the next vertex is not the initial vertex, add it to the list.
            next_vtx = pointer_dict[vtx]  # Now step to the next vertex, and look at the new next vertex past it.
            del pointer_dict[vtx]  # Delete the directed edge we just looked at.
            vtx = next_vtx 
        outer_regions.append(new_region)
    
    # We now have the all the vertices on the open regions, but we care only about the roots.
    roots = {v for v in graph.vertices() if graph.degree(v)==2}  # This includes the open spine ends, if there are any.
    return [[z for z in region if z in roots] for region in outer_regions]












def stem_structure_generator(roots_by_region):
    # roots_by_region is a list of lists, each holding a cyclic ordering of distinct roots in a region.  We generate all possible ways of adding new edges or new stems of degree 1, 2, or 3 that preserve planarity within each region.
    # We loop through the possibilities for identifying the last root in the last region with earlier roots in the last region.  For each possibility, we fragment the last region into the resulting smaller regions and recursively call the generator.
    # The yielded output ("structure") is a list of three lists:  the pairs (from all the regions) which form a new edge, the pairs which solely share a stem (2-stem), and the triples which share a stem (3-stem).
    
    # Note:  We do not bother to check if a "new" edge already exists in the core subgraph (for example, if two roots are adjacent on a face, then there is already an edge between them).  If the edge is already present, then when we add it in core_square_graph there is no effect.  So we are generating redundant realizations in exchange for simplicity.  (In fact, we generate redundant realizations as a result of certain 2-stems and 3-stems as well.)
    
    # If there are no regions with any available roots left, then yield the empty structure.
    if len(roots_by_region) == 0:
        yield [[],[],[]]
        return
    
    # Otherwise, examine the last region.
    else:
        roots = roots_by_region[-1]
        
        # If the last region has at most one root, then there are no possible pairs or triples in that last region.  This means that the structures we can obtain from roots_by region are exactly the structures we can obtain from roots_by_region[:-1].
        if len(roots) < 2:
            for substructure in stem_structure_generator(roots_by_region[:-1]):
                yield substructure
            return
        
        # If the last region has at least two roots, then we can try to form pairs or triples.  We do so conditioned on the last root.
        else:
            
            # Case 1:  Condition on the last root being in a 3-stem.
            for two in combinations(range(len(roots)-1),2):
                # Break up the last region into three regions.
                new_roots_by_region = roots_by_region[:-1]
                new_roots_by_region.append(roots[:two[0]])
                new_roots_by_region.append(roots[two[0]+1:two[1]])
                new_roots_by_region.append(roots[two[1]+1:-1])
                # For each substructure on new_roots_by_region, append our triple to the list of 3-stems and yield.
                for substructure in stem_structure_generator(new_roots_by_region):
                    structure = [[set(x) for x in y] for y in substructure]
                    structure[2].append({roots[two[0]],roots[two[1]],roots[-1]})
                    yield structure
            
            # Case 2:  Condition on the last root being in a new edge.
            for one in range(len(roots)-1):
                # Break up the last region into two regions.
                new_roots_by_region = roots_by_region[:-1]
                new_roots_by_region.append(roots[:one])
                new_roots_by_region.append(roots[one+1:-1])
                # For each substructure on new_roots_by_region, append our pair to the list of edges and yield.
                for substructure in stem_structure_generator(new_roots_by_region):
                    structure = [[set(x) for x in y] for y in substructure]
                    structure[0].append({roots[one],roots[-1]})
                    yield structure
            
            # Case 3:  Condition on the last root being in a 2-stem.  (Virtually identical to Case 2.)
            for one in range(len(roots)-1):
                # Break up the last region into two regions.
                new_roots_by_region = roots_by_region[:-1]
                new_roots_by_region.append(roots[:one])
                new_roots_by_region.append(roots[one+1:-1])
                # For each substructure on new_roots_by_region, append our pair to the list of 2-stems and yield.
                for substructure in stem_structure_generator(new_roots_by_region):
                    structure = [[set(x) for x in y] for y in substructure]
                    structure[1].append({roots[one],roots[-1]})
                    yield structure
            
            # Case 4:  Condition on the last root being in a 1-stem (by itself).
            # Remove the last root from the last region.
            new_roots_by_region = roots_by_region[:-1]
            new_roots_by_region.append(roots[:-1])
            # In this case, the structures we can obtain from roots_by region are exactly the structures we can obtain from new_roots_by_region.
            for substructure in stem_structure_generator(new_roots_by_region):
                yield substructure
            
            # After considering all the possibilities for the last root in the last region, we are finished.
            return

























#------------------------------------------
#  Section 4:  Checking Realizations
#------------------------------------------




def check_all_stem_structures_for_given_core_subgraph(graph,open_region_roots):
    # Given a core subgraph (as a graph and a list of root orderings in the open regions), we generate all stem structures among the roots in their respective open regions.  For each resulting realization, we check core-choosability.
    
    # Initialize summary variables
    begin = time.clock()
    count = 0  # Total number of stem structures (and thus realizations) for this core subgraph
    good_count = 0  # Number of core-choosable realizations
    bad_count = 0  # Number of non-core-choosable realizations
    
    print " "*15+"Index      Good   Bad   Time (Cum.)   Stem Structure"
    
    # Now we check them all.
    i = 0  # i is the index of each stem structure, for printing purposes
    for SS in stem_structure_generator(open_region_roots):
        i += 1
        G_sq,f = core_square_graph(graph,SS)
        
        # fChoosableNoPrint takes a graph and a list size function and returns a tuple of two elements.  The first element is boolean:  True means that the graph is f-choosable, and False means it is not.  The second element is a string indicating the method used for its evaluation:  'error', 'greedy', 'CNS', or 'brute'
        y = ch.fChoosableNoPrint(G_sq,f,inprocess=4,print_mod=100,maxIS=True,rev=False)
        if y[0]:
            good_count += 1
        else:
            bad_count += 1
        ti = timestring(time.clock()-begin)
        print " "*(20-len(str(i)))+str(i)+" "*(10-len(str(good_count)))+str(good_count)+" "*(6-len(str(bad_count)))+str(bad_count)+" "*(14-len(ti))+ti+"   "+str(SS)
    
    # Print a concluding statement.
    if bad_count:
        print "      >>> Uh-oh!  Problem!"
    else:
        print "      >>> This core subgraph is good with all stem structures!"
    print
    
    # Return the total number of realizations and the number of bad realizations.
    return good_count,bad_count














def check_all_realizations_from_initial_plane_graph(NCS,forbid_identifications={}):
    # NCS is an object from the NaturalCoreSubgraph class which holds all the information we need to pass on to core_subgraph_generator.  For Lemma 20, we do not need to change anything about NCS from what the contructor did.  However, for Lemma 29 we will have already made some changes manually.
    # forbid_identifications is a dictionary holding vertex:{set of vertices which cannot be identified with vertex}.  In addition to the forbidden identifications from Lemma 9 that we use in core_subgraph_generator, this allows us to manually insert other restrictions for Lemma 29.
    # We use core_subgraph_generator to generate all core subgraphs which can be formed from our initial plane graph.  For each of these, we find the open regions and list the roots in each region in cyclic order.  Using this information, for each core subgraph we generate all possible stem structures and check the resulting realizations for core-choosability.  The results are printed.
    
    # Initialize summary variables
    begin = time.clock()
    core_subgraph_count = 0  # Counts the total number of core subgraphs
    realization_count = 0  # Counts the total number of realizations
    bad_count = 0  # Counts the total number of realizations which are not core-choosable
    
    # totals_str will keep track of the individual totals for each core subgraph
    totals_str = "      CS#               #SS       #Bad\n"
    
    # For each core subgraph, we generate each stem structure and check the resulting realization for core-choosability.
    for CS in core_subgraph_generator(NCS,forbid_identifications):
        # CS is a tuple of (a graph, a list of faces in the graph).
        
        # Get the necessary information for the stem structure generator.
        # open_region_roots is a list containing orderings of roots from each open region from the core subgraph.
        open_region_roots = get_open_region_roots(CS[0],CS[1])
        
        # Check the stem structures
        good_stem_structure_count,bad_stem_structure_count = check_all_stem_structures_for_given_core_subgraph(CS[0],open_region_roots)
        
        # Update the summary info
        core_subgraph_count += 1
        realization_count += good_stem_structure_count+bad_stem_structure_count
        bad_count += bad_stem_structure_count
        core_subgraph_count_str = str(core_subgraph_count)
        stem_structure_count_str = str(good_stem_structure_count+bad_stem_structure_count)
        bad_stem_structure_count_str = str(bad_stem_structure_count)
        totals_str += " "*(max(9-len(core_subgraph_count_str),0))+core_subgraph_count_str+" "*(max(18-len(stem_structure_count_str),0))+stem_structure_count_str+" "*(max(11-len(bad_stem_structure_count_str),0))+bad_stem_structure_count_str+"\n"
    
    # Wrap up and print all that information!
    end = time.clock()
    
    # First, print the summary for each core subgraph.
    print "Done!\nSummary:"
    print totals_str
    
    # Second, print the totals from all core subgraphs.
    print "Totals:"
    print "      #CS     #Realizations       #Bad"
    core_subgraph_count_str = str(core_subgraph_count)
    realization_count_str = str(realization_count)
    bad_count_str = str(bad_count)
    print " "*(max(9-len(core_subgraph_count_str),0))+core_subgraph_count_str+" "*(max(18-len(realization_count_str),0))+realization_count_str+" "*(max(11-len(bad_count_str),0))+bad_count_str
    
    # Finally, print a concluding statement and total runtime
    if bad_count > 0:
        print "\nNope!  Some realizations are not core-choosable."
    else:
        print "\nGood!  All realizations are core-choosable!"
    print "Time: "+timestring(end-begin)

