#Some things to-do (before going through line-by-line):
#    Adjust printing output to have readable final format.  (also correct "instances" vs "realizations")
#    Make sure that check_c7a4x5x5 stuff is actually checking the right restrictions --- if there is more than one FAS, do the stem restrictions get passed on correctly?  (whatever that means)
#    Revise Sections 4 and 5.
#    Check imported modules.


#------------------------------------------------------
#CONTENTS
#------------------------------------------------------


#SECTION 0:  Libraries and Comments


#SECTION 1:  Graph Essentials
#    graph_power
#    makeGraph
#    core_square_graph


#SECTION 2:  Checking Core-choosability of Subgraph Realizations
#    NS_generator_no_enforced_planarity
#    check_all_subgraph_realizations


#SECTION 3:  Checking Core-choosability of Configuration Realizations
#
#  Subsection A:  Facial Adjacency Structure Generation
#    restricted_partitions
#    nonidentified_faces
#    configuration_FAS_generator
#
#  Subsection B:  Neighborhood Structure Generation
#    get_outer_region_roots
#    NS_generator_with_enforced_planarity
#
#  Subsection C:  Checking All Realizations
#    check_all_neighborhood_structures_for_FAS
#    check_all_realizations_from_initial_plane_graph
#    check_all_configuration_realizations
#    check_all_realizations_from_expanded_c7a4x5x5_case


#SECTION 4:  The Target Set
#    HRCs
#    UsedInHappy3
#    UsedInHappy4
#    UsedInHappy5
#    UsedInDischarging4
#    UsedInDischarging5
#    initializeRedConList


#SECTION 5:  Automated Discharging
#    nonreducible_block_generator
#    reducible_check_last_face_large_face
#    str_converter
#    descriptor_match
#    nonreducible_block_stack_generator
#    reducible_check_last_face
#    stack2str
#    final_charge
#    charge_pull
#    simmer_once
#    simmer
#    less_than
#    fnmatch_gen
#    run_discharging_analysis
#    run_discharging_analysis_no_print




















#------------------------------------------------------
#SECTION 0:  Libraries and Comments
#------------------------------------------------------


import copy
from sage.graphs.graph import Graph
from fractions import Fraction
from sage.rings.rational import Rational

from itertools import combinations
#Used in:
#    NS_generator_no_enforced_planarity
#    NS_generator_with_enforced_planarity
#another?










import choosability as ch
import wegner as wg
import time
import copy



import fnmatch



























#------------------------------------------------------
#SECTION 1:  Graph Essentials
#------------------------------------------------------


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




















#Not made for cxa33 or cxa34 (but perhaps made for c3a3 and c3a4?)
#Not made for full chain around a specified face.  (c4a34x5 should be c4a534.)
def makeGraph(configuration):
    a = configuration.index('a')
    
    if configuration[1] == 'x':
        b = 'x'
    else:
        b = int(configuration[1:a])
    
    chain = []
    for s in configuration[a+1:]:
        if s == 'x':
            chain.append('x')
        else:
            chain.append(int(s))
    
    E = []#E is our set of edges.
    faces=[]#Each entry is a clockwise ordering of the boundary of a specified face.
    #(Except for the first one.  If cxa-type, then the first "face" is just the spine.)
    roots = []#2-vertices
    
    if b == 'x':
        for j in range(len(chain)):
            E.append((j,j+1))
        spine = range(len(chain)+1)
        faces.append(spine[:])
        vtx = spine[-1]
        #vtx is the current vertex we're building edges to and from.  We'll increment it before using it.
    else:
        for j in range(b-1):
            E.append((j,j+1))
        E.append((0,b-1))
        spine = range(len(chain)+1)
        faces.append(range(b))
        vtx = b-1
    
    for i in range(len(chain)):
        if chain[i]=='x':
            if chain[i-1]=='x':
                roots.append(i)
            #Not the first entry, so there is potentially a previous face to tie off.
            else:
                E.append((i,vtx))
                if chain[i-1]==3 and i>1 and chain[i-2]!='x':
                    del roots[-1]
            continue
        
        #Otherwise, set up first edge and then middle outer edges.
        if i==0 or chain[i-1]=='x':
            vtx += 1
            roots.append(vtx)
        else:#i>0 and previous/current two faces are not 'x'
            del roots[-1]#Then vtx is a 3-vertex, not a root

        #First edge.
        E.append((i,vtx))
        
        #Special cases with a previous 3-face:
        if i>0 and chain[i-1]==3:
            #del roots[-1]
            if i==1 or chain[i-2]=='x':
                pass
            else:
                vtx += 1
                E.append((vtx-2,vtx))
                faces.append([i+1,i,vtx-1,vtx-2,vtx])
                del roots[-1]#Remove vtx-2 from roots
                roots.append(vtx)
                for j in range(chain[i]-5):
                    vtx += 1
                    faces[-1].append(vtx)
                    E.append((vtx-1,vtx))
                    roots.append(vtx)
                continue
        
        #Normal case:
        faces.append([i+1,i,vtx])
        for j in range(chain[i]-3):
            vtx += 1
            faces[-1].append(vtx)
            E.append((vtx-1,vtx))
            roots.append(vtx)
    
    E.append((len(spine)-1,vtx))
    
    if b=='x':
        roots.insert(0,spine[0])
        roots.append(spine[-1])
    else:
        for j in range(len(chain)+1,b):
            roots.append(j)
    
    return vtx+1,E,spine,roots,faces












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


























































#------------------------------------------------------
#SECTION 2:  Checking Core-choosability of Subgraph Realizations
#------------------------------------------------------













#roots is a list of distinct objects.
#restrictions is a list of three dictionaries of the form {v:S}:
#    In restrictions[0], S is the set of vertices that v cannot make a new edge with.
#    In restrictions[1], S is the set of vertices that v cannot make a 2-stem with.
#    In restrictions[2], S is the set of 2-sets of vertices that v cannot make a 3-stem with.
def NS_generator_no_enforced_planarity(roots,restrictions=False):
    if not restrictions:
        edge_restrictions = {x:set([]) for x in roots}
        ident_2_restrictions = {x:set([]) for x in roots}
        ident_3_restrictions = {x:set([]) for x in roots}
        restrictions = [edge_restrictions,ident_2_restrictions,ident_3_restrictions]
    
    if len(roots) < 2:
        yield [[],[],[]]
        return
    else:
        #We recursively proceed by cases based on the last root.
        
        #Case 1:  Condition on the last root being in a three-stem.
        for two in combinations(range(len(roots)-1),2):
                
            if frozenset([roots[x] for x in two]) in restrictions[2][roots[-1]]:
                continue
            
            #close_count = 0
            #close_couples = set(frozenset(x) for x in dist1) | set(frozenset(x) for x in dist2)
            #triple_couples = set([frozenset([roots[two[0]],roots[two[1]]])])
            #triple_couples |= set([frozenset([roots[two[0]],roots[-1]])])
            #triple_couples |= set([frozenset([roots[two[1]],roots[-1]])])
            #if len(triple_couples & close_couples) > 2:
                #continue
                
            new_roots = roots[:-1]
            del new_roots[two[1]]
            del new_roots[two[0]]
            for substructure in NS_generator_no_enforced_planarity(new_roots,restrictions=restrictions):
                structure = [x[:] for x in substructure]
                structure[2].append([roots[two[0]],roots[two[1]],roots[-1]])
                yield structure
        
        #Case 2:  Condition on the last root being in a new edge.
        for one in [x for x in range(len(roots)-1) if not roots[x] in restrictions[0][roots[-1]]]:
        #too_close = {x for y in dist1 if roots[-1] in y for x in y}
        #for one in [x for x in range(len(roots)-1) if not roots[x] in too_close]:
            new_roots = roots[:-1]
            del new_roots[one]
            for substructure in NS_generator_no_enforced_planarity(new_roots,restrictions=restrictions):
                structure = [x[:] for x in substructure]
                structure[0].append([roots[one],roots[-1]])
                yield structure
        
        #Case 3:  Condition on the last root being in a two-stem.
        for one in [x for x in range(len(roots)-1) if not roots[x] in restrictions[1][roots[-1]]]:
        #too_close |= {x for y in dist2 if roots[-1] in y for x in y}
        #for one in [x for x in range(len(roots)-1) if not roots[x] in too_close]:
            new_roots = roots[:-1]
            del new_roots[one]
            for substructure in NS_generator_no_enforced_planarity(new_roots,restrictions=restrictions):
                structure = [x[:] for x in substructure]
                structure[1].append([roots[one],roots[-1]])
                yield structure
            
        
        #Case 4:  Condition on the last root being in a one-stem.
        new_roots = roots[:-1]
        for substructure in NS_generator_no_enforced_planarity(new_roots,restrictions=restrictions):
            structure = [x[:] for x in substructure]
            yield structure
            
        #Done.
        return











def check_all_subgraph_realizations(subgraph):
    begin = time.clock()
    roots = [v for v in subgraph.vertices() if subgraph.degree(v)==2]
    dist1 = []
    dist2 = []
    for pair in combinations(roots,2):
        d = subgraph.distance(pair[0],pair[1])
        if d == 1:
            dist1.append(set(pair))
        elif d == 2:
            dist2.append(set(pair))
    count_NS = 0
    count_r = 0
    count_good = 0
    count_bad = 0
    
    
    
    edge_restrictions = {x:set([]) for x in subgraph.vertices()}
    ident_2_restrictions =  {x:set([]) for x in subgraph.vertices()}
    ident_3_restrictions =  {x:set([]) for x in subgraph.vertices()}
    
    close_couples = set([])
    for pair in dist1:
        for v in pair:
            edge_restrictions[v] |= pair
            ident_2_restrictions[v] |= pair
        close_couples |= set([frozenset(pair),])
    for pair in dist2:
        for v in pair:
            ident_2_restrictions[v] |= pair
        close_couples |= set([frozenset(pair),])
    for triple in combinations(roots,3):
        triple_couples = set([frozenset([triple[0],triple[1]])])
        triple_couples |= set([frozenset([triple[0],triple[2]])])
        triple_couples |= set([frozenset([triple[1],triple[2]])])
        if len(triple_couples & close_couples) > 2:
            ident_3_restrictions[triple[0]] |= set([frozenset(triple[1:]),])
            ident_3_restrictions[triple[1]] |= set([frozenset([triple[0],triple[-1]]),])
            ident_3_restrictions[triple[2]] |= set([frozenset(triple[:2]),])
    
    for structure in NS_generator_no_enforced_planarity(roots,restrictions=[edge_restrictions,ident_2_restrictions,ident_3_restrictions]):
        count_NS += 1
        planar_check = Graph(subgraph.edges())
        for e in structure[0]:
            planar_check.add_edge(e)
        for pair in structure[1]:
            planar_check.add_edge(pair)
        i = subgraph.order()
        for triple in structure[2]:
            planar_check.add_edge([i,triple[0]])
            planar_check.add_edge([i,triple[1]])
            planar_check.add_edge([i,triple[2]])
            i += 1
        if not planar_check.is_planar():
            continue
        else:
            count_r += 1
            core_square,f = core_square_graph(subgraph,structure)
            y = ch.fChoosableNoPrint(core_square,f,inprocess=4,print_mod=100,maxIS=True,rev=False)
            if y[0]:
                count_good += 1
            else:
                count_bad += 1
    print "Done!  Time:",ch.timestring(time.clock()-begin)
    print "#neighborhood structures attempted:",count_NS
    print "#realizations:",count_r
    print "#good:",count_good
    print "#bad:",count_bad
    if count_bad:
        print "Since not all realizations are core-choosable, the subgraph is not necessarily reducible."
    else:
        print "Since all realization are core-choosable, the subgraph is reducible!"



































































#------------------------------------------------------
#SECTION 3:  Checking Core-choosability of Configuration Realizations
#------------------------------------------------------




#--------------------------
#  Subsection A:  Facial Adjacency Structure Generation
#--------------------------




def restricted_partitions(li,di):
    if len(li) == 0:
        yield []
        return
    if len(li) == 1:
        yield [li[:],]
        return
    else:
        for partition in restricted_partitions(li[:-1],di):
            for i in [x for x in range(len(partition)) if not set(partition[x]) & di[li[-1]]]:
                new_partition = [x[:] for x in partition]
                new_partition[i].append(li[-1])
                yield new_partition
            new_partition = partition[:]
            new_partition.append([li[-1],])
            yield new_partition
        return









def nonidentified_faces(faces,id_dict):
    new_faces = []
    for face in faces[1:]:#Skip the central face; insert later.
        #(Central face can't be identified with outer faces since they share an edge in opposite directions.)
        flag = True
        #flag is whether this face will be represented with its own vertex
        #flag = False means we don't add it to new_faces.
        for prev_face in new_faces:
            if len(prev_face) != len(face):
                continue
            prev_parts = [id_dict[x] for x in prev_face]
            face_parts = [id_dict[x] for x in face]
            li = prev_parts[:]
            li.extend(prev_parts[:-1])
            for i in range(len(prev_parts)):
                if face_parts == li[i:i+len(prev_parts)]:
                    flag = False
                    break
                #We need to not check reverse because faces must have same orientation.
        if flag:
            new_faces.append(face[:])
    new_faces.insert(0,faces[0][:])#central face cannot be identified with specified outer faces for configurations in the target set.
    return new_faces













def configuration_FAS_generator(g,open_spine=False,include_restrictions={}):
    if open_spine:
        spine_ends = [g[4][0][0],g[4][0][-1]]
    else:
        spine_ends = []
    di = {x:set([]) for x in range(g[0])}
    for x in include_restrictions.keys():
        di[x] |= set(include_restrictions[x])
    if open_spine:
        init_spec_face = 1
        the_spine = g[4][0]
        di[the_spine[0]] |= set(the_spine [0:3])
        di[the_spine[1]] |= set(the_spine[0:4])
        di[the_spine[-2]] |= set(the_spine[-4:])
        di[the_spine[-1]] |= set(the_spine[-3:])
        for v in the_spine[2:-2]:
            di[v] |= set(the_spine[v-2:v+3])
    else:
        init_spec_face = 0
    for face in g[4][init_spec_face:]:
        for v in face:
            di[v] |= set(face)
    count = 0
    countc = 0
    countt = 0
    countp = 0
    print "#Partitions     #FAS     Time (Cum.)"
    print "..........0        0     begin!"
    begin = time.clock()
    for partition in restricted_partitions(range(g[0]),di):
        count += 1
        id_dict = {y:x for x in range(len(partition)) for y in partition[x]}
        new_edges = set([])
        for e in g[1]:
            a = id_dict[e[0]]
            b = id_dict[e[1]]
            if a > b:
                a,b = b,a
            new_edges.add((a,b))
            #If a==b, don't add anything.
        flag = True
        for x in range(len(partition)):
            if len([0 for e in new_edges if x in e])>3:
                flag = False
                break
        if flag:
            countc += 1
            new_G = Graph(list(new_edges))
            nonid_faces = nonidentified_faces(g[4],id_dict)
            nonid_faces = [[id_dict[x] for x in f] for f in nonid_faces]
            if nonid_faces[0][0]==nonid_faces[0][-1]:
                del nonid_faces[0][-1]
            #Check for trapped 2-vertices
            for v in new_G.vertices():
                if new_G.degree(v) != 2:
                    continue
                #The following condition for spinal ends also holds even if the two are identified.
                if open_spine and v in (id_dict[spine_ends[0]],id_dict[spine_ends[1]]):
                    if len([0 for f in nonid_faces[1:] if v in f]) > 1:
                        flag = False
                        break
                else:
                    if len([0 for f in nonid_faces if v in f]) > 1:
                        flag = False
                        break
            if flag:      
                countt += 1
                n = new_G.order()
                i = n
                
                newer_G = Graph(new_G)
                if open_spine:
                    a,b = id_dict[spine_ends[0]],id_dict[spine_ends[1]]
                    if a > b:
                        a,b = b,a
                    if a!=b and not (a,b) in new_edges:
                        newer_G.add_edge((a,b))
                initial_edges = newer_G.edges()[:]
                edge_di = {}
                for e in initial_edges:
                    edge_di[e[0:2]] = i
                    edge_di[(e[1],e[0])] = i
                    newer_G.subdivide_edge(e[0:2],1)
                    i += 1
                for face in nonid_faces:
                    for j in range(len(face)):
                        newer_G.add_edge((i,edge_di[(face[j-1],face[j])]))
                    i += 1
                    
                if newer_G.is_planar():
                    countp += 1
                    print
                    print "      >>> "+"Realization #%d"%(countp),"(Partition #%d)"%(count)
                    print " "*10+"Partition: %s"%(str(partition))
                    ed = [e[:2] for e in new_G.edges()]
                    ed_str_li = []
                    rows = len(ed)/10
                    for x in range(rows):
                        ed_str_li.append(str(ed[10*x:10*(x+1)])[1:-1]+",")
                    ed_str_li.append(str(ed[10*rows:])[1:-1])
                    print " "*10+"Edges: %s"%(ed_str_li[0])
                    for s in ed_str_li[1:]:
                        print " "*17+s
                    yield new_G,nonid_faces,partition,id_dict
        if count % 10000 == 0:
            print "."*max(11-len(str(count)),0)+str(count)+" "*max(9-len(str(countp)),0)+str(countp)+"     "+ch.timestring(time.clock()-begin)
    end = time.clock()
    t = end-begin
    #print "Total partitions:   ",count
    #print "Passed cubic test:  ",countc
    #print "Passed trapped test:",countt
    #print "Passed planar test: ",countp
    #print "Time:",t











#--------------------------
#  Subsection B:  Neighborhood Structure Generation
#--------------------------
    















#The fact that every vertex is on at most one outer region also implies that each edge is on at most one outer region.  So we can list the directed edges which don't show up in the specified faces and follow them around the cycle boundaries of the outer regions.
def get_outer_region_roots(graph,specified_faces):
    #graph.show()
    #print specified_faces
    R = {v for v in graph.vertices() if graph.degree(v)==2}
    all_diedges = {(e[0],e[1]) for e in graph.edges()} | {(e[1],e[0]) for e in graph.edges()}
    #Assume spine ends are connected by some closed curve (acts like a diedge).
    all_diedges |= {(specified_faces[0][0],specified_faces[0][-1]),(specified_faces[0][-1],specified_faces[0][0])}
    facial_diedges = {(face[j-1],face[j]) for face in specified_faces for j in range(len(face))}
    pointer_di = {e[0]:e[1] for e in all_diedges-facial_diedges}
    #print pointer_di
    outer_regions = []
    while pointer_di:
        v = pointer_di.keys()[0]
        new_region = [v,]
        vtx = pointer_di[v]
        del pointer_di[v]
        #print v
        while vtx != v:
            #print vtx
            new_region.append(vtx)
            next_vtx = pointer_di[vtx]
            del pointer_di[vtx]
            vtx = next_vtx
        outer_regions.append(new_region[:])
    return [[z for z in x if z in R] for x in outer_regions]












#roots_by_region is a list of lists, each of distinct objects.
#restrictions is a list of three dictionaries of the form {v:S}:
#    In restrictions[0], S is the set of vertices that v cannot make a new edge with.
#    In restrictions[1], S is the set of vertices that v cannot make a 2-stem with.
#    In restrictions[2], S is the set of 2-frozensets of vertices that v cannot make a 3-stem with.
def NS_generator_with_enforced_planarity(roots_by_region,restrictions=False):
    if not restrictions:
        edge_restrictions = {x:set([]) for y in roots_by_region for x in y}
        ident_2_restrictions = {x:set([]) for y in roots_by_region for x in y}
        ident_3_restrictions = {x:set([]) for y in roots_by_region for x in y}
        restrictions = [edge_restrictions,ident_2_restrictions,ident_3_restrictions]
    
    if len(roots_by_region) == 0:
        yield [[],[],[]]
        return
    else:
        #We examine the last region.
        roots = roots_by_region[-1][:]
        
        if len(roots) < 2:
            for substructure in NS_generator_with_enforced_planarity(roots_by_region[:-1],restrictions=restrictions):
                structure = [[x[:] for x in y] for y in substructure]
                yield structure
            return
        
        else:
            #We examine the last root.

            #Case 1:  Condition on the last root being in a three-stem.
            for two in combinations(range(len(roots)-1),2):
                
                if frozenset([roots[x] for x in two]) in restrictions[2][roots[-1]]:
                    continue
                
                new_roots_by_region = [x[:] for x in roots_by_region[:-1]]
                new_roots_by_region.append(roots[:two[0]])
                new_roots_by_region.append(roots[two[0]+1:two[1]])
                new_roots_by_region.append(roots[two[1]+1:-1])
                for substructure in NS_generator_with_enforced_planarity(new_roots_by_region,restrictions=restrictions):
                    structure = [[x[:] for x in y] for y in substructure]
                    structure[2].append([roots[two[0]],roots[two[1]],roots[-1]])
                    yield structure

            #Case 2:  Condition on the last root being in a new edge.
            for one in [x for x in range(len(roots)-1) if not roots[x] in restrictions[0][roots[-1]]]:
                new_roots_by_region = [x[:] for x in roots_by_region[:-1]]
                new_roots_by_region.append(roots[:one])
                new_roots_by_region.append(roots[one+1:-1])
                for substructure in NS_generator_with_enforced_planarity(new_roots_by_region,restrictions=restrictions):
                    structure = [[x[:] for x in y] for y in substructure]
                    structure[0].append([roots[one],roots[-1]])
                    yield structure

            #Case 3:  Condition on the last root being in a two-stem.
            for one in [x for x in range(len(roots)-1) if not roots[x] in restrictions[1][roots[-1]]]:
                new_roots_by_region = [x[:] for x in roots_by_region[:-1]]
                new_roots_by_region.append(roots[:one])
                new_roots_by_region.append(roots[one+1:-1])
                for substructure in NS_generator_with_enforced_planarity(new_roots_by_region,restrictions=restrictions):
                    structure = [[x[:] for x in y] for y in substructure]
                    structure[1].append([roots[one],roots[-1]])
                    yield structure


            #Case 4:  Condition on the last root being in a one-stem.
            new_roots_by_region = [x[:] for x in roots_by_region[:-1]]
            new_roots_by_region.append(roots[:-1])
            for substructure in NS_generator_with_enforced_planarity(new_roots_by_region,restrictions=restrictions):
                structure = [[x[:] for x in y] for y in substructure]
                yield structure

            return












#--------------------------
#  Subsection C:  Checking All Realizations
#--------------------------
















#specific to configurations rather than subgraphs

#Now outer_lists are multiple regions (list of lists).
def check_all_neighborhood_structures_for_FAS(edges,outer_lists,include_restrictions=False):
    begin = time.clock()
    G = Graph(edges)
    if not include_restrictions:
        edge_restrictions = {x:set([]) for x in G.vertices()}
        ident_2_restrictions =  {x:set([]) for x in G.vertices()}
        ident_3_restrictions =  {x:set([]) for x in G.vertices()}
    else:
        edge_restrictions,ident_2_restrictions,ident_3_restrictions = include_restrictions
    count = 0
    good_count = 0
    bad_count = 0
    keys = ['error', 'greedy', 'CNS', 'brute']
    count_dict = {s:0 for s in keys}
    bad_li = []
    
    
    li = [x[:] for x in outer_lists]
    twos = [x for y in li for x in y]
    dist1 = []
    dist2 = []
    for i in range(len(twos)):
        u = twos[i]
        for j in range(i+1,len(twos)):
            v = twos[j]
            d = G.distance(u,v)
            if d == 1:
                dist1.append(set([u,v]))
            elif d == 2:
                dist2.append(set([u,v]))
    
    close_couples = set([])
    for pair in dist1:
        for v in pair:
            edge_restrictions[v] |= pair
            ident_2_restrictions[v] |= pair
        close_couples |= set([frozenset(pair),])
    for pair in dist2:
        for v in pair:
            ident_2_restrictions[v] |= pair
        close_couples |= set([frozenset(pair),])
    for roots in li:
        for triple in combinations(roots,3):
            triple_couples = set([frozenset([triple[0],triple[1]])])
            triple_couples |= set([frozenset([triple[0],triple[2]])])
            triple_couples |= set([frozenset([triple[1],triple[2]])])
            if len(triple_couples & close_couples) > 2:
                ident_3_restrictions[triple[0]] |= set([frozenset(triple[1:]),])
                ident_3_restrictions[triple[1]] |= set([frozenset([triple[0],triple[-1]]),])
                ident_3_restrictions[triple[2]] |= set([frozenset(triple[:2]),])
    
    
    
    for idents in NS_generator_with_enforced_planarity(li,restrictions=[edge_restrictions,ident_2_restrictions,ident_3_restrictions]):
        #restriction_flag = False
        #for edge in idents[0]:
            #if edge[1] in edge_restrictions[edge[0]]:
                #restriction_flag = True
                #break
        #if restriction_flag:
            #continue
        #for pair in idents[1]:
            #if pair[1] in ident_2_restrictions[pair[0]]:
                #restriction_flag = True
                #break
        #if restriction_flag:
            #continue
        #for triple in idents[2]:
            #if set(triple[1:]) in [set(x) for x in ident_3_restrictions[triple[0]]]:
                #restriction_flag = True
                #break
        #if restriction_flag:
            #continue
        count += 1
    print " "*10+"There are "+str(count)+" neighborhood structures to check.  (Took "+ch.timestring(time.clock()-begin)+" to count.)"
    print " "*15+"Index      Good   Bad   Time (Cum.)   Identification"
    
    i = 0
    for idents in NS_generator_with_enforced_planarity(li,restrictions=[edge_restrictions,ident_2_restrictions,ident_3_restrictions]):
        #restriction_flag = False
        #for edge in idents[0]:
            #if edge[1] in edge_restrictions[edge[0]]:
                #restriction_flag = True
                #break
        #if restriction_flag:
            #continue
        #for pair in idents[1]:
            #if pair[1] in ident_2_restrictions[pair[0]]:
                #restriction_flag = True
                #break
        #if restriction_flag:
            #continue
        #for triple in idents[2]:
            #if set(triple[1:]) in [set(x) for x in ident_3_restrictions[triple[0]]]:
                #restriction_flag = True
                #break
        #if restriction_flag:
            #continue
        i += 1
        G_sq,f = core_square_graph(G,idents)
        y = ch.fChoosableNoPrint(G_sq,f,inprocess=4,print_mod=100,maxIS=True,rev=False)
        if y[0]:
            good_count += 1
        else:
            bad_count += 1
            bad_li.append(idents)
        count_dict[y[1]] += 1
        ti = ch.timestring(time.clock()-begin)
        print " "*(20-len(str(i)))+str(i)+" "*(10-len(str(good_count)))+str(good_count)+" "*(6-len(str(bad_count)))+str(bad_count)+" "*(17-len(ti))+ti+"   "+str(idents)
    
    #print "\n"*7
    #print "Total:",count
    #print "Bad:",bad_count
    #print "Good:",good_count
    if count == good_count:
        print "      >>> Good for all stem identifications!"
    else:
        print "      >>> Uh-oh!  Problem!"
    print
    #for s in keys:
        #print s+":",count_dict[s]
    end = time.clock()
    #print "\nTime:",ch.timestring(end-begin)
    #if good_count < count:
        #print "Problems:"
        #for x in bad_li:
            #print x
    #print
    #if good_count == count:
        ##print "Complete:  "+rc_str+" is GOOD."
        #return True
    #else:
        ##print "Complete:  "+rc_str+" is BAD."
        #return False
    return count,bad_count










#g is same as what gets yielded by makeGraph, but is entered manually
#g has form:  order,edges,spine,roots,faces
#open_spine is whether configuration is cxa-type.
#partition_restrictions is a dictionary holding vertex:{set of forbidden identification for vertex}.
#stem_restrictions is a list of three dictionaries:
#stem_restrictions[0] holds vert:{set of forbidden vertices to add edge to vert}
#stem_restrictions[1] holds vert:{set of forbidden vertices to form 2-stem-identification with vert}
#stem_restrictions[2] holds vert:{set of forbidden vertices (as 2-sets) to form 3-stem-identification with vert}
def check_all_realizations_from_initial_plane_graph(g,open_spine=False,partition_restrictions={},stem_restrictions=False):
    begin = time.clock()
    FAS_count = 0
    realization_count = 0
    bad_count = 0
    totals_str = "     FAS#               #NS       #Bad\n"
    for realization in configuration_FAS_generator(g,open_spine=open_spine,include_restrictions=partition_restrictions):
        root_lists = get_outer_region_roots(realization[0],realization[1])
        partition,id_dict = realization[2:4]

        #if stem_restrictions:#if coming from the c7a4x5x5 cases
            ##We need to make sure that the stem restrictions we were given are preserved, since an FAS might have identified some core vertices and shifted the labels.  (Even though the original base vertices stay in order, one new vertex from the additional face might displace a higher-labeled new vertex from the additional face.  And we are passing on restrictions with these vertices.)
            ##print "id_dict:",id_dict
            ##print "stem_restrictions:"
            ##print "stem_1"
            ##for v in stem_restrictions[0].keys():
                ##print v,":",stem_restrictions[0][v]
            ##print "stem_2"
            ##for v in stem_restrictions[1].keys():
                ##print v,":",stem_restrictions[1][v]
            ##print "stem_3"
            ##for v in stem_restrictions[2].keys():
                ##print v,":",stem_restrictions[2][v]
            ##print
            
            ##stem_rest_1 = {}
            ##for x in range(len(partition)):
                ##print "x:",x
                ##for y in partition[x]:
                    ##print "    y:",y
                    ##for z in stem_restrictions[1][y]:
                        ##print "        z:",z,"   id_dict[z]:",id_dict[z]
            
            #stem_rest_1 = {x:{id_dict[z] for y in partition[x] for z in stem_restrictions[0][y]} for x in range(len(partition))}
            #stem_rest_2 = {x:{id_dict[z] for y in partition[x] for z in stem_restrictions[1][y]} for x in range(len(partition))}
            #stem_rest_3 = {x:{frozenset([id_dict[w] for w in z]) for y in partition[x] for z in stem_restrictions[2][y]} for x in range(len(partition))}
            #new_stem_restrictions = [stem_rest_1,stem_rest_2,stem_rest_3]
        #else:
            #new_stem_restrictions = False
            
        new_stem_restrictions = stem_restrictions
        
        NS_count,bad_NS_count = check_all_neighborhood_structures_for_FAS(edges=realization[0].edges(),outer_lists=root_lists,include_restrictions=new_stem_restrictions)
        FAS_count += 1
        realization_count += NS_count
        bad_count += bad_NS_count
        FAS_count_str = str(FAS_count)
        NS_count_str = str(NS_count)
        bad_NS_count_str = str(bad_NS_count)
        totals_str += " "*(max(9-len(FAS_count_str),0))+FAS_count_str+" "*(max(18-len(NS_count_str),0))+NS_count_str+" "*(max(11-len(bad_NS_count_str),0))+bad_NS_count_str+"\n"
    end = time.clock()
    print "Done!\nSummary:"
    print totals_str
    print "Totals:"
    print "     #FAS     #Realizations       #Bad"
    FAS_count_str = str(FAS_count)
    realization_count_str = str(realization_count)
    bad_count_str = str(bad_count)
    print " "*(max(9-len(FAS_count_str),0))+FAS_count_str+" "*(max(18-len(realization_count_str),0))+realization_count_str+" "*(max(11-len(bad_count_str),0))+bad_count_str
    
    
    
    if bad_count > 0:
        print "\nNope!  Some realizations are not core-choosable."
    else:
        print "\nGood!  All realizations are core-choosable!"
    print "Time: "+ch.timestring(end-begin)











def check_all_configuration_realizations(config_str):
    print "Configuration:",config_str
    print "Checking all realizations for core-choosability.\n"
    g = makeGraph(config_str)
    if config_str[1]=='x':
        open_spine=True
    else:
        open_spine=False
    check_all_realizations_from_initial_plane_graph(g,open_spine=open_spine,partition_restrictions={},stem_restrictions=False)
















def check_all_realizations_from_expanded_c7a4x5x5_case(case):
    print "Expanding scope around c7a4x5x5:"

    order,edges,spine,roots,faces = makeGraph('c7a4x5x5')#spine won't be important here.
    open_spine = False
    roots.sort()

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
        stem_1[v] = set([])
        stem_2[v] = set([])
        stem_3[v] = set([])
    
    

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

    del roots[roots.index(base_border[0])]
    del roots[roots.index(base_border[-1])]
    roots.extend(range(old_order,order))

    new_face = base_border[:]
    new_face.extend(range(old_order,order))
    faces.append(new_face)

#     Graph(edges).show()
    g = (order,edges,spine,roots,faces)
#     print g
    #Initialize dictionaries for the new vertices:
    for v in range(old_order,order):
        partition_di[v] = set(faces[-1])
        #We could further restrict things that new vertices can identify with, but we're not going to.
        stem_1[v] = set([])
        stem_2[v] = set([])
        stem_3[v] = set([])
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

    check_all_realizations_from_initial_plane_graph(g,open_spine=open_spine,partition_restrictions=partition_di,stem_restrictions=[stem_1,stem_2,stem_3])





















































































































































#------------------------------------------------------
#SECTION 4:  The Target Set
#------------------------------------------------------











HRCs =[
'c3a3',
'cxa3x3',
'cxa3xx3',
'cxa3xxx3',
'cxa3xx5x3',
'c3a4',
'cxa3x4',
'cxa3x54',
'cxa3xx54',
'cxa3xxx54',
'cxa375',
'cxa3x555',
'c3a6',
'c4a4',
'cxa454',
'cxa464',
'cxa474',
'cxa4x54',
'cxa4x55x4',
'cxa4x4x4x4',
'cxa4x45',
'cxa456',
'cxa555555',
'cxa546',
'c4a55',
'c4a56',
'c4a57',
'c4a66',
'c4a585',
'c4a676',
'c4a686',
'c5a555',
'c5a556',
'c5a565',
'c5a566',
'c5a575',
'c5a656',
'c5a666',
'c5a5x65',
'c5a5585',
'c5a5x66',
'c5a55x6',
'c5a56x6',
'c5a6x66',
'c7a3xx4',
'c7a3xx5',
'c7a4xx4',
'c7a4x55',
'c7a4xx55',
'c7a4x5xx5',
'c7a55x5',
'c8a3xx4',
'c8a3xxx4',
'c8a3x55x55',
'c8a45xx4',
'c8a455x5',
'c8a4x5x4x4',
'c9a3xx4',
'c9a3xxx4',
'c9a455x4',
'c9a455x5',
'c9a545x5',
'cxa4x555x4',
'cxa4x5555',
'cxa4555',
'c5a4x55',
'cxa355',
'cxa356',
'cxa535',
'cxa35x4',
'c3a55',#Same as cxa535.  Need to fix RC dictionary maker.
'c3a57',
'c8a35x5',
'c8a35xx5',
'c8a35xxx5',
'c8a35xxxx5',
'c9a35x5'
    ]







UsedInHappy3 = set(['c3a3','c3a4','c3a55','c3a57','c3a6'])
UsedInHappy4 = set(['c3a4','c4a4','c4a55','c4a56','c4a57','c4a585','cxa546','c4a66','c4a676','c4a686'])
UsedInHappy5 = set(['c3a3','c3a4','cxa3x3','cxa3x4','c3a55','cxa355','c3a6','cxa356','c4a4','cxa454','c4a55','c4a56','c5a4x55','cxa456','c5a555','c5a556','c5a575','c5a5585','c5a55x6','c5a565','c5a5x65','c5a666','c5a56x6','c5a656','c5a566','c5a5x66','c5a6x66'])
UsedInDischarging4 = set(['cxa464'])
UsedInDischarging5 = set(['cxa355','cxa356'])


















def initializeRedConList(huge,exclusions):
    RedConList = []
    for s in HRCs:
        if not s in exclusions:
            RedConList.append(s)

    #max_chain_length is used primarily for storing cxa-type configurations -- it gives a cap on where to stop.
    #Fixed central face lengths give a natural, but an x-face doesn't.
    #We can just choose a universal cap (like 12) instead if we want to.
    max_chain_length = max([len(x[x.index('a')+1:]) for x in RedConList])

    #def makeInventories(StandardRedConList):
    #------------------------------------------------------------------------------------------------
    central_options = ['x',]
    central_options.extend(range(3,huge+1))
    contain_options = range(3,huge+1)
    # length_options = range(1,max_chain_length+1)

    ChainsDict = {}
    for cent in central_options:
        ChainsDict[cent] = {}
        for contain in contain_options:
            ChainsDict[cent][contain] = {}
            if cent=='x':
                for length in range(1,max_chain_length+1):
                    ChainsDict[cent][contain][length] = []
            else:
                for length in range(1,cent+1):
                    ChainsDict[cent][contain][length] = []

    for s in RedConList:
        if s[1]=='x':
            #Special cases:  if chain is small enough, it can go through the central face.
            if len(s)==5:
                y,z = s[3],s[4]
                ChainsDict[int(y)][int(z)][1].append(z)
                if y != z:
                    ChainsDict[int(z)][int(y)][1].append(y)
            elif len(s)==6 and not s[4]=='x':
                for z in {s[3],s[5]}:
                    ChainsDict[int(s[4])][int(z)][3].append(s[3]+'x'+s[5])
            #The normal case (chain is entirely on outside) is performed on s even if a special case applied.
            for z in set(s[3:]):
                if z != 'x':
                    ChainsDict['x'][int(z)][len(s[3:])].append(s[3:])
        else:
            #Special cases:  if there are only one or two (one-apart) outer faces, they can be exterior chains.
            #If there are at most three (adjacent) outer faces, they can be rotated to have different central faces.
            if len(s)==4:
                y,z = s[1],s[3]
                ChainsDict[int(y)][int(z)][1].append(z)
                if y != z:
                    ChainsDict[int(z)][int(y)][1].append(y)
                for w in {int(y),int(z)}:
                    ChainsDict['x'][w][2].append(y+z)
            elif len(s)==5:
                w,y,z = s[1],s[3],s[4]
                if w==y and y==z:
                    ChainsDict[int(w)][int(w)][2].append(w+w)
                elif w==y:
                    ChainsDict[int(w)][int(w)][2].append(w+z)
                    ChainsDict[int(w)][int(z)][2].append(w+z)
                    ChainsDict[int(z)][int(w)][2].append(w+w)
                elif w==z:
                    ChainsDict[int(w)][int(w)][2].append(w+y)
                    ChainsDict[int(w)][int(y)][2].append(w+y)
                    ChainsDict[int(y)][int(w)][2].append(w+w)
                elif y==z:
                    ChainsDict[int(y)][int(w)][2].append(w+y)
                    ChainsDict[int(y)][int(y)][2].append(w+y)
                    ChainsDict[int(w)][int(y)][2].append(y+y)
                else:
                    ChainsDict[int(w)][int(y)][2].append(y+z)
                    ChainsDict[int(w)][int(z)][2].append(y+z)
                    ChainsDict[int(y)][int(w)][2].append(w+z)
                    ChainsDict[int(y)][int(z)][2].append(w+z)
                    ChainsDict[int(z)][int(w)][2].append(w+y)
                    ChainsDict[int(z)][int(y)][2].append(w+y)
            elif len(s)==6:
                if s[4]=='x':
                    w,y,z = s[1],s[3],s[5]
                    ChainsDict[int(w)][int(y)][3].append(y+'x'+z)
                    if y!=z:
                        ChainsDict[int(w)][int(z)][3].append(y+'x'+z)
                    for v in {int(w),int(y),int(z)}:
                        ChainsDict['x'][v][3].append(y+w+z)
                else:
                    v,w = s[1],s[4]
                    y,z = s[3],s[5]
                    for u in {int(w),int(y),int(z)}:
                        ChainsDict[int(v)][u][3].append(y+w+z)
                    if v!=w:
                        for u in {int(v),int(y),int(z)}:
                            ChainsDict[int(w)][u][3].append(y+v+z)
            #The normal case is when the outer faces span at least four spots.  No rotations or interior stuff.
            else:
                for z in set(s[s.index('a')+1:]):
                    if z!='x':
                        ChainsDict[int(s[1:s.index('a')])][int(z)][len(s[s.index('a')+1:])].append(s[s.index('a')+1:])
    #------------------------------------------------------------------------------------------------


    #Make a dictionary to put things back in equivialence class.
    RCID = {}
    for s in RedConList:
        if len(s)==4:
            #cJaK
            RCID[s] = s
            J = s[1]
            K = s[3]
            RCID['c'+K+'a'+J] = s
            RCID['cxa'+J+K] = s
            RCID['cxa'+K+J] = s
        elif len(s)==5:
            #cJaKL
            #(cxaJK is not the form used in StandardRedCon)
            RCID[s] = s
            J = s[1]
            K = s[3]
            L = s[4]
            RCID['c'+J+'a'+L+K] = s
            RCID['c'+K+'a'+J+L] = s
            RCID['c'+K+'a'+L+K] = s
            RCID['c'+L+'a'+J+K] = s
            RCID['c'+L+'a'+K+J] = s
        elif len(s)==6:
            #cxaJKL or cKaJxL, or cHaJKL
            H = s[1]
            J = s[3]
            K = s[4]
            L = s[5]
            if 'x' in [H,K]:
                RCID['cxa'+J+K+L] = s
                RCID['cxa'+L+K+J] = s
                RCID['c'+K+'a'+J+'x'+L] = s
                RCID['c'+K+'a'+L+'x'+J] = s
            else:
                RCID['c'+H+'a'+J+K+L] = s
                RCID['c'+H+'a'+L+K+J] = s
                RCID['c'+K+'a'+J+H+L] = s
                RCID['c'+K+'a'+L+H+J] = s
        else:
            RCID[s] = s
            pre = s[:3]
            post = s[3:]
            post = post[::-1]
            RCID[pre+post] = s
        
    return RedConList, ChainsDict, RCID, max_chain_length




























































#------------------------------------------------------
#SECTION 5:  Automated Discharging
#------------------------------------------------------






#Generates blocks of contiguous (5-)-faces which contain no reducible cxa-type configuration.
def nonreducible_block_generator(ChainsDict,RCID):
    used_basket = set([])
    stack = [3,]
    i = 1
    low = 3
    thresh = 5 #All face lengths in the stack will be at least low and at most thresh.
    
    while stack[0] <= thresh:
        #Does the current configuration have any reducible pieces?
        x = reducible_check_last_face_large_face(stack,ChainsDict,RCID)
        if not x[0]:
            #If no reducible pieces, then yield and keep going.
            #(One option is to check that the last face is not smaller than the first -- this would be a 
            #repeat block, just given in reverse.  However, we might actually want to use these repeats,
            #so we are leaving them here.)
            yield stack
            stack.append(low)
            i += 1
            continue
            
        else:
            #If reducible pieces, then make note in the list of used configurations and then backtrack.
            used_basket.add(x[1])
            
        #Backtrack
        while stack and stack[-1] >= thresh:
            del stack[-1]
            i -= 1
        if stack:
            stack[-1] += 1
        else:
            break
    
    yield None
    yield used_basket
    return








#Checking only an external stack, assuming no wrap-around.
def reducible_check_last_face_large_face(stack,ChainsDict,RCID):
    stack_copy = stack[:]
    
    for i in range(1,len(stack)+1):#i will be length of exterior chain.
        sub = stack_copy[-i:]
#         #First, make sure neither endpoint of the substack is an 'x'.  If one is, move on to the next substack.
#         if 'x' in [sub[0],sub[-1]]:
#             continue
# ^^Won't happen, because the setting I'm using this for is generating all outside strings for large faces -- no 'x's.
        #Check the reducibility of only the outside faces.
        try:
            for rc in ChainsDict['x'][stack[-1]][i]:
                if descriptor_match(sub,rc):
                    return (True,str_converter('x',rc,RCID))
        except KeyError:
            pass
    
    return (False,None)







def str_converter(cent,rc,RCID):
    s = "c"+str(cent)+"a"
    for x in rc:
        s += str(x)
    return RCID[s]




def descriptor_match(described_object,descriptor_candidate):
    if len(described_object)!=len(descriptor_candidate):
        print "Error!  Lengths of input do not match!  Input:",described_object,descriptor_candidate
        return
    bad = False
    for i in range(len(described_object)):
        if not descriptor_candidate[i] in ['x',str(described_object[i])]:
            bad = True
            break
    if bad:
        bad = False
        rev_object = described_object[:]
        rev_object.reverse()
        for i in range(len(rev_object)):
            if not descriptor_candidate[i] in ['x',str(rev_object[i])]:
                bad = True
                break
        if bad:
            return False
    return True









#Generates full stacks around a given face length, but only where small faces are specified and (6+)-faces are marked with 'x'.  
#(Ignores reducibility of configurations when 6,7,8-faces are present.)
def nonreducible_block_stack_generator(central_face_length,nonreducibles,max_chain_length,ChainsDict,RCID):
    used_basket = set([])
    pieces = nonreducibles[:]
    thresh = len(pieces)
    pieces.append(['x'])
    sizes = {p:len(pieces[p]) for p in range(thresh+1)}
    stack = [0,thresh]
    actual_stack = pieces[0][:]
    actual_stack.append('x')
#     i = 1
    
    while stack[0] < thresh:
        if len(actual_stack) <= central_face_length:#Might actually be too long.  If so, just backtrack.
            #Does the current configuration have any reducible pieces?
    #         print "> > > >",actual_stack
            #If the last piece is not just 'x' (if stack[-2]<thresh), then run sizes[stack[-2]] many reducible_check_last_face runs.
            #If the last piece is just 'x' (stack[-2]==thresh), then there is nothing to check that we haven't already.  Add and continue.
            back = False
            if stack[-2] < thresh:
                ell = sizes[stack[-2]]
                for i in range(ell):
                    x = reducible_check_last_face(central_face_length,actual_stack[:i-ell],max_chain_length,ChainsDict,RCID)
                    if x[0]:
                        used_basket.add(x[1])
                        back = True
                        break
            if not back:#If no new reducible parts showed up.
                #If no reducible pieces, is the configuration full?
                if len(actual_stack) < central_face_length:#If not full, then add a new face and go back to top of loop.
                    stack.append(stack[0])
                    actual_stack.extend(pieces[stack[0]])
                    stack.append(thresh)
                    actual_stack.append('x')
    #                 i += 1
                    continue
                else:
                    yield actual_stack#If full, then yield and then backtrack.
            
        #Backtrack
#         try:
        while stack[-1] >= thresh:
            del actual_stack[-1]#The last entry in actual stack is 'x'.
            del stack[-1]
#             i -= 1
#         except IndexError:#stack is empty
#             print 'meow'
#             yield None
#             yield used_basket
#             return
        for i in range(sizes[stack[-1]]):
            del actual_stack[-1]
        stack[-1] += 1
        actual_stack.extend(pieces[stack[-1]])
        if stack[-1] < thresh:
            stack.append(thresh)
            actual_stack.append('x')
#     print 'woof'
    yield None
    yield used_basket
    return
#     return










#Could bring over ChainsDict['x'] into all of ChainsDict[k]'s (when defining ChainsDict).
#This would make the case distinction unnecessary here.
#^^Above comment would be true if the cxa chains were also reducible when wrapped around and 
#the end faces were actually adjacent.  But if this is the case, we would just put in the 
#configuration with the appropriate central face.  So we should probably keep it like this.
def reducible_check_last_face(central_face_length,stack,max_chain_length,ChainsDict,RCID):
    spot = len(stack)-1
    stack_copy = stack[:]
    if len(stack) < central_face_length:
        for i in range(central_face_length-len(stack)):
            stack_copy.append('x')
    for i in range(spot):
        stack_copy.append(stack_copy[i])
        #The highest position we ever use below is central_face_length+spot-2, except in the full stack case.
        #Then the highest position is central_face_length+spot-1.
    
#     print stack
#     print stack_copy
    
    #The next portion of code checks for reducible subchains when the stack is not full.
    #(We need to handle this separately because chains of the form cxaJ...K have not technically
    #been checked for the case when the J- and K-faces are adjacent.)
    for i in range(1,min(central_face_length,max_chain_length+1)):#i will be length of exterior chain.
        for j in range(max(0,i-1-spot),i):
            sub = stack_copy[spot-i+1+j:spot+j+1]
            #First, make sure neither endpoint of the substack is an 'x'.  If one is, move on to the next substack.
            if 'x' in [sub[0],sub[-1]]:
                continue
            first_nonx = -1
            for z in sub:
                if z != 'x':
                    first_nonx = z
                    break
            if first_nonx == -1:
                continue
            #Check the reducibility of only the outside faces.
#             print "Red1"
            for rc in ChainsDict['x'][first_nonx][i]:
                if descriptor_match(sub,rc):
#                     print "meow1\n"
                    return (True,str_converter('x',rc,RCID))
            #Now check the reducibility of the outside faces conditioned on the central face.
#             print "Red2"
            for rc in ChainsDict[central_face_length][first_nonx][i]:
                if descriptor_match(sub,rc):
#                     print "meow2\n"
                    return (True,str_converter(central_face_length,rc,RCID))
    
    #Now we handle the full stack case.
    if len(stack)==central_face_length:
        #i=central_face_length
        for j in range(central_face_length):
            sub = stack_copy[j:j+central_face_length]
            #stack_copy might have 'x's.
            first_nonx = -1
            for z in sub:
                if z != 'x':
                    first_nonx = z
                    break
            if first_nonx == -1:
                continue
            #No cxa-type chains here.
#             print "Red3"
            for rc in ChainsDict[central_face_length][first_nonx][central_face_length]:
                if descriptor_match(sub,rc):
#                     print "meow4\n"
                    return (True,str_converter(central_face_length,rc,RCID))
    
#     print "meow3\n"
    return (False,None)












def stack2str(stack):
    s = ""
    for a in stack:
        s += str(a)
    return s










#for large faces; assumes length at least 7.
#for full stacks.
#A 3-face pulls 1 charge unless it is next to a 5-face -- then 3/2 charge.
#A 4-face pulls 1 charge if the central face length is at least 9 and its neighbors are (6-)-faces;
# otherwise,it pulls at most 2/3 charge.  Since (6+)-face-neighbors are not specified except by 'x's, we 
# check for small faces two faces away.
#A 5-face pulls at most 1/2 charge if the central face length is at least 9 and both its neighbors 
# are 5-faces.  If we use crowns as reducible configurations (i.e., c5a55(5)55), then we can further 
# restrict that this 5-face does not also have a 5-face two spots away.  Otherwise, pulls at most 1/3.
#An alternative way of tabulating charge is by calculating charge by hand for nonreducible blocks.
def final_charge(central_face_length,stack):
    charge = Fraction(0)
    
    ##Checks if the configuration contains (is) c7a4x5x5.  If it is, then we have special
    ##knowledge:  4 pulls 1/2 and 5s pull 1/4.  So the 7-face is happy.
    #if central_face_length == 7 and 4 in stack and 5 in stack:
        #longstack = stack[:]
        #longstack.extend(longstack[:-1])
        #fours = [x for x in range(central_face_length) if stack[x] == 4]
        #fives = [x for x in range(len(longstack)) if longstack[x] == 5]
        #for i in fours:
            #if i+2 in fives and i+4 in fives:
                #return charge
            #if i+3 in fives and i+5 in fives:
                #return charge
    
    charge += central_face_length-6
    
    threes = [x for x in range(central_face_length) if stack[x]==3]
    for i in threes:
        if 5 in [stack[i-1],stack[i+1-central_face_length]]:
            charge -= Fraction(3,2)
        else:
            charge -= 1
    
    fours = [x for x in range(central_face_length) if stack[x]==4]
    for i in fours:
        if i-2 in fours or i+2 in fours or i-2+central_face_length in fours or i+2-central_face_length in fours:
            charge -= Fraction(2,3)
        elif central_face_length<9:# or stack[i-1]>6 or stack[i+1-central_face_length]>6:#Now only using 'x' for (6+)-faces.
            charge -= Fraction(2,3)
        else:
            charge -= 1
    
    fives = [x for x in range(central_face_length) if stack[x]==5]
    for i in fives:
        nbrs = {y:len([x for x in [stack[i-1],stack[i+1-central_face_length]] if x==y]) for y in [3,4,5]}
        if nbrs[3]>0:
            charge -= Fraction(1,4)
        elif nbrs[5]>=2 and central_face_length>=9:
            charge -= Fraction(1,2)
        else:
            charge -= Fraction(1,3)
    
    return charge











#for huge faces; assumes length at least 9.
#for blocks of (5-)-faces.
def charge_pull(stack):
    charge = Fraction(0)
    
    threes = [x for x in range(len(stack)) if stack[x]==3]
    for i in threes:
        boo = False
        if len(stack)==1:
            pass
        elif i==0 and stack[1]==5:
            boo = True
        elif i==len(stack)-1 and stack[-2]==5:
            boo = True
        elif 5 in [stack[i-1],stack[i+1]]:
            boo=True
        if boo:
            charge += Fraction(3,2)
        else:
            charge += 1
    
    fours = [x for x in range(len(stack)) if stack[x]==4]
    for i in fours:
        charge += 1
    
    fives = [x for x in range(len(stack)) if stack[x]==5]
    for i in fives:
        if len(stack)==1:
            nbrs = {3:0,4:0,5:0}
        elif i==0:
            nbrs = {y:len([x for x in [stack[1]] if x==y]) for y in [3,4,5]}
        elif i==len(stack)-1:
            nbrs = {y:len([x for x in [stack[-2]] if x==y]) for y in [3,4,5]}
        else:
            nbrs = {y:len([x for x in [stack[i-1],stack[i+1]] if x==y]) for y in [3,4,5]}
        if nbrs[3]>0:
            charge += Fraction(1,4)
        elif nbrs[5]>=2:
            charge += Fraction(1,2)
        else:
            charge += Fraction(1,3)
    
    return charge











#size="s":  sm=6, la=7
#size="L":  sm=5, la=6
#Want to find the classes of configurations (given by a single minimal element) which 
# are equal on small faces (length at most sm) and majorize the minimal element on large
# faces (length at least la).
#
def simmer_once(input_list,size="s"):
    basket = []
    for s in input_list:
        new = True
        for t in basket:
            if less_than(t,s,size=size):#t<s
                #done with s; no need to add it
                new = False
                break
            if less_than(s,t,size=size):#s<t
                #remove t
                basket.remove(t)
                #add s
                basket.append(s)
                new = False
                break
        if new:
            #add s
            basket.append(s)
    return basket






def simmer(input_list,size="s"):
    current_list = input_list[:]
    while True:
        a = len(current_list)
        new_list = list(simmer_once(current_list))
        if len(new_list)==a:
            return new_list
        else:
            current_list = new_list








#assumes single-digit face lengths
#not strict.  should return True if equal
def less_than(li_1,li_2,size="s"):
    if size=="s":
        sm,la = 6,7
    if size=="L":
        sm,la = 7,8
    
    #Outside round:  check to see if objects are equal on the small faces (up to rotation).
    li_1_ext = ""
    for x in li_1:
        li_1_ext += str(x)
    for x in li_1[:-1]:
        li_1_ext += str(x)
    
    li_2_ext = ""
    for x in li_2:
        li_2_ext += str(x)
    for x in li_2[:-1]:
        li_2_ext += str(x)


    smalls = [x for x in range(len(li_1)) if li_1[x]<=sm]
    larges = [x for x in range(len(li_1)) if li_1[x]>=la]
    
    key = ""
    for i in range(len(li_1)):
        if i in smalls:
            key += str(li_1[i])
        else:
            key += "?"
    
    for i in fnmatch_gen(key,li_2_ext):
        #Inside round:  for each match-up of the small faces, see if the large faces are majorized.
        maj = True
        for j in larges:
            if li_1_ext[j] > li_2_ext[i+j]:#string comparison, not int.
                maj = False
                break
        if maj:
            return True
    #Do it again for reversals as well.
    li_2_ext = li_2_ext[::-1]
    for i in fnmatch_gen(key,li_2_ext):
        #Inside round:  for each match-up of the small faces, see if the large faces are majorized.
        maj = True
        for j in larges:
            if li_1_ext[j] > li_2_ext[i+j]:#string comparison, not int.
                maj = False
                break
        if maj:
            return True
    
    return False














def fnmatch_gen(substring, string):
    index = 0
    while index <= len(string)-len(substring):
        if fnmatch.fnmatch(string[index:index+len(substring)],substring):
            yield index
        index += 1
    return












def run_discharging_analysis(huge,exclude=[]):
    #Initialize list and dictionaries of reducible configurations.
    RedConList, ChainsDict, RCID, max_chain_length = initializeRedConList(huge,exclude)
    UsedForSmallFacesHappy = UsedInHappy3 | UsedInHappy4 | UsedInHappy5
    UsedForDischargingCalculations = UsedInDischarging4 | UsedInDischarging5
    UsesDict = {s:[] for s in set(RedConList) | UsedForSmallFacesHappy | UsedForDischargingCalculations}
    for s in UsedInHappy3:
        UsesDict[s].append("3-Faces Are Happy")
    for s in UsedInHappy4:
        UsesDict[s].append("4-Faces Are Happy")
    for s in UsedInHappy5:
        UsesDict[s].append("5-Faces Are Happy")
    for s in UsedInDischarging4:
        UsesDict[s].append("4-Faces Pull Limited Charge When Close")
    for s in UsedInDischarging5:
        UsesDict[s].append("5-Faces Pull Certain Charge")
    
    #Initialize nonreducible_blocks.
    begin = time.clock()
    nonreducible_blocks = []
    it = nonreducible_block_generator(ChainsDict,RCID)
    x = it.next()
    while x != None:
        nonreducible_blocks.append(x[:])
        x = it.next()
    x = it.next()
    #print x
    for s in x:
        UsesDict[s].append("Nonreducible Blocks for All Faces")

    print "Nonreducible blocks of small faces (computed without a central face, \nwith reversals distinct):"
    for block in nonreducible_blocks:
        print stack2str(block)
    print "(There are %d.  Time:  %s.)"%(len(nonreducible_blocks),ch.timestring(time.clock()-begin))
    print

    unhappy_count = 0
    for k in range(7,huge):
        print "Now analyze central face length k = %d:"%(k)

        #Which blocks are still nonreducible for this length of central face?
        print "  Nonreducible blocks around a %d-face:"%(k)
        lap = time.clock()
        new_blocks = []
        used = set([])
        for block in nonreducible_blocks:
            boo = True
            for i in range(len(block)):
                x = reducible_check_last_face(k,block[:i+1],max_chain_length,ChainsDict,RCID)
                if x[0]:
                    boo = False
                    used.add(x[1])
                    break
            if boo:
                new_blocks.append(block[:])
                print "    "+stack2str(block)
            #if block == [4,5]:
                #print "boo:",boo
                #print "x:",reducible_check_last_face(k,block[:i+1],max_chain_length,ChainsDict,RCID)
        for s in used:
            UsesDict[s].append("Nonreducible Blocks for %d-Faces"%(k))
        print "    (There are %d.  Time: %s)"%(len(new_blocks),ch.timestring(time.clock()-lap))

        #How can these blocks form full stacks around this central face length?
        bag = []
        lap1 = time.clock()
        it = nonreducible_block_stack_generator(k,new_blocks,max_chain_length,ChainsDict,RCID)
        for x in it:
            if x==None:
                break
            bag.append(x[:])
        x = it.next()
        for s in x:
            UsesDict[s].append("Nonreducible Block Stacks for %d-Faces"%(k))
        print "  Number of full stacks (unsimmered) formed by the nonreducible \n  blocks: %d.  (Time to Generate: %s)"%(len(bag),ch.timestring(time.clock()-lap1))

        #Which ones pull too much charge?
        print "  Full stacks which pull too much charge (simmered):"
        lap3 = time.clock()
        unhappy_bag = []
        for stack in bag:
            if final_charge(k,stack) < 0:
                unhappy_bag.append(stack)
        sim_bag = simmer(unhappy_bag)
        for stack in sim_bag:
            print "    "+stack2str(stack)+" (Final Charge =",final_charge(k,stack),")"
        print "   (There are %d unhappy stacks.  Time (including simmer): %s)"%(len(sim_bag),ch.timestring(time.clock()-lap3))
        print
        unhappy_count += len(sim_bag)


    a = int(huge-6)
    b = int(huge)
    avgchg = Fraction(a,b)#This was very persnickety -- wouldn't accept Fraction(huge-6,huge).
    print "Finally, analyze faces of length at least %d."%(huge)
    print "  Such faces have at least %d charge to give away, and thus can afford \n  to give away an average charge of %s across each edge."%(int(huge-6),str(avgchg))
    print "  We calculate an upper bound on the total charge pull of each \n  nonreducible block, as well as the amount it would have to exceed the \n  average available charge for the central face:"
    print "    Block ... Pull ... Allowed ... Good?"
    bad_count = 0
    for block in nonreducible_blocks:
        s = stack2str(block)
        take = charge_pull(block)
        t = str(take)
        le = len(block)+1
        ma = le*avgchg
        u = str(ma)
        if take > ma:
            w = "Bad!  <---------"
            bad_count += 1
        else:
            w = "Good"
        print "    "+s+" "+"."*(8-len(s))+" "+t+" "+"."*(7-len(t))+" "+u+" "+"."*(10-len(u))+" "+w
    print "    (So we have %d bad blocks.)"%(bad_count)
    print
    
    

    print "-"*20+" SUMMARY "+"-"*20
    good = True
    s = ""
    t = ""
    if unhappy_count > 0:
        s = "NOT "
        good = False
    if bad_count > 0:
        t = "NOT "
        good = False
    print "Since there are %d unhappy stacks for faces of lengths 7 through %d, \nthese faces are "%(unhappy_count,huge-1)+s+"happy."
    print "Since there are %d bad blocks for faces of length at least %d, these \nfaces are "%(bad_count,huge)+t+"happy."
    
    small_happy_flag = False
    missing1 = []
    for s in UsedForSmallFacesHappy:
        if s not in RedConList:
            small_happy_flag = True
            missing1.append(s)
    discharge_calc_flag = False
    missing2 = []
    for s in UsedForDischargingCalculations:
        if s not in RedConList:
            discharge_calc_flag = True
            missing2.append(s)
    if small_happy_flag or discharge_calc_flag:
        print "*However, there is at least one problem:"
        if small_happy_flag:
            print "   The following configurations were used to show small faces are happy\n   (and maybe calculate charge pull) but are missing from the list:"
            for s in missing1:
                print "       "+s
        if discharge_calc_flag:
            print "   The following configurations were used to calculate charge pull,\n   but are missing from the list:"
            for s in missing2:
                print "       "+s
        good = False
    
    if good:
        print "So the discharging argument with this set of reducible configurations is complete!"
    else:
        print "So the discharging argument with this set of reducible configurations is NOT complete."
    print "Total time:  "+ch.timestring(time.clock()-begin)
    
    
    
    rebegin = time.clock()
    NotUsed = []
    print "\n"*7
    print "*"*20+" INFORMATION ON CONFIGURATIONS "+"*"*20
    for s in RedConList:
        #print s+":  ",UsesDict[s]
        if len(UsesDict[s])==0:
            NotUsed.append(s)
    print "Of the %d configurations in the list fed to this analysis (given \nbelow), %d were used."%(len(RedConList),len(RedConList)-len(NotUsed))
    print "       LIST:"
    for s in RedConList:
        print "           "+s
    print "\n"*3
    print "Information about each configuration from the list is given below:"
    for s in RedConList:
        print "\n"
        print "*"*10+" "+s+" "+"*"*10
        if len(UsesDict[s])>0:
            print s+" was used in the discharging analysis to show the following:"
            for t in UsesDict[s]:
                print "    * "+t
        else:
            print s+" was not used in the discharging analysis."

    return NotUsed











def run_discharging_analysis_no_print(huge,exclude=[]):
    #Initialize list and dictionaries of reducible configurations.
    RedConList, ChainsDict, RCID, max_chain_length = initializeRedConList(huge,exclude)
    UsedForSmallFacesHappy = UsedInHappy3 | UsedInHappy4 | UsedInHappy5
    UsedForDischargingCalculations = UsedInDischarging4 | UsedInDischarging5
    
    #Initialize nonreducible_blocks.
    nonreducible_blocks = []
    it = nonreducible_block_generator(ChainsDict,RCID)
    x = it.next()
    while x != None:
        nonreducible_blocks.append(x[:])
        x = it.next()

    unhappy_count = 0
    for k in range(7,huge):
        #Which blocks are still nonreducible for this length of central face?
        new_blocks = []
        for block in nonreducible_blocks:
            boo = True
            for i in range(len(block)):
                x = reducible_check_last_face(k,block[:i+1],max_chain_length,ChainsDict,RCID)
                if x[0]:
                    boo = False
                    break
            if boo:
                new_blocks.append(block[:])

        #How can these blocks form full stacks around this central face length?
        bag = []
        it = nonreducible_block_stack_generator(k,new_blocks,max_chain_length,ChainsDict,RCID)
        for x in it:
            if x==None:
                break
            bag.append(x[:])

        #Which ones pull too much charge?
        unhappy_bag = []
        for stack in bag:
            if final_charge(k,stack) < 0:
                return False
                unhappy_bag.append(stack)
        sim_bag = simmer(unhappy_bag)
        unhappy_count += len(sim_bag)


    a = int(huge-6)
    b = int(huge)
    avgchg = Fraction(a,b)#This was very persnickety -- wouldn't accept Fraction(huge-6,huge).
    bad_count = 0
    for block in nonreducible_blocks:
        take = charge_pull(block)
        le = len(block)+1
        ma = le*avgchg
        u = str(ma)
        if take > ma:
            return False
            bad_count += 1
    
    

    good = True
    if unhappy_count > 0:
        good = False
    if bad_count > 0:
        good = False
    
    small_happy_flag = False
    for s in UsedForSmallFacesHappy:
        if s not in RedConList:
            small_happy_flag = True
    discharge_calc_flag = False
    for s in UsedForDischargingCalculations:
        if s not in RedConList:
            discharge_calc_flag = True
    if small_happy_flag or discharge_calc_flag:
        good = False
    
    return good











