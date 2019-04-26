#------------------------------------------------------
# CONTENTS
#------------------------------------------------------


#SECTION 0:  Libraries and Comments


#SECTION 1:  Graph Essentials
#    graph_power
#    makeGraph
#    core_square_graph


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





#------------------------------------------------------
# Libraries, Comments, etc.
#------------------------------------------------------

import copy
# We use deepcopy in graph_power.

from sage.graphs.graph import Graph
# We use the Graph class to work with our realizations.

from itertools import combinations
#Used in:
#    NS_generator_no_enforced_planarity
#    NS_generator_with_enforced_planarity
#another?

import choosability as ch
# We use fChoosable as a black box.

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

        #In this case, vertices which are joining two specified outer faces cannot be identified with anything else on the spine.
        spine_list = open_spine[3:]
        forbidden = [x for x in range(1,len(spine_list)) if spine_list[x-1]!='x' and spine_list[x]!='x']
        for v in the_spine:
            if v in forbidden:
                di[v] |= set(the_spine)
            else:
                di[v] |= set(forbidden)

        #In addtion, vertices can't be identified which are too close.
        di[the_spine[0]] |= set(the_spine[0:3])
        di[the_spine[1]] |= set(the_spine[0:4])
        di[the_spine[-2]] |= set(the_spine[-4:])
        di[the_spine[-1]] |= set(the_spine[-3:])
        for i in range(2,len(the_spine)-2):
            di[the_spine[i]] |= set(the_spine[i-2:i+3])
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
            print "."*max(11-len(str(count)),0)+str(count)+" "*max(9-len(str(countp)),0)+str(countp)+"     "+timestring(time.clock()-begin)
    end = time.clock()
    t = end-begin
    print "Total #partitions:   ",count
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






















#remove close pair.triple restrictions entirely

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
    
    
    ##################################################3
    li = [x[:] for x in outer_lists]
    #twos = [x for y in li for x in y]
    #dist1 = []
    #dist2 = []
    #for i in range(len(twos)):
        #u = twos[i]
        #for j in range(i+1,len(twos)):
            #v = twos[j]
            #d = G.distance(u,v)
            #if d == 1:
                #dist1.append(set([u,v]))
            #elif d == 2:
                #dist2.append(set([u,v]))
    
    #close_couples = set([])
    #for pair in dist1:
        #for v in pair:
            #edge_restrictions[v] |= pair
            #ident_2_restrictions[v] |= pair
        #close_couples |= set([frozenset(pair),])
    #for pair in dist2:
        #for v in pair:
            #ident_2_restrictions[v] |= pair
        #close_couples |= set([frozenset(pair),])
    #for roots in li:
        #for triple in combinations(roots,3):
            #triple_couples = set([frozenset([triple[0],triple[1]])])
            #triple_couples |= set([frozenset([triple[0],triple[2]])])
            #triple_couples |= set([frozenset([triple[1],triple[2]])])
            #if len(triple_couples & close_couples) > 2:
                #ident_3_restrictions[triple[0]] |= set([frozenset(triple[1:]),])
                #ident_3_restrictions[triple[1]] |= set([frozenset([triple[0],triple[-1]]),])
                #ident_3_restrictions[triple[2]] |= set([frozenset(triple[:2]),])
    ##################################################3
    
    
    
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
    print " "*10+"There are "+str(count)+" neighborhood structures to check.  (Took "+timestring(time.clock()-begin)+" to count.)"
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
        ti = timestring(time.clock()-begin)
        print " "*(20-len(str(i)))+str(i)+" "*(10-len(str(good_count)))+str(good_count)+" "*(6-len(str(bad_count)))+str(bad_count)+" "*(14-len(ti))+ti+"   "+str(idents)
    
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
    #print "\nTime:",timestring(end-begin)
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













#remove close pairs/triples entirely

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
    print "Time: "+timestring(end-begin)












def check_all_configuration_realizations(config_str):
    print "Configuration:",config_str
    print "Checking all realizations for core-choosability.\n"
    g = makeGraph(config_str)
    if config_str[1]=='x':
        open_spine=config_str
    else:
        open_spine=False
    check_all_realizations_from_initial_plane_graph(g,open_spine=open_spine,partition_restrictions={},stem_restrictions=False)


















#currently blocking out stem restrictions for testing

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


#testing begins
#overwrite stem123
    #stem_1 = {x:set([]) for x in roots}
    #stem_2 = {x:set([]) for x in roots}
    #stem_3 = {x:set([]) for x in roots}
#testing ends


    check_all_realizations_from_initial_plane_graph(g,open_spine=open_spine,partition_restrictions=partition_di,stem_restrictions=[stem_1,stem_2,stem_3])







