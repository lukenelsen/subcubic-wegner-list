import copy
from sage.graphs.graph import Graph

mapping=reduce(lambda x,y: x+y,
    [chr(ord('0')+i) for i in range(10)]+
    [chr(ord('A')+i) for i in range(26)]+
    [chr(ord('a')+i) for i in range(26)]+
    ["@#"],
    "")
inverse_mapping=[-1]*256
for i,x in enumerate(mapping):
    inverse_mapping[ord(x)]=i


def encode_6bits(x):
    return mapping[x]


def decode_6bits(c):
    return inverse_mapping[ord(c)]


def encode_fgraph_string(G,f):
    n=G.num_verts()
    if n>=63:
        print "This graph has >=63 vertices and cannot be encoded in fgraph_string."
        raise NotImplementedError

    s=encode_6bits(n)+"_"

    # output f, each entry encoded in 6-bits
    for i in range(n):
        s+=encode_6bits(f[i])

    s+="_"

    # output the adjacency matrix
    val=0
    bits_avail_in_val=6
    for j in range(n):  # adj matrix is bit packed in colex order, just like in graph6 format
        for i in range(j):
            val<<=1
            val|=(1 if G.has_edge(i,j) else 0)
            bits_avail_in_val-=1

            if bits_avail_in_val==0:  # val has been filled with 6 bits
                s+=encode_6bits(val)
                val=0
                bits_avail_in_val=6

    if (bits_avail_in_val<6):
        # some bits are in val that still need to be outputted
        val<<=bits_avail_in_val  # shift bits so they are in the high position (bigendian)
        s+=encode_6bits(val)

    return s


def decode_fgraph_string(s):
    data=s.split("_")
    n=decode_6bits(data[0])
    f=[decode_6bits(data[1][i]) for i in range(n)]

    G=Graph()
    G.add_vertices(range(n))

    # read in the adjacency matrix
    cur=0
    val=decode_6bits(data[2][cur])
    mask=1<<5  # start with the high bit
    for j in range(n):  # adj matrix is bit packed in colex order, as in graph6 format
        for i in range(j):
            if val&mask!=0:  # test whether that bit is nonzero
                G.add_edge(i,j)
            mask>>=1
            if mask==0:  # mask has become 0
                cur+=1
                val=decode_6bits(data[2][cur])
                mask=1<<5

    return G,f













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




#edges is the list of edges in the original graph, before it is squared.
#name is the configuration code as a string.
class RedCon:

    def __init__(self,edges,name):
        G = Graph(edges)
        self.underlying_graph = G
        self.graph = graph_power(G,2)

        twos = [x for x in G.vertices() if G.degree(x)==2]
        f = []
        for v in range(G.order()):
            fval = 7
            if v in twos:
                fval -= 3
            for u in G.neighbors(v):
                if u in twos:
                    fval -=1
            f.append(fval)
        self.f = f[:]
        self.name = name
        self.info = "Name: " + name + "\nfGraph string: "+encode_fgraph_string(G,f) + "\nOrder: "+str(G.order())+",  Size: "+str(self.graph.size()) + "\nf: "+str(self.f) + "\nDisparity: "+str(sum(f)-G.order()-self.graph.size())

        self.underlying_graph_stems = Graph(self.underlying_graph.adjacency_matrix())
        for u in twos:
            v = self.underlying_graph_stems.order()
            self.underlying_graph_stems.add_edge([u,v])
            self.underlying_graph_stems.add_edge([v,v+1])
            self.underlying_graph_stems.add_edge([v,v+2])
        #No point in having an f for the stems, since all interior vertices have f=7, stems have f=3, and leaves have f=1.





#Take a configuration string and parses it into type_string, central_face_length, and remaining_list.  Then uses these to generate an edge list of the underlying graph.
#If prefix is 'cxa', does not take any central_face_length but proceeds similarly to type "c".
#type_string can be "d" (distance), "f" (facial chain), or "c" (central face).
#central_face_length is an int, the length of the first/central face in the configuration.
#remaining list is a list of ints and possibly the character 'e' or 'x' (if the configuration type is 'd' or 'c', respectively).
#Type 'd':  remaining_list should begin with 'e', continue with some number of 'e's, and then end with an int.  ('e' represents an edge increasing the distance.)
#Type 'f':  remaining_list should be a list of ints.
#Type 'c':  remaining_list should begin and end with ints, with other entries being ints or 'x'.  ('x' represents a skipped face.)
#Since all ints are face lengths, they should all be at least 3.
#Note:  type 'c' currently assumes all non-central faces have length less than 10.
def makeRedCon(configuration):
    
    if configuration[:3] == 'cxa':
        remaining_list = []
        for a in configuration[3:]:
            if a=='x':
                remaining_list.append(a)
            else:
                remaining_list.append(int(a))
        i=0
        E=[]
        last_base = 0
        last_stem = False  #last_stem marks whether or not we are building the next face adjacent to another non-central face.
        last_3 = False #last_3 marks whether or not the previous face was a 3-face with a previous face adjacent to it.
        while i < len(remaining_list):
            if remaining_list[i]=='x':
                E.append((last_base,last_base+1))
                last_base += 1
                last_stem = False
                last_3 = False
            elif not last_stem:#If we skipped a face (or are building our first face)
                #build length-1 new vertices
                j=0
                while j < remaining_list[i]-1:
                    E.append((last_base+j,last_base+j+1))
                    j += 1
                E.append((last_base,last_base+j))
                last_base += j
                last_stem = last_base-1
                #last_3 remains False
            else:#If we are building a face adjacent to a previous face
                if last_3:#assumes the next face being built is not a (6-)-face
                    #build first edge out from previous face
                    E.append((last_stem,last_stem+4))
                    last_stem += 4
                    #build length-4 new vertices
                    j=0
                    while j < remaining_list[i]-4:
                        E.append((last_stem+j,last_stem+j+1))
                        j += 1
                    #build last edge (the one on the ignored central face)
                    E.append((last_base,last_stem+j))
                    #adjust guides
                    last_base = last_stem+j
                    last_stem = last_base-1
                    last_3 = False
                else:
                    if remaining_list[i]==3:#assumes the previous face built was not a (6-)-face
                        #build first edge out from previous face
                        E.append((last_stem,last_stem+2))
                        #build last edge (the one on the ignored central face)
                        E.append((last_base,last_stem+2))
                        #adjust guides
                        last_base = last_stem+2
                        last_stem -= 1
                        last_3 = True
                    else:
                        #build first edge out from previous face
                        E.append((last_stem,last_stem+2))
                        last_stem += 2
                        #build length-3 new vertices
                        j=0
                        while j < remaining_list[i]-3:
                            E.append((last_stem+j,last_stem+j+1))
                            j += 1
                        #build last edge (the one on the ignored central face)
                        E.append((last_base,last_stem+j))
                        #adjust guides
                        last_base = last_stem+j
                        last_stem = last_base-1
                        last_3 = False
            i += 1
        return RedCon(E,configuration)
    
    
    #else
    type_string = configuration[0]
    d = {'d':'e','f':'v','c':'a'}
    s = configuration.index(d[type_string])
    central_face_length = int(configuration[1:s])
    remaining_list = []
    
    if type_string=='d':
        while 'e' in configuration[s:]:
            remaining_list.append('e')
            s += 1
        remaining_list.append(int(configuration[s:]))
    elif type_string=='f':
        try:
            while True:
                t = configuration[s+1:].index('v')
                remaining_list.append(int(configuration[s+1:s+1+t]))
                s += t + 1
        except ValueError:
            remaining_list.append(int(configuration[s+1:]))
    elif type_string=='c':
        for a in configuration[s+1:]:
            if a=='x':
                remaining_list.append(a)
            else:
                remaining_list.append(int(a))
        
    
    #Make the first/central face.
    E = [(0,1),]
    z = 2
    while z < central_face_length:
        E.append((z-1,z))
        z += 1
    E.append((0,z-1))
    
    
    #Now the rest.
    if type_string=='d':
        i=0
        while remaining_list[i]=='e':
            v = 1+max([max(x) for x in E])
            E.append((v-1,v))
            i += 1
        j = 0
        while j+1 < remaining_list[i]:
            E.append((v+j,v+j+1))
            j += 1
        E.append((v,v+j))
        
    elif type_string=='f':
        for i in range(len(remaining_list)):
            a = 1+max([max(x) for x in E])
            E.append((a-2,a))
            j=0
            while j<remaining_list[i]-3:
                E.append((a+j,a+j+1))
                j+= 1
            E.append((a-1,a+j))
        
    elif type_string=='c':
        i=0
        last_stem = False  #last_stem marks whether or not we are building the next face adjacent to another non-central face.
        last_3 = False #last_3 marks whether or not the previous face was a 3-face with a previous face adjacent to it.
        while i < len(remaining_list):
            if remaining_list[i]=='x':
                last_stem = False
            elif i < central_face_length-1:
                if last_3:
                    E.append((last_stem-1,last_stem+1))
                    j=0
                    while j < remaining_list[i]-5:#Assumes 3-faces are not adjacent to 3- or 4-faces.
                        E.append((last_stem+1+j,last_stem+1+j+1))
                        j+=1
                    E.append((i+1,last_stem+j+1))
                    last_stem += j+1
                    last_3 = False
                else:
                    if not last_stem:#If we skipped a face, then set up the first edge (which acts like the last edge of a previous face).
                        last_stem = 1+max([max(x) for x in E])
                        E.append((i,last_stem))
                    j=0
                    while j < remaining_list[i]-3:
                        E.append((last_stem+j,last_stem+j+1))
                        j+=1
                    E.append((i+1,last_stem+j))
                    last_stem += j
                if remaining_list[i] == 3 and i > 0 and remaining_list[i-1] != 'x':
                    last_3 = True
            else:#If we have all possible adjacencies and have reached the last one, then we finish off slightly differently.  Assume no 'x's in configuration.
                if remaining_list[0] == 3:#Assumes that if there is a 3-face, it will be listed first.  Also assumes that 3-faces are not near each other.
                    j=0
                    while j < remaining_list[i]-5:
                        E.append((last_stem+j,last_stem+j+1))
                        j+=1
                    E.append((last_stem+j,central_face_length+1))
                else:
                    j=0
                    while j < remaining_list[i]-4:
                        E.append((last_stem+j,last_stem+j+1))
                        j+=1
                    E.append((last_stem+j,central_face_length))
            i += 1

    return RedCon(E,configuration)



global ReducibleListD
global ReducibleListF
global ReducibleListC
global ReducibleList
ReducibleListD = ['d3e3','d3ee3','d3eee3','d3e4']
ReducibleListF = ['f3v3','f3v4','f3v5','f3v6','f4v4',
                 'f3v7v5','f3v7v6',
                 'f4v5v4','f4v5v5','f4v5v6','f4v6v4','f4v6v5','f4v6v6','f4v7v4',
                 'f5v4v6',
                 'f4v5v7v4',
                 'f5v5v5v5v5v5']
ReducibleListC = ['c3a777',
                 'c4a55','c4a56','c4a57','c4a66',
                 'c4a585','c4a676','c4a677','c4a686',
                 'c5a555','c5a556','c5a557','c5a565','c5a566','c5a567','c5a575','c5a576',
                 'c5a656','c5a657','c5a666','c5a667','5a676',
                 'c5a55x6',
                 'c6a4xx4',
                 'c7a3xx4','c7a3xx5','c7a4xx4','c7a55x5',
                 'c8a3xx4','c8a3xxx4','c8a3x55x55','c8a4xx4','c8a4xxx4',
                 'c9a3xx4','c9a3xxx4',
                 'c10a3xxx4','c10a3xxxx3','c10a3xxxx4',
                 'c11a3xxxx3',
                 'c12a3xxxxx3']
ReducibleList = ReducibleListD
ReducibleList.extend(ReducibleListF)
ReducibleList.extend(ReducibleListC)

global RedDict
RedDict = {}

#This makes all of the RedCon objecs for the configurations listed in ReducibleList, and associates them with their string in RedDict.
def createLibRCs():
    global RedDict
    for x in ReducibleList:
        RedDict[x] = makeRedCon(x)







#Returns problems (nonreducible and not happy) majorized by ceiling and with at least one (5-)-face.
def findProblems(central_length,ceiling):
    nonreducibles_unsimmered = findNonReducible(central_length,ceiling)
    nonreducibles = simmer(nonreducibles_unsimmered)
    problems = []
    for t in nonreducibles:
        if not happy(central_length,t):
            problems.append(t)
    return problems

#Returns non-reducible configurations whose outer faces are majorized by ceiling and with at least one (5-)-face.
def findNonReducible(central_length,ceiling):
    nonreducible = []
    stack=[3,]
    i=0
    if central_length < 6:
        a = 100#artificial infinity
    else:
        a = 6
    while stack[0] < a:
        full =  (i == central_length-1)
        if reducible(central_length,stack,full):
            back = True
        else:
            back = False
        
        if full and not back:
            nonreducible.append(tuple(stack))
            back = True
        
        if back:
            while stack[i] == ceiling[i]:
                del stack[i]
                i -= 1
                if i < 0:
                    return nonreducible
            stack[i] += 1
        else:
            stack.append(stack[0])
            i += 1
    return nonreducible



#Generate and check possible subconfigurations.
#If any configurations are added to type 'd' or 'f', or to 'c' with 'x's, then the code will need to be adjusted accordingly.
def reducible(central_length,stack,full):
    global ReducibleListF
    global ReducibleListC
    global ReducibleList
    
    #First, check for subconfigurations with two faces (type 'd' and type 'fAvB' where A is 3 or 4)
    #Central face included
    if central_length == 3 and [x for x in stack if x < 7]:
        return True
    if central_length == 4 and 4 in stack:
        return True
    if central_length < 7 and 3 in stack:
        return True
    #Outer faces only
    spots3 = [i for i,x in enumerate(stack) if x==3]
    if len(spots3) > 1:
        if min([spots3[j+1]-spots3[j] for j in range(len(spots3)-1)]) < 5:
            return True
        if full and spots3[0] + central_length - spots3[-1] < 5:
            return True
        spots4 = [i for i,x in enumerate(stack) if x==4]
        for i in spots4:
            for j in spots3:
                if -3 < i-j < 3:
                    return True
        if full and ( spots3[0] + central_length - spots4[-1] < 3 or spots4[0] + central_length - spots3[-1] < 3 ):
            return True
    
    #Second, check for other type 'f' subconfigurations.
    if len(stack) > 2:
        #Three faces -- central face included
        if central_length == 4:
            if stack[0] == 5 and stack[2] == 6:
                return True
            if stack[0] == 6 and stack[2] == 5:
                return True
            if full:
                if stack[1] == 5 and stack[3] == 6:
                    return True
                if stack[1] == 6 and stack[3] == 5:
                    return True
        if central_length in [5,6]:
            for i in range(len(stack)-2):
                if stack[i] == 4 and stack[i+2] in [4,5,6]:
                    return True
                if stack[i+2] == 4 and stack[i] in [4,5,6]:
                    return True
            if len(stack) > central_length-2:
                if stack[0] == 4 and stack[central_length-2] in [4,5,6]:
                    return True
                if stack[central_length-2] == 4 and stack[0] in [4,5,6]:
                    return True
            if full:
                if stack[1] == 4 and stack[central_length-1] in [4,5,6]:
                    return True
                if stack[central_length-1] == 4 and stack[1] in [4,5,6]:
                    return True
        if central_length == 7:
            for i in range(len(stack)-2):
                if stack[i] == 3 and stack[i+2] in [5,6]:
                    return True
                if stack[i+2] == 3 and stack[i] in [5,6]:
                    return True
                if stack[i] == 4 and stack[i+2] == 4:
                    return True
            if len(stack) > 5:
                if stack[0] == 3 and stack[5] in [5,6]:
                    return True
                if stack[5] == 3 and stack[0] in [5,6]:
                    return True
                if stack[0] == 4 and stack[5] == 4:
                    return True
            if full:
                if stack[0] == 3 and stack[6] in [5,6]:
                    return True
                if stack[6] == 3 and stack[0] in [5,6]:
                    return True
                if stack[0] == 4 and stack[6] == 4:
                    return True
        
        #Three faces -- outer faces only
        for i in range(len(stack)-2):
            a = str(stack[i])
            b = str(stack[i+1])
            c = str(stack[i+2])
            name = 'f'+a+'v'+b+'v'+c
            if name in ReducibleListF:
                return True
            name = 'f'+c+'v'+b+'v'+a
            if name in ReducibleListF:
                return True
        if full:
            a = str(stack[-2])
            b = str(stack[-1])
            c = str(stack[0])
            d = str(stack[1])
            name = 'f'+a+'v'+b+'v'+c
            if name in ReducibleListF:
                return True
            name = 'f'+c+'v'+b+'v'+a
            if name in ReducibleListF:
                return True
            name = 'f'+b+'v'+c+'v'+d
            if name in ReducibleListF:
                return True
            name = 'f'+d+'v'+c+'v'+b
            if name in ReducibleListF:
                return True
        
    #Four faces (outer faces only) ('f4v5v7v4' only)
    if len(stack) > 3:
        if full:
            s = str(stack[-3])+str(stack[-2])+str(stack[-1])
        else:
            s = ''
        for i in stack:
            s += str(i)
        if '4574' in s or '4754' in s:
            return True

    #Lastly, check for type 'c' configurations.
    if len(stack) < 2:
        return False
    c = str(central_length)
    #Start with using an outer face as a central face of a subconfiguration.  (Only forms cAaBC and cAaBCD are necessary to check here, which also means no 'x's.)
    if full:
        for i in range(-1,len(stack)-1):
            a = str(stack[i])
            b = str(stack[i-1])
            d = str(stack[i+1])
            li = ['c'+a+'a'+b+c,'c'+a+'a'+c+b,'c'+a+'a'+c+d,'c'+a+'a'+d+c,'c'+a+'a'+b+c+d,'c'+a+'a'+d+c+b]
            for name in li:
                if name in ReducibleListC:
                    return True
    else:
        if 'c'+str(stack[0])+'a'+str(stack[1])+c in ReducibleListC:
            return True
        if 'c'+str(stack[0])+'a'+c+str(stack[1]) in ReducibleListC:
            return True
        for i in range(1,len(stack)-1):
            a = str(stack[i])
            b = str(stack[i-1])
            d = str(stack[i+1])
            li = ['c'+a+'a'+b+c,'c'+a+'a'+c+b,'c'+a+'a'+c+d,'c'+a+'a'+d+c,'c'+a+'a'+b+c+d,'c'+a+'a'+d+c+b]
            for name in li:
                if name in ReducibleListC:
                    return True
        if 'c'+str(stack[-1])+'a'+str(stack[-2])+c in ReducibleListC:
            return True
        if 'c'+str(stack[-1])+'a'+c+str(stack[-2]) in ReducibleListC:
            return True
    #Now use the central face as the central face of the subconfiguration.  (For non-'x's -- which means two or three outer faces in the subconfiguration.)
    if full:
        for i in range(len(stack)):
            b = str(stack[i-2])
            a = str(stack[i-1])
            d = str(stack[i])
            li = ['c'+c+'a'+a+d,'c'+c+'a'+d+a,'c'+c+'a'+b+a+d,'c'+c+'a'+d+a+b]
            for name in li:
                if name in ReducibleListC:
                    return True
    else:
        if 'c'+c+'a'+str(stack[0])+str(stack[1]) in ReducibleListC:
            return True
        if 'c'+c+'a'+str(stack[1])+str(stack[0]) in ReducibleListC:
            return True
        for i in range(2,len(stack)):
            b = str(stack[i-2])
            a = str(stack[i-1])
            d = str(stack[i])
            li = ['c'+c+'a'+a+d,'c'+c+'a'+d+a,'c'+c+'a'+b+a+d,'c'+c+'a'+d+a+b]
            for name in li:
                if name in ReducibleListC:
                    return True
    #Now use the central face as the central face of the subconfiguration.  (For 'x's.  Each configuration manually coded.)
    if len(stack) < 4:
        return False
    stack_list = stack[:]
    for i in range(central_length-len(stack)):
        stack_list.append('x')
    stack_list.extend(stack[:])
    if central_length == 5:#c5a55x6
        spots5 = [i for i,x in enumerate(stack) if x==5]
        spots6 = [i for i,x in enumerate(stack) if x==6]
        for i in [x for x in spots5 if x<5]:
            if i+1 in spots5 and i+3 in spots6:
                return True
        for i in [x for x in spots6 if x<5]:
            if i+2 in spots5 and i+3 in spots5:
                return True
    elif central_length == 6:#c6a4xx4
        spots4 = [i for i,x in enumerate(stack) if x==4]
        for i in [x for x in spots4 if x<6]:
            if i+3 in spots4:
                return True
    elif central_length == 7:#c7a3xx4, c7a3xx5, c7a4xx4, c7a55x5
        spots3 = [i for i,x in enumerate(stack) if x==3]
        spots4 = [i for i,x in enumerate(stack) if x==4]
        spots5 = [i for i,x in enumerate(stack) if x==5]
        for i in [x for x in spots3 if x<7]:
            if i+3 in spots4 or i+3 in spots5:#3xx4, 3xx5
                return True
        for i in [x for x in spots4 if x<7]:
            if i+3 in spots3 or i+3 in spots4:#4xx3, 4xx4
                return True
        for i in [x for x in spots5 if x<7]:
            if i+3 in spots3:#5xx3
                return True
            if i+3 in spots5 and (i+1 in spots5 or i+2 in spots5):#55x5
                return True
    elif central_length == 8:#c8a3xx4, c8a3xxx4, c8a3x55x55, c8a4xx4, c8a4xxx4
        spots3 = [i for i,x in enumerate(stack) if x==3]
        spots4 = [i for i,x in enumerate(stack) if x==4]
        spots5 = [i for i,x in enumerate(stack) if x==5]
        for i in [x for x in spots3 if x<8]:
            if i+3 in spots4 or i+4 in spots4:#3xx4, 3xxx4
                return True
            if i+2 in spots5 and i+3 in spots5 and i+5 in spots5 and i+6 in spots5:#3x55x55
                return True
        for i in [x for x in spots4 if x<8]:
            if i+3 in spots3 or i+4 in spots3 or i+3 in spots4 or i+4 in spots4:#4xx3, 4xxx3, 4xx4, 4xxx4
                return True
        for i in [x for x in spots5 if x<8]:
            if i+6 in spots3 and i+1 in spots5 and i+3 in spots5 and i+4 in spots5:#55x55x3
                return True
    elif central_length == 9:#c9a3xx4, c9a3xxx4
        spots3 = [i for i,x in enumerate(stack) if x==3]
        spots4 = [i for i,x in enumerate(stack) if x==4]
        for i in [x for x in spots3 if x<9]:
            if i+3 in spots4 or i+4 in spots4:#3xx4, 3xxx4
                return True
        for i in [x for x in spots4 if x<9]:
            if i+3 in spots3 or i+4 in spots3:#4xx3, 4xxx3
                return True
    elif central_length == 10:#c10a3xxx4, c10a3xxxx3, c10a3xxxx4
        spots3 = [i for i,x in enumerate(stack) if x==3]
        spots4 = [i for i,x in enumerate(stack) if x==4]
        for i in [x for x in spots3 if x<10]:
            if i+4 in spots4 or i+5 in spots3 or i+5 in spots4:#3xxx4, 3xxxx3, 3xxxx4
                return True
        for i in [x for x in spots4 if x<10]:
            if i+4 in spots3 or i+5 in spots3:#4xxx3, 4xxxx3
                return True
    elif central_length == 11:#c11a3xxxx3
        spots3 = [i for i,x in enumerate(stack) if x==3]
        for i in [x for x in spots3 if x<11]:
            if i+5 in spots3:
                return True
    elif central_length == 12:#c12a3xxxxx3
        spots3 = [i for i,x in enumerate(stack) if x==3]
        for i in [x for x in spots3 if x<12]:
            if i+6 in spots3:
                return True
    
    #After all that, we've got nothing left.
    return False


#Takes a list of tuples (of lengths of outer faces) and removes duplicates (under cyclic rotation and reversal).
#Does not change in place -- creates new list from scratch and leaves input unaltered.
def simmer(input_list):
    keys = []
    new_list = []
    for t in input_list:
        s = ''
        for x in t:
            s += str(x)
        r = ''
        for x in reversed(t):
            r += str(x)
        
        dup = False
        for key in keys:
            if s in key or r in key:
                dup = True
                break
        
        if not dup:
            t_key = s+s[:-1]
            keys.append(t_key)
            new_list.append(t)
        
    return new_list




#Current Rules:
#Every 3-face takes 1 charge from each adjacent face.
#If a 4-face is adjacent to at least three (7+)-faces, then it takes 2/3 charge from each adjacent (7+)-face.
#If a 4-face is not adjacent to at least three (7+)-faces, then it takes 1 charge from each adjacent (9+)-face.
#If a 5-face is adjacent to at least three (7+)-faces, then it takes 1/3 charge from each adjacent (7+)-face.
#If a 5-face is not adjacent to at least three (7+)-faces, then it takes 1/2 charge from each adjacent (8+)-face.
def happy(central_length,outer_lengths):
    if central_length == 3:
        #Every 3-face takes 1 charge from each adjacent face.
        return True
    
    if central_length == 4:
        n7p = len([i for i,x in enumerate(outer_lengths) if x >= 7])
        #If a 4-face is adjacent to at least three (7+)-faces, then it takes 2/3 charge from each adjacent (7+)-face.
        if n7p >= 3:
            return True
        #If a 4-face is not adjacent to at least three (7+)-faces, then it takes 1 charge from each adjacent (9+)-face.
        else:
            n9p = len([i for i,x in enumerate(outer_lengths) if x >= 9])
            if n9p >= 2:
                return True
            else:
                return False
    
    if central_length == 5:
        n7p = len([i for i,x in enumerate(outer_lengths) if x >= 7])
        #If a 5-face is adjacent to at least three (7+)-faces, then it takes 1/3 charge from each adjacent (7+)-face.
        if n7p >= 3:
            return True
        #If a 5-face is not adjacent to at least three (7+)-faces, then it takes 1/2 charge from each adjacent (8+)-face.
        else:
            n8p = len([i for i,x in enumerate(outer_lengths) if x >= 8])
            if n8p >= 2:
                return True
            else:
                return False
    
    if central_length == 6:
        return True
    
    if central_length == 7:
        n3 = len([i for i,x in enumerate(outer_lengths) if x == 3])
        n4 = len([i for i,x in enumerate(outer_lengths) if x == 4])
        n5 = len([i for i,x in enumerate(outer_lengths) if x == 5])
        if 3*n3+2*n4+n5 <= 3:#Not an if and only if test.
            return True
        else:
            return False
    
    if central_length == 8:
        n3 = len([i for i,x in enumerate(outer_lengths) if x == 3])
        n4 = len([i for i,x in enumerate(outer_lengths) if x == 4])
        n5 = len([i for i,x in enumerate(outer_lengths) if x == 5])
        if 6*n3+4*n4+3*n5 <= 12:#Not an if and only if test.
            return True
        else:
            return False
    
    if central_length >= 9:
        n3 = len([i for i,x in enumerate(outer_lengths) if x == 3])
        n4 = len([i for i,x in enumerate(outer_lengths) if x == 4])
        n5 = len([i for i,x in enumerate(outer_lengths) if x == 5])
        if 2*n3+2*n4+n5 <= 2*(central_length-6):#Not an if and only if test.
            return True
        else:
            return False







##Checks if li_1 is majorized by li_2.  Assumes identical lengths and int entries.
#def maj(li_1,li_2):
    #for i in range(len(li_1)):
        #if li_1[i] > li_2:
            #return False
    #return True



#c5a55x6 = RedCon([(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),(11,12),(12,13),(2,7),(8,0),(11,0),(13,1)],'c5a55x6')




























































#Reducible configurations.
#---------------------------------------

#First entry is filler to shift indices to match with SGE task IDs.
StandardRedConList = [None,
'c3a3',
'cxa3x3',
'cxa3xx3',
'cxa3xxx3',
'cxa3xx4x3',
'cxa3xx5x3',
'cxa3xx6x3',
'cxa3xxx73',
'cxa3xx4xx3',
'cxa3xxx4x3',
'cxa3xxx5x3',
'cxa3xxx6x3',
'cxa3xxxx73',
'cxa3xxx4xx3',
'cxa3xxxx4xx3',
'c3a4',
'cxa3x4',
'cxa3x54',
'cxa3x64',
'cxa3x74',
'cxa37x4',
'cxa3xx54',
'cxa3xx64',
'cxa3xx74',
'cxa3x5x4',
'cxa3x6x4',
'cxa37xx4',
'cxa3xxx54',
'cxa3xxx64',
'cxa3xxx74',
'cxa3xx4x4',
'cxa3xx6x4',
'cxa3x55x4',
'cxa3x57x4',
'cxa3x65x4',
'cxa3x75x4',
'cxa3xx4xx4',
'c3a5',
'cxa375',
'cxa3775',
'cxa3xx45',
'cxa3x555',
'cxa3x575',
'cxa3x655',
'cxa3x675',
'cxa3x765',
'cxa3xxx45',
'c3a6',
'cxa376',
'cxa3x56',
'cxa3x66',
'cxa3776',
'cxa3xx46',
'cxa3xx66',
'cxa3x556',
'cxa3x576',
'cxa3x656',
'cxa3x676',
'cxa3xxx46',
'c4a4',
'cxa454',
'cxa464',
'cxa474',
'cxa4x54',
'cxa4x64',
'cxa4x74',
'cxa4xx54',
'cxa4x55x4',
'cxa4x56x4',
'cxa4x66x4',
'cxa45xx54',
'cxa45xxx54',
'cxa4x4x4x4',
'cxa455',
'cxa465',
'cxa4x45',
'cxa4575',
'cxa4755',
'cxa4775',
'cxa4x555',
'cxa45xxxx545',
'cxa456',
'cxa466',
'cxa476',
'cxa4x46',
'cxa4576',
'cxa4756',
'cxa5475',
'cxa5565',
'cxa555555',
'cxa546',
'cxa5476',
'cxa5746',
'cxa5666',
'c4a55',
'c4a56',
'c4a57',
'c4a66',
'c4a585',
'c4a676',
'c4a677',
'c4a686',
'c5a555',
'c5a556',
'c5a557',
'c5a565',
'c5a566',
'c5a567',
'c5a575',
'c5a576',
'c5a577',
'c5a656',
'c5a657',
'c5a666',
'c5a667',
'c5a676',
'c5a757',
'c3a777',
'c4a6787',
'c4a7777',
'c5a5x65',
'c5a5585',
'c5a5775',
'c5a5785',
'c5a5x66',
'c5a55x6',
'c5a56x6',
'c5a57x6',
'c5a5777',
'c5a5857',
'c5a6x66',
'c5a6776',
'c5a7577',
'c5a67777',
'c5a77777',
'c6a4xx4',
'c7a3xx4',
'c7a3xx5',
'c7a4xx4',
'c7a4x55',
'c7a4xx55',
'c7a4x5x5',
'c7a4x5xx5',
'c7a55x5',
'c8a3xx4',
'c8a3xxx4',
'c8a4xx4',
'c8a4xxx4',
'c8a3x55x55',
'c9a3xx4',
'c9a3xxx4',
'c9a545x5',
'c9a555x55x5',
'c10a3xxxx3',
'c10a3xxx4',
'c10a3xxxx4',
'c11a3xxxx3',
'c12a3xxxxx3']































































#Stuff for stem identifications.
#---------------------------------------




#Generates the sublists of size 3 and 2 which contain the first element of li and whose elements 
#are not separated by a single element from li.  Yields the sublist and the remaining contiguous 
#sublists.  Assumes that len(li) is at least 2.
def single_partition_generator(li,restrictions):
    n = len(li)
    
    #First, generate the lists of size 3.
    if n >= 5:
        a_options = [1,]
        a_options.extend(range(3,n-3))
        a_options.append(n-2)
    elif n==3:
        a_options = [1,]
    else:#n in [2,4]
        a_options = []
    
    for a in a_options:
        if n-1-a >= 3:
            b_options = [a+1,]
            b_options.extend(range(a+3,n-2))
            b_options.append(n-1)
        else:#n-1-a==1
            b_options = [n-1]
        
        for b in b_options:
            if {li[0],li[a]} in restrictions[1] and {li[0],li[b]} in restrictions[1] and {li[a],li[b]} in restrictions[1]:
                continue
            part = [li[x] for x in [0,a,b]]
            res = []
            if a > 1:
                res.append(li[1:a])
            if b-a > 1:
                res.append(li[a+1:b])
            if n-b > 1:
                res.append(li[b+1:])
            yield (part,res)
    
    #Next, generate the lists of size 2.
    if n >= 4:
        a_options = [1,]
        a_options.extend(range(3,n-2))
        a_options.append(n-1)
    elif n==2:
        a_options = [1,]
    else:#n ==3
        a_options = []
    
    for a in a_options:
        if {li[0],li[a]} in restrictions[0]:
            continue
        part = [li[x] for x in [0,a]]
        res = []
        if a > 1:
            res.append(li[1:a])
        if n-a > 1:
            res.append(li[a+1:])
        yield (part,res)
    
    #Send note of completion.
    yield None
    return

from collections import deque

#Update comments below for res_deque & res_stack.

#Given a list li, generates all partitions of li such that each part is of size 2 or 3 and such that
#the parts can be nested in order.
#Example: li = [1,2,3,4,5]
#[[1,2,3],[4,5]]
#[[1,2,5],[3,4]]
#[[1,4,5],[2,3]]
#[[1,2],[3,4,5]]
#[[1,5],[2,3,4]]
#'Nested' means that if the elements are put on a horizontal line in order, then a point for each part can 
#be plotted above the horizontal line so that if straight lines are drawn from the point to each of the part's 
#elements, then none of these lines intersect.
def full_partition_generator(li,restrictions):
    #At a given step of the generation, the currently determined parts of the partition are given by partition,
    #  the remaining contiguous pieces of li are given by res_stack[-1] (in order), and gen_stack[-1] gives 
    # the single-piece generator which generated partition[-1] and the first parts of res_stack[-1].
    gen_stack = [single_partition_generator(li,restrictions),]
    x = gen_stack[-1].next()
    if x == None:
        return
    partition = [x[0],]
    res_stack = [deque(x[1])]
    
    #Every pass through the loops takes us to the next partition to yield, then backtracks 
    # through any exhausted generators.
    while True:
#         print "\n"*4+"top:"+"---"*20
#         print 'partition:',partition
#         print 'res_stack:',res_stack
        good = True#boolean to indicate skipping to the backtrack stage -- necessary when single generators are empty from the start.
        #From wherever we are now, build all the way up to a full partition.
#         print "entering forward loop"
        while res_stack[-1]:
#             print "---"
            #If res_deque[0] is a restricted 2-list, then we need to backtrack now.
            #Otherwise, the single_partition_generator will yield at least once.
            #The only other possible issue would be if len(res_deque)==4 and every pair is restricted -- but
            #this cannot happen since K_4 is not half-plane-embeddable.
            next_interval = res_stack[-1].popleft()
#             if len(next_interval) < 3 and set(next_interval) in restrictions:
#                 good = False
#                 break
            gen_stack.append(single_partition_generator(next_interval,restrictions))
            x = gen_stack[-1].next()
            #If gen_stack[-1] is already exhausted, then we need to backtrack now.
            #This could happen if next_interval induces only restricted edges in the relevant spots.
            #Example:  len(next_interval)==2 and set(next_interval) in restrictions[0].
            #Example:  len(next_interval)==3 and each 2-set of next_interval in restrictions[0]|restrictions[1].
            if x == None:
                del gen_stack[-1]
                good = False
                break
            partition.append(x[0])
            res_stack.append(deque(res_stack[-1]))
            res_stack[-1].extendleft(x[1][::-1])
#             print 'partition:',partition
#             print 'res_stack:',res_stack
#         print "exiting forward loop"
        
#         print "good:",good
        #Yield
        if good:
            yield partition
        
#         print "entering backward loop"
        #Move back through all exhausted generators.  They are set to yield 'None' when finished.
        x = gen_stack[-1].next()#If good, the top generator must be exhausted since the last piece has size 2 or 3.
        while x == None:
#             print "---"
#             print "x:",x
            del gen_stack[-1]
            del partition[-1]
            del res_stack[-1]
#             print 'partition:',partition
#             print 'res_stack:',res_stack
            try:
                x = gen_stack[-1].next()
            except IndexError:#If gen_stack is empty, we are finished!
                return
        
#         print "exiting backward loop"
        
        #At this point, we have come back to the top-most generator and prodded it.  Adjust stacks accordingly.
        #Need to swap current top pieces for partition and res_stack with new things.
        del partition[-1]
        partition.append(x[0])
        del res_stack[-1]
        if len(res_stack) > 0:
            res_stack.append(deque(res_stack[-1]))
            res_stack[-1].extendleft(x[1][::-1])
        else:
            res_stack.append(deque(x[1]))


        

#Assumes a unique outer boundary for the configuration.
#
def get_roots_cyclic(rc_str):
    rc = makeRedCon(rc_str)
    G = Graph(rc.underlying_graph)
    
    #Calculate 2-vertices before the core is deleted
    twos = [x for x in G.vertices() if G.degree(x)==2]
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
#     print "DISTS:"
#     print dist1
#     print dist2
#     print dist12
    
    pre_threes = [x for x in G.vertices() if G.degree(x)==3]
    #If the configuration is type 'cxa', then the first and last vertices are roots that bridge the cyclic order.
    if rc_str[1]=='x':
        G.add_edge([0,len(rc.f)-1])
    #Delete the core so that we can get a cyclic ordering of the vertices (li).
    core = [x for x in pre_threes if not [y for y in G.neighbors(x) if G.degree(y) < 3]]
    for v in core:
        G.delete_vertex(v)
#     G.show()
#     print "core:",core
    G = G.hamiltonian_cycle()
    mi = min(twos)
    li = [mi,]
    next_one = G.neighbors(mi)[0]
    last_one = G.neighbors(mi)[1]
    G.delete_vertex(mi)
    while next_one != last_one:
        if next_one in twos:
            li.append(next_one)
        next_nbr = G.neighbors(next_one)[0]#There should be only one neighbor.
        G.delete_vertex(next_one)
        next_one = next_nbr
    if next_one in twos:
        li.append(next_one)
    return li,dist1,dist2
    

def identification_generator(li,dist1,dist2):
    dist12 = dist1[:]
    dist12.extend(dist2)
    
#     print "li:",li
    #Now iterate through the subsets to obtain all planar partitions on li into 2/3-lists for each subset.
    #Note:  If there is a 2-identification, it must not be in dist1.
    it = powerset(li)
    #Move past empty set and singletons.
    it.next()
    for x in range(len(li)):
        it.next()
    
    for x in it:
#         print "\n\nvertex subset:",x
        #If our subset only partitions into 2-sets of adjacent roots, then we don't feed it to the partition generator.
        if len(x) in [2,4]:
            empty_gen = True
            for edge in combinations(x, 2):
                if not set(x) in dist1:
                    empty_gen = False
                    break
            if empty_gen:
                continue
        for y in full_partition_generator(list(x),[dist1,dist12]):
#             print "\ny:",y
            #y is a particular partition of list(x) into 2- and 3-lists.
            #We want to reformat it to be a list of three lists:
            #0.  2-lists for stem-to-root identification,
            #1.  2-lists for stem-to-stem identification,
            #2.  3-lists for stem-to-stem identification.
            #The 3-lists of y go into (2).  The 2-lists of y may go into either (0) or (1) unless they are in dist2,
            #in which case they may go only into (0).
            L0_certain = []
            L1_potential = []
            L2 = []
            for z in y:
                if len(z) > 2:
                    L2.append(z)
                elif set(z) in dist2:
                    L0_certain.append(z)
                else:
                    L1_potential.append(z)
#             print "L0_certain:",L0_certain
#             print "L0_potential:",L1_potential
#             print "L2:",L2
            for subset in powerset(L1_potential):
#                 print ">>subset:",subset
                L0 = L0_certain[:]
                L1 = list(subset)
                L0.extend([pair for pair in L1_potential if not pair in subset])
                #Now, yield this particular choice for 2-identifications for this partition
#                 print [L0,L1,L2]
                yield [L0,L1,L2]
    
    return


from itertools import *

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))



#identifications is a list of lists:  
#   [stem-to-root 2-identifications, stem-to-stem 2-identifications, stem-to-stem 3-identifications].
#
def configuration_identifier(underlying_graph,identifications):
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
    
    
    
    
    
    
    
import choosability_8 as ch
import time

def check_stem_identifications(rc_str):
    print "-"*30+rc_str+"-"*30
    rc = makeRedCon(rc_str)
    begin = time.clock()
    count = 0
    good_count = 0
    bad_count = 0
    keys = ['error', 'greedy', 'CNS', 'brute']
    count_dict = {s:0 for s in keys}
    bad_li = []
    
    li,dist1,dist2 = get_roots_cyclic(rc_str)
    
    for idents in identification_generator(li,dist1,dist2):
        count += 1
    print "There are "+str(count)+" identifications to check.  (Took "+ch.timestring(time.clock()-begin)+" to count.)  Starting:"
    print "Index   Good   Bad   Time"
    
    i = 0
    for idents in identification_generator(li,dist1,dist2):
        i += 1
        G,f = configuration_identifier(rc.underlying_graph,idents)
        y = ch.fChoosableNoPrint(G,f,inprocess=4,print_mod=100,maxIS=True,rev=False)
        if y[0]:
            good_count += 1
        else:
            bad_count += 1
            bad_li.append(idents)
        count_dict[y[1]] += 1
        print " "*(len(str(count))-len(str(i)))+str(i)+" "*(len(str(count))-len(str(good_count)))+"   "+str(good_count)+"   "+" "*(len(str(count))-len(str(bad_count)))+str(bad_count)+"   "+ch.timestring(time.clock()-begin)
    
    print "\n"*7
    print "Total:",count
    print "Bad:",bad_count
    print "Good:",good_count
    if count == good_count:
        print "-->  Good for all stem identifications!"
    else:
        print "-->  Uh-oh!  Problem!"
    print
    for s in keys:
        print s+":",count_dict[s]
    end = time.clock()
    print "\nTime:",ch.timestring(end-begin)
    if good_count < count:
        print "Problems:"
        for x in bad_li:
            print x
    print
    if good_count == count:
        print "Complete:  "+rc_str+" is GOOD."
    else:
        print "Complete:  "+rc_str+" is BAD."

