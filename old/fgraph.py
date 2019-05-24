
# Package to import and export fgraph6 format into Sage

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


def leading_coefficient(G,f):
    n=G.num_verts()
    num_edges=G.num_edges()
    
    x=[var("x%d" % i) for i in range(n)]
    poly=prod((x[e[0]]-x[e[1]]) for e in G.edges())
    print "poly=",poly
    h=prod(x[i]^(f[i]-1) for i in range(n))
    print "h=",h
    coeff=poly.coefficient(h)  # this does not do what I think it should do
    print "coeff=",coeff
    print expand(poly)
    return coeff


if __name__=="__main__":
    # running directly from sage
    from sage.all_cmdline import *
    
    #print "mapping=",mapping
    #print "inverse_mapping=",inverse_mapping
    
    G=graphs.CycleGraph(5)
    G.add_vertices([5])
    G.add_edges([(3,5),(4,5)])
    #G.plot().save("temp.pdf")
    print "graph6:",G.graph6_string()
    f=[2]*G.num_verts()
    f[4]=3  # increase the list size
    s=encode_fgraph_string(G,f)
    print s
    G,f=decode_fgraph_string(s)
    print "f=",f
    print "graph6:",G.graph6_string()
    #print G.adjacency_matrix()
    
    coeff=leading_coefficient(G,f)
    print "coeff=",coeff
