

import choosability_four_planar as ch4
import time
#import copy

from sage.graphs.graph_input import from_graph6



#from subprocess import call


#def print_it(n,res,mod):
    #n_s = str(n)
    #res_mod_s = str(res)+"/"+str(mod)
    #count = 0
    #for x in call(["./plantri", "-c3m4", n_s, res_mod_s]):
        #print x
        #count += 1
    #print count







import sys

begin = time.clock()
count = 0
bad_count = 0
#CNS_count = 0
brute_count = 0
#bad_list = [] #Since chances are low we'll find a problem at all, we'll go ahead and risk storing any.
s = "."

for line in sys.stdin:
    count += 1
    G = Graph()
    from_graph6(G,line)#It looks like I don't need to strip a newline character.
    f = [4]*G.order()
    x = ch4.fChoosable_plantri_c3m4(G,f)
    if x[1] != 'CNS':
        brute_count += 1
        s += "brute"
        if not x[0]:
            s += "flag"
            bad_count += 1
            #bad_list.append((int(count),line[:]))
    if count %100 == 0:
        print ch4.timestring(time.clock()-begin)+" Total: "+str(count)+s
print "Done."
print "non-CNS:  "+str(brute_count)
print "bad:  "+str(bad_count)
print "#"+str(count)








#sage script:

#while:

#take from standard in

#convert to sage graph from graph6

#run fChoosable on it (adjust to not check first for 4-degeneracy?)

#print count, graph, yes/no, no_count, time








#bash script:  (make plantri on server)

# plantri -gc3m4 [n] [res]/[mod] | [sage path?] [sage script] >> plantric3m4fourchoosable[res]of[mod].out





