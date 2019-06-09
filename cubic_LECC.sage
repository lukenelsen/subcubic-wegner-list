import choosability as ch
import time
#import sage.all

import sys

k = int(sys.argv[1])
r = int(sys.argv[2])
m = int(sys.argv[3])

def LECC(G):
    max_degree = max(G.degree())
    LG = G.line_graph()
    LG.relabel()
    chi_prime = LG.chromatic_number()
    f = [chi_prime]*LG.order()
    return ch.fChoosable(LG,f,inprocess=4), max_degree, chi_prime


def cubic_LECC(k,res=0,mod=1):
    print "Check the List Edge Coloring Conjecture for all cubic connected graphs of order %d:"%(k)
    print "Res/Mod: "+str(res)+"/"+str(mod)
    gen = graphs.nauty_geng(str(k)+" -c -d3 -D3 "+str(res)+"/"+str(mod))
    total_count = 0
    bad_count = 0
    greedy_count = 0
    CNS_count = 0
    brute_count = 0
    greedy_time = 0.0
    CNS_time = 0.0
    brute_time = 0.0
    class1_count = 0
    class2_count = 0
    time_zero = time.clock()
    for g in gen:
        total_count += 1
        begin = time.clock()
        x = LECC(g)
        end = time.clock()
        if not x[0][0]:
            print "!!!!!!!!!!"
            print "#",total_count
            g.show()
            break
        if x[0][1]=='greedy':
            greedy_count += 1
            greedy_time += end-begin
        elif x[0][1]=='CNS':
            CNS_count += 1
            CNS_time += end-begin
        else:
            brute_count += 1
            brute_time += end-begin
        if x[1]<x[2]:
            class2_count += 1
        else:
            class1_count += 1
        if total_count%1000==0:
            print "        >>> %d checked, %s cumulative."%(total_count,ch.timestring(time.clock()-time_zero))
    time_end = time.clock()
    print "    Total Count:",total_count
    print "        Total Time:",ch.timestring(time_end-time_zero)
    print "    Bad Count:",bad_count
    print "    Greedy Count:",greedy_count
    print "        Greedy Time:",ch.timestring(greedy_time)
    print "    CNS Count:",CNS_count
    print "        CNS Time:",ch.timestring(CNS_time)
    print "    Brute Count:",brute_count
    print "        Brute Time:",ch.timestring(brute_time)
    print "    Class 1 Count:",class1_count
    print "    Class 2 Count:",class2_count



if __name__ == "__main__":
    cubic_LECC(k,r,m)
