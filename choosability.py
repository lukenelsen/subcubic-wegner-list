##Goals:
##    () Clean up.  (This includes consolidating inprocess CNS-pruning into one routine, which checks just the next subset with the current PLA.)

# To-Do:
    # After paper is finished, check references to Lemmas/Observations/Definitions.



### Is it true that if a connected component of an induced subgraph is G,f-reducible, then the entire induced subgraph is G,f-reducible?  Or is it that if *all* connected components of an induced subgraph are G,f-reducible, then the entire induced subgraph is G,f-reducible?


### order vertices such that first is in the fewest eligible (current, with multiplicity) vs first has the least sum of eligible second option over all first options
### when counting multiplicity above, should multiplicity be squared instead of just added?
### Shouldn't we be putting vertices with smaller f-values to be earlier as well?  Maybe not.  But it feels like this would be a factor *in addition to* how many colorability classes contain them.  (Perhaps this also depends on the LAscheme?)

### What to do when f-values are small (relative to the graph, but not relative to computer preprocessing) and clearly will fail?  We wanted to have something that does well at *verifying* choosability, but should we also have a quick check for failing choosability?  Can we try to greedily build a bad LA?


### Instead of an inprocess depth, maybe always try the first, say 100, IS-combinations?  Then we won't sink too much time into resistant colorability classes but might still trim out a few.

### Are denser or sparser subgraphs more likely to be CNS-resistant?  Would it at all help to know?

### Future:  Pre-check if a reduction results in disconnected graph.
### Future:  Choose vertex based on remaining subsets?

### *Can this code check all pairs of subgraphs?  If so, think about.  If not, make note about future possible change.



### Interesting phenomenon:  when I switched from doing RootedConnectedSubgraphs in initial_restriction_list_maker to the naive powerset + is_connected generation, there was not much difference in terms of the first-level generation (for the standard c7a4x5x5 square, a couple seconds both ways).  And as expected, there was not an impact on the list assignment building process since we weren't using RootedConnectedSubgraphs for that and the subsets all get reordered anyway.  However, it does appear to make an impact on the vetting process, even though we don't use RootedConnectedSubgraphs there, either.  The difference is the ordering we feed into the vetting process.  If we were directly checking each subgraph with the CNS, then it should not make a difference.  But we aren't, because we have a precheck to see if the net subgraph fits into our previous findings.  Here are the vetting times for the square of c7a4x5x5 with the {6,8} 2-stem (with Wegner f-vector):
#
###                           powerset      RootedConnectedSubgraphs
### with subgraph precheck    7m 55s (0)    0m 41s (22344)
### without the precheck      5m 59s        5m 56s
#
### With the RootedConnectedSubgraphs ordering, the precheck makes a big difference and save s us time by not checking the CNS.  But with the powerset ordering, it is actually WORSE to precheck.  Apparently, the powerset ordering doesn't take advantage of the precheck AT ALL because every subsequent subset is not a subset of any previous subset.  (The numbers in parentheses above indicated how many times the precheck resulted in a skip.  For reference, there are 22750 connected subgraphs in the example.)  So I made the generator rev_powerset to reverse the ordering so that for every subset, its subsets follow it in the ordering.
def rev_powerset(n):  ### reverse ordering of colex on sets.
    if n > 0:
        for subset in rev_powerset(n-1):
            subset.append(n-1)
            yield subset
        for subset in rev_powerset(n-1):
            yield subset
    else:
        yield []
    return
### Results:
###                           rev_powerset
### with subgraph precheck    0m 48s (22344)
### without the precheck      6m 58s
#
### Since I believe RootedConnectedSubgraphs was built on a similar idea as rev_powerset (with skipping disconnected subgraphs in mind), these results fit.




### Recalibrate num_prunes.
### Take out num_prunes?







#------------------------------------------------------
# Contents
#------------------------------------------------------

#    SubgraphsDealer
#    CNS
#    IndSets
#    restriction_list_single_check
#    set2tup
#    reduce
#    timestring
#    subsetstring
#    fChoosable
#    restrictedListAssignments
#    assignmentCheck



#------------------------------------------------------
# Call Structure
#------------------------------------------------------

# (Ignoring timestring, subsetstring, set2tup.)

#fChoosable
    #reduce
    #CNS_smart_print
    #restrictedListAssignments
        #rev_powerset
        #restriction_list_single_check
            #IndSets
            #reduce
            #CNS_smart
        #relabeling
        #SubgraphDealer
            #restriction_list_single_check
        #restriction_list_single_check
    #Gf_relabeling
    #assignmentCheck





#------------------------------------------------------
# Libraries, Comments, Global Variables
#------------------------------------------------------

# These scripts were written with the intent of using Python 2.7 before the switch to Python 3.  There should not be many differences to implement our code in Python 3 once this change is made.

#from sage.misc.misc import powerset
import copy
from sage.graphs.graph import Graph
import time
from sage.rings.polynomial.polynomial_ring_constructor import *
from sage.rings.finite_rings.finite_field_constructor import *
ZZ = sage.rings.integer_ring.IntegerRing()
#from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.graphs.independent_sets import IndependentSets
#from sage.combinat.subset import Subsets
from itertools import combinations

#from fractions import Fraction
from numpy import prod





#------------------------------------------------------
# The Routines
#------------------------------------------------------








### Go back through and add comments.

class SubgraphsDealer:
    
    def __init__(self,subgraphs_list,LAscheme="E"):
        self.subgraphs_list = copy.copy(subgraphs_list)
        self.LAscheme = LAscheme
        
        #subgraph is the current connected subgraph (set of vertices).
        self.subgraph = None
        
        self.required_vertex = None
    
    
    
    def copy(self):
        clone = SubgraphsDealer(copy.copy(self.subgraphs_list),self.LAscheme)
        clone.subgraph = copy.copy(self.subgraph)
        return clone
    
    
    
    def prune_subgraphs_list_with_E(self,E):
        i = 0
        while i < len(self.subgraphs_list):
            A = self.subgraphs_list[i]
            ###if A not subset E, then toss A from the list
            if not A <= set(E):
                del self.subgraphs_list[i]
                continue
            else:
                i += 1
    
    
    ### call if small condition met
    def prune_subgraphs_list_with_CNS(self,G_copy, ordering, f, CA):
        i = 0
        while i < len(self.subgraphs_list):
            A = self.subgraphs_list[i]
            ###if A can be be shown to be unnecessary via CNS, then toss A from the list
            if restriction_list_single_check(G_copy, ordering, f, CA, A):
                del self.subgraphs_list[i]
                continue
            else:
                i += 1
    
    
    
    
    ### E should be nonempty to call this
    def choose_required_vertex(self,L,E):
        
        if self.LAscheme == "E":
            self.required_vertex = E[0]
        
        elif self.LAscheme == "L":
            m = min([(L[y],y) for y in E])
            self.required_vertex = m[1]
    
    
    
    ### before calling this, call choose root.
    ### before calling this for the first time, call prune by E (and possibly by CNS)
    def next_with_required_vertex(self):
        ###go through list to find A containing required_vertex
        i = 0
        while i < len(self.subgraphs_list):
            ###if required_vertex not in A, then keep A in the list
            A = self.subgraphs_list[i]
            if not self.required_vertex in A:
                i += 1
                continue
            ###if A works, then set .subgraph to A, remove A, and return True
            else:
                self.subgraph = A
                del self.subgraphs_list[i]
                return True
        ###if no such A, then return False (v can't be filled)
        else:
            return False
    
    
    ### before calling this, call choose root
    def next_with_required_vertex_not_containing_set(self,D):
        while self.next_with_required_vertex():
            if D <= self.subgraph:
                continue  ### try again -- we can toss that subgraph
            else:
                return True
        else:
            return False  ### no subgraph remaining
    
    
    




















def CNS(G, f, verbosity=0):
    
    if verbosity > 7:
        print_flag = True
    else:
        print_flag = False
    
    #INITIALIZE
    
    #n is the number of vertices.
    n = G.order()
    if n != len(f):
        print "Error:  Length of f does not match order of G."
        return

    if [x for x in f if x<1]:
        if print_flag:  print "Error: f-vector contains non-positive entries."
        return False
    
    #gapleft is the remaining number of steps below the sum of (f-1)-values for our monomial.
    #gapleft begins as the disparity and should drop down to 0.
    #In general, it is given by disparity - sum([f[x]-1-deg[x] for x in range(v)]).
    gapleft = sum(f) - len(f) - G.size()
    if gapleft < 0:
        if print_flag:  print "    >> Disparity:  %d.  Cannot apply Nullstellensatz."%(gapleft)
        return False
    elif print_flag:  print "    >> Disparity:  %d.  Proceeding."%(gapleft)
    
    #Relabel graph so that expansion steps have smaller memory requirements.
    perminv = []
    ma = G.order()
    while len(perminv) < ma:
        degrees = [len(set(G.neighbors(x))-set(perminv)) for x in G.vertices()]
        for j in perminv:
            degrees[j] = ma
        mi = min(degrees)
        i = degrees.index(mi)
        perminv.append(i)
    
    if print_flag:  print "    >> Relabeling via",perminv
    #Want perminv[i] to be relabeled as i.
    perm = [perminv.index(x) for x in range(len(perminv))]
    G.relabel(perm)
    f = [f[x] for x in perminv]
    
    #Polynomial stuff.
    stringlist = ['x'+str(x) for x in range(n)]
    P=PolynomialRing(ZZ, n, stringlist)
    varlist = P.objgens()[1]
    
    #v is the vertex/variable we are inspecting.
    v = 0
    
    #deg, mins, and coeff all act as stacks.
    deg = []    #deg is a list holding the degrees we're currently using.
    mins = []   #mins is a list holding the lower bounds for the degrees we can use in deg.
    coeff = []  #coeff is a list holding the coefficients we've found for varlist[i]**deg[i].
    #Entries in deg start high, but go down by 1 to the corresponding mins value as we backtrack.
    
    #The coeff stack is always one longer than deg and mins; the first entry should be the polynomial 1.
    coeff.append(varlist[0]**0)
    
    
    #START LOOP
    
    while True:
        
        if coeff[-1] == 0 or gapleft < 0:  #If we have already encountered an issue with expansion.
            
            #BACKTRACK
            #If the current monomial expansion has already resulted in a coefficient of 0,
            #then the coefficient is 0 in the entire polynomial and we need to backtrack.
            #Lower the previous variable's degree to see if we get higher degrees afterward.
            
            while True:
                
                del coeff[-1]
                v -= 1
                
                #If the current degree is more than the minimum, we can break out of the loop to drop down now.
                if deg[-1] > mins[-1]:
                    break
                #Otherwise, march backwards.  We'll keep going until we can drop without going below mins.
                
                gapleft += f[v]-1-deg[-1]
                del deg[-1]
                del mins[-1]
                
                #If we've marched back through everything, then we've exhausted our options.
                if v < 1:
                    if print_flag:  print "    >> The CNS is inconclusive; no relevant monomials have nonzero coefficient."
                    G.relabel(perminv)
                    f = [f[x] for x in perm]
                    return False
            
            #Drop down the degree (and note this in diff).
            deg[-1] -= 1
            gapleft -= 1
            #Now that we've backtracked, we try moving forward again by finding the next coefficient piece.
        
        
        else:  #If there is hope for expanding the current monomial.
            
            #FORWARD
            
            #If we have built the entire monomial, then stop.
            if v > n-1:
                if print_flag:  print "    >> The graph is choosable!  \n    >>Certificate:  The monomial with power list "+str(deg)+" has coefficient "+str(coeff[-1])+"."
                G.relabel(perminv)
                f = [f[x] for x in perm]
                return True
            
            #Otherwise, append next values to deg and mins.
            
            #f[v]-1-gapleft is a lower bound for deg[v].  (Any power less than this number 
            #will cause the monomial to already have too low of total degree, even if the 
            #remaining powers are all the f-1 values.)
            mins.append(max(0,f[v]-1-gapleft))
            
            a = coeff[-1].degree(varlist[v])  #a is the degree of xv in the residual coefficient.
            b = len([x for x in G.neighbors(v) if v < x])  #b is the number of forward neighbors of v.
            
            #The monomial we have built so far cannot have a power of xv higher than a+b.
            if a+b < f[v]-1:
                deg.append(a+b)
                gapleft -= f[v]-1-a-b
            else:
                deg.append(f[v]-1)
            #Now find the next coefficient piece and try to keep going.
        
        
        #FIND THE NEXT COEFFICIENT PIECE
        
        #Expand remaining polynomial in the vth variable (times the previous coefficient, which might contain 
        #the vth variable) to find coefficient of varlist[:v+1]**deg[:v+1].
        
        prod = varlist[0]**0
        for factor in [varlist[v]-varlist[u] for u in G.neighbors(v) if v < u]:
            prod *= factor
        
        newco = 0
        for x in range(0,deg[-1]+1):
            newco += coeff[-1].coefficient({varlist[v]:x})*prod.coefficient({varlist[v]:(deg[-1]-x)})
        coeff.append(newco)
        
        
        #Increment and go back to the top of the loop!
        v += 1
        























### ????????

def IndSets(graph,maximal=True,complement=False):
    for I in IndependentSets(graph,maximal=maximal,complement=complement):
        yield I











#Checks whether the subset A is CNS-reducible when adding it to the current list assignment CA.  Returns True if the subset is reducible for the given CA; returns False otherwise (which means it must stay in the restriction list).
def restriction_list_single_check(G_copy, ordering, f, CA, A):
    #we test colorings from CA.append(A) by using backtracking with a stack of color classes
    ISgens = []      #elements of ISgens will be IndependentSets iterators
    coloring = []    #elements of coloring will be pairs (maximal independent set, union of all first entries up to this pair)
    c = 0            #at the top of the following loop, c is the color we are about to assign a color class (an IndependentSets item)
    while True:
        if c==len(CA):
            ordering_template = list(ordering)
            newf_template = list(f)
            for C in CA:
                for v in C:
                    newf_template[v] -= 1
            subG_template = Graph(G_copy)
            S = A    #S is going to be A, but without vertices which have already been colored
            if c>0:
                S = S - coloring[c-1][1]
                for v in coloring[c-1][1]:
                    del newf_template[ordering_template.index(v)]
                    del ordering_template[ordering_template.index(v)]
                    subG_template.delete_vertex(v)
            H = Graph(G_copy)
            for v in H.vertices():
                if v not in S:
                    H.delete_vertex(v)
            disparity_temp = sum(newf_template) - len(newf_template) - subG_template.size()
            #even if the disparity is negative, coloring certain independent sets in certain classes might still make the CNS applicable to the remaining graph.
            for I in IndSets(H):#check if the coloring given by 'coloring' and I shows via CNS that A is reducible
                if sum([newf_template[ordering_template.index(v)]-1 for v in I]) + len(S)-len(I) - sum([subG_template.degree(v) for v in I]) <= disparity_temp:
                    newf = list(newf_template)
                    newordering = list(ordering_template)
                    subG = Graph(subG_template)
                    for v in I:
                        del newf[newordering.index(v)]
                        del newordering[newordering.index(v)]
                        subG.delete_vertex(v)
                    for v in S-set(I):
                        newf[newordering.index(v)] -= 1
                    subG.relabel()
                    newf = [newf[newordering.index(x)] for x in range(G_copy.order()) if x in newordering]

                    subG,newf,same,low = reduce(subG,newf)
                    
                    if low:#In this case, the remaining subG,newf are NOT f-choosable.  Try another I.
                        continue

                    if len(newf) < 1 or CNS(subG,newf):
                        return (set(I),A)

                     #for fillvec in IntegerListsLex(sum(newf)-len(newf) - subG.size(),length=len(newf),ceiling=[newf[x]-1 for x in range(len(newf))]):
                        #tvec = [newf[x]-1-fillvec[x] for x in range(len(newf))]
                        #if not tvec or coeff_extract(subG,tvec):#If tvec is [], then we've actually just colored all of G.
                            #return True
            #otherwise, backtrack as far as necessary and take one step forward
            while c>0:
                try:
                    del coloring[c-1]
                    new = set(ISgens[c-1].next())#this throws the exception if the generator is exhausted
                    if c==1:
                        coloring.append((new,new))
                    else:
                        coloring.append((new,new|coloring[c-2][1]))
                    break
                except StopIteration:
                    del ISgens[c-1]
                    c -= 1
            if c==0:
                break
            else:
                continue

        #if c < len(CA), we come here and add a fresh generator to our stack
        if c==0:
            S = CA[0]
        else:
            S = CA[c] - coloring[c-1][1]
        H = Graph(G_copy)
        for v in H.vertices():
            if not v in S:
                H.delete_vertex(v)
        ISgens.append(IndSets(H,maximal=True))#Using maximal works since they maximal with respect to the *remaining* graph
        new = set(ISgens[c].next())
        if c==0:
            coloring.append((new,new))
        else:
            coloring.append((new,new|coloring[c-1][1]))
        c+=1

    return False
























#To handle the unorderedness of sets.  Takes a set S of nonnegative integers and returns a tuple of the same elements, ordered by <.  Could easily be adjusted to respect ordering by passing ordering in.
def set2tup(S):
    return tuple(x for x in range(max(S)+1) if x in S)






#Takes G,f as input and returns the pruned version, having removing any vertices with f-values higher than degree (iteratively) and having removed precolored vertices (f-values of 1) and adjusted the neighborhoods accordingly (also iteratively).
def reduce(G_copy,f_copy):
    G = Graph(G_copy)
    f = f_copy[:]
    redundance = [x for x in G.vertices() if f[x] > G.degree(x)]
    while redundance:
        for v in redundance:
            G.delete_vertex(v)
        G.relabel()
        for v in reversed(redundance):
            del f[v]
        redundance = [x for x in G.vertices() if f[x] > G.degree(x)]
    ones = [x for x in G.vertices() if f[x]==1]
    while ones:
        v = ones[0]
        for u in G.neighbors(v):
            f[u] -= 1
        del f[v]
        G.delete_vertex(v)
        G.relabel()
        ones = [x for x in G.vertices() if f[x]==1]
    same = (len(f_copy)==len(f)) #same indicates whether the returned configuration is the same as the input configuration
    low = (len([x for x in f if x<1])>0) #low indicates whether there are nonpositive f-values in the output configuration
    ### move same and low to checks outside of reduce
    return G,f,same,low









def timestring(time,setting="m"):
    t = int(time)
    d = str(int((10*(time-t))%10))
    if setting=="s":
        string = str(t)+"."+d+"s"
    else:
        s = str(t%60)
        m = int(t/60)
        if setting=="m":
            string = str(m)+"m "+s+"."+d+"s"
        else:
            h = int(m/60)
            m %= 60
            if setting=="h":
                string = str(h)+"h "+str(m)+"m "+s+"."+d+"s"
            else:
                d = int(h/24)
                h %= 24
                string = str(d)+"d "+str(h)+"h "+str(m)+"m "+s+"."+d+"s"
    return string







def subsetstring(A,n):
    string = " "*8+"["
    for i in range(n):
        string += " "*3
        if i in A:
            string += str(i)
        else:
            string += " "*len(str(i))
    string += "   ]"
    return string
















#Identical to fChoosableReturnBad, but does not store bad list assignments.  No bad_nums, z.
def fChoosable(G_copy, f_copy, inprocess=3, print_mod=1000, LAscheme="E", rev=False, verbosity=0, search_all=0):
    
    if not LAscheme in ["L","E"]:
        print 'Error:  LAscheme must be one of "L" or "E".'
        return

    if verbosity>0:
        start = time.clock()
        num_pruned = 0

    if len(f_copy) != G_copy.order():
        print "Error: f-vector length does not match graph order."
        return

    if [x for x in f_copy if x<1]:###use min(f_copy) instead?
        print "Error: f-vector contains non-positive entries.  (If correct f, graph is not f-choosable.)"
        return

    G,f,same,low = reduce(G_copy,f_copy)
    
    if not same:
        if verbosity>0:  print "The input configuration was reduced for having f-values either higher than the degree or equal to 1."
        if not f:
            if verbosity>0:  print "In fact, it reduced entirely; the graph is f-choosable via a greedy argument."
            return (True,'greedy')
        if verbosity>0:  print "f:",f
    if low:
        if verbosity>0:  print "Reduced f-vector contains non-positive entries.  Graph is not f-choosable."
        return (False,'precoloring')  ###different term?  pass the ordering?
    
    if verbosity>0:  print "Testing:  Is the graph (%d vertices, %d edges) "%(G.order(),G.size())+str(f)+"-choosable? \n\nFirst, we consider whether the Combinatorial Nullstellensatz applies."
    if verbosity>0:  boo = CNS(G,f,verbosity=10)
    else:  boo = CNS(G,f,verbosity=0)
    if verbosity>0:  print "    (Time so far:  "+timestring(time.clock()-start)+")"
    if boo:
        return (True,'CNS')  ### pass the monomial and coefficient?

    if verbosity>0:  start2 = time.clock()

    if verbosity>0:  print "\nNow we move on to our exhaustive search.  \nTo start up our process, we find and remove connected subgraphs which cannot be part of a bad list assignment, and then check the multiplicities of those that remain."

    it = restrictedListAssignments(G,f,inprocess_depth=inprocess,LAscheme=LAscheme,rev=rev,verbosity=verbosity)
    
    bad_list_assignments = []

    if verbosity>0:
        nodes_count = 0
        full_count = 0
        bottom_count = 0
        node_code = [0,]
        bad_LA_count = 0

    try:
        x = it.next()           #The first call to a Python iterator must be .next().
        #^^Recall: x is a tuple (CA,E,L,num_subsets,back,num_pruned), where CA is a colorability class assignment, E is the list of vertices whose lists aren't full, L is the list of the sizes of the current lists of the vertices, num_subsets is the list of the number of remaining subsets at each level of the subtree, and back is a boolean indicating whether we just got to this node after backtracking.  num_pruned is the sum of the number of subsets which were pruned inprocess at each step.  perminv is the reordering of vertices -- we use it only once to relabel the graph here.
        
        #G,f = Gf_relabeling(G,f,x[6])
        #perminv = x[6]
        ##Want perminv[i] to be relabeled as i.
        #perm = [perminv.index(z) for z in range(len(perminv))]
        #G_new = Graph(G)  ### why are we copying?
        #G_new.relabel(perm)
        #f_new = [f[x] for x in perminv]  
        
        if verbosity>0:
            print "    >> Start-up complete.  Total time so far:  "+timestring(time.clock()-start)
        
        #return
        
            print "\nNow, begin building and checking list assignments with the in-process pruning parameter set to %d. \nInformation will print every %d nodes we visit in the search space."%(inprocess,print_mod)
            print "( a | b ) means that colorability class is the a-th of b subgraphs on that level of the stack."
            start3 = time.clock()
            lap = start3

            nodes_count += 1
            node_code[0] += 1

        #At the top of this loop, we have obtained an x from the listAssignment iterator.  We check its representability -- unless it is representable, uncolorable, and full, we take the next x from the iterator.
        while True:

            if not x[1][-1] and verbosity>0:
                full_count += 1

            #If there is an unfilled vertex which has no supergraph candidates remaining, then turn back!
            #(If the list assignment is full, then this cannot happen.)
            ###flag = True
            #print "E:",x[1][-1]
            #print "->",[sum([len(y) for y in x[3][-1].containment_dict[v]]) for v in x[1][-1]]
            #The top subgraph has already been removed from the containment_dict.  If it can be used again, then we need not check its vertices.  Otherwise, we check them.
            ### stuff below used:
            ### ticket = yield (colClasses,E_stack,L,generator_stack,back,num_pruned,perminv,current_mult,max_mult)
            ###begin
            #if x[7][set2tup(x[3][-1].subgraph)] < x[8][set2tup(x[3][-1].subgraph)]:
                #for v in x[1][-1]:
                    #if not v in x[3][-1].subgraph and sum([len(y) for y in x[3][-1].containment_dict[v]]) < 1:
                        #flag = False
                        #break
            #else:
                #for v in x[1][-1]:
                    #if sum([len(y) for y in x[3][-1].containment_dict[v]]) < 1:
                        #flag = False
                        #break
            ###end
            
            ###if flag:
            c = assignmentCheck(G,x[0])

            #If the list assignment is not colorable and full, then report it and stop.
            if not c and not x[1][-1]:
                if verbosity>0:
                    bad_LA_count += 1
                    print " "*10+"*"*10+" BAD LIST ASSIGNMENT #%d:"%(bad_LA_count)+"*"*10
                    for i in range(len(x[0])):
                        print " "*10+subsetstring(x[0][i],G.order()),"(",node_code[i],")"#|",x[3][i],")"
                    print
                if search_all>0:
                    #Then switch c to be nonempty so that the iterator backtracks.
                    c = [0,]
                    if search_all>1:
                        bad_list_assignments.append(copy.copy(x[0]))
                else:
                    bad_list_assignments.append(copy.copy(x[0]))
                    break
            ###else:
                ###c = [0,]


            if verbosity>0 and x[3]:
                bottom_count += 1
            
            #Tell the listAssignment iterator to go to the next relevant list assignment (either moving forward or backtracking).
            x = it.send(c)
            
            if verbosity>0:
                if len(x[0]) > len(node_code):
                    node_code.append(1)
                else:
                    node_code[len(x[0])-1] += 1
                    for i in range(len(x[0]),len(node_code)):
                        del node_code[-1]

                num_pruned += x[4]

                nodes_count += 1
                if nodes_count%print_mod == 0:
                    print "-----> ",nodes_count,"partial list assignments inspected so far.  # full LAs:",full_count,".  # bottom LAs:",bottom_count,".  # bad LAs:",bad_LA_count,".  # inprocess prunes:",num_pruned,"."
                    print "^^^^^>  Time taken on this batch of PLAs:",timestring(time.clock()-lap),"    Time taken on the whole (subgraph) job so far:",timestring(time.clock()-start)
                    print "^^^^^>  Current partial list assignment:"
                    print subsetstring(x[0][0],G.order()),"(",node_code[0],"|",len(x[2][0].subgraphs_list),")"
                    for i in range(1,len(x[0])):
                        #print subsetstring(x[0][i],G.order()),"(",node_code[i],"|",x[3][i],")"
                        print subsetstring(x[0][i],G.order()),"(",node_code[i],"|",len(x[2][i].subgraphs_list),")"
                    print "\n"*3
                    lap = time.clock()
                    #print "E_stack:",x[1]
                    #for g in x[3]:
                        #print g.generator_cascade
                        #print

    except StopIteration:
        if verbosity>0:  print "\nNo bad list assignments!\n"

    if verbosity>0:
        if nodes_count < 1:#If CNS inconclusive but a vertex is not contained by any remaining subgraphs, then the normal start3 gets skipped.
            start3 = time.clock()
        end = time.clock()
        print "Finished.  %d PLAs visited which weren't pruned via CNS methods."%(nodes_count)
        print "# full LAs:  %d.  # bottom LAs:  %d.  # bad LAs:  %d.  # inprocess prunes:  %d"%(full_count,bottom_count,bad_LA_count,num_pruned)
        print "Time taken on initial CNS:       "+timestring(start2-start)
        print "Time taken on subgraph pruning:  "+timestring(start3-start2)
        print "Time taken on exhaustive search: "+timestring(end-start3)
        print "Total time taken on entire job:  "+timestring(end-start)
    if not bad_list_assignments:
        return (True,'brute')
    elif search_all>1:  ### stored all bad list assignments
        return (False,'brute',bad_list_assignments)
    elif search_all<1:  ###found one bad list assignment and then stopped
        return (False,'brute',bad_list_assignments[0])
    else:  ###printed all bad list assignments
        return (False,'brute')





### order so that least vertex indexed in minimizer of f(v)*(# subsets containing v -- with multiplicity)?  or (...)**f(v)?  or just by min f, min (#...)?








### For commenting:  make it clear that E_stack is appended with the E that is used for the next colClasses subset

### Make sure that E_stack and prune with E is updated immediately after adding a new colorability class.

### Before checking for a coloring, check for an unfillable vertex and keep track of the firsts.


#inprocess_depth controls how many partial list assignments are checked for the reducibility when extending during the normal routine.
def restrictedListAssignments(G,f,inprocess_depth=3,LAscheme="E",rev=False,verbosity=0):
    #moved from parameters
    ordering = range(G.order())
    reducers=[]

    #Initialize the things!
    generator_stack = []
    colClasses = []                                    #colClasses[i] is {vertices with color i in their list}.  We yield this.
    L = [0]*G.order()                                  #Keeps track of how many colors each vertex is currently assigned.
    D = None
    E_stack = []
    E_stack.append([x for x in ordering if L[x]<f[x]])       #Eligible vertices for colorability classes.
    col = 0                                            #Keeps track of next colorability class('s color).
    ticket = None

    reducers = [set(x) for x in reducers]
    
    
    
    CA = []#CA is from the inprocess version of this routine.  Here, it can be removed entirely.
    
    if verbosity>1:
        begin = time.clock()
        
        ###used to be stuff between here.
        
        lap = time.clock()
        
        print "    >> Now vetting all connected subgraphs."

    initial_restriction_list=[]
    
    reducibleList = []
    
    if verbosity>1:
        tenpercent = int(2**G.order()/100)+1###onepercent?
        counter = 0
        lap = time.clock()
        begin = lap###
    
    it = rev_powerset(G.order())
    A = set(it.next())    
    while A:  ### set up this way to avoid the empty set.
    
        if not G.subgraph(A).is_connected():
            A = set(it.next())
            continue
    
        if verbosity>1:
            if counter %tenpercent==tenpercent-1:
                print "        %d vetted.  Time on batch:  "%(counter)+timestring(time.clock()-lap)
                lap = time.clock()
            counter += 1

        red = False

        #if A is reducible for CA, continue/break to move on to next subset

        #first see if A is reducible by previous observations
        for x in reducibleList:
            if x[0] <= A and A <= x[1]:
                red = True
                break
        if red:
            A = set(it.next())
            continue

        #now see if A can be shown reducible by CNS
        red = restriction_list_single_check(G, ordering, f, CA, A)  ### red is either False or a tuple (I,A) where I is the independent set that worked with the CNS to show that A is "reducible" with CA

        if not red:
            initial_restriction_list.append(A)  ### start changing into one-dimensional list here.
        else:
            reducibleList.append(red)
        
        A = set(it.next())

    if verbosity>1:
        print "        Time vetting all:  "+timestring(time.clock()-begin)

        #num_subsets is going to keep track of how many nodes we have at each level of the subtree.  Entries are overwritten every time we go down to start the level.
        #num_subsets = [sum([len(x) for x in initial_restriction_list_dict[()]]),]
        #for i in range(sum(f)):
            #num_subsets.append(0)

        #print "    >>",num_subsets[0],"remaining subgraphs can potentially be used to build list assignments."
        print "    >>",len(initial_restriction_list),"remaining subgraphs can potentially be used to build list assignments."
        
        if len(initial_restriction_list) == 0:
            print "    >> So we can't build any bad list assignments!"
            return

        tenpercent = int(len(initial_restriction_list)/100)+1#int(num_subsets[0]/10)+1
        lap = time.clock()
        counter = 1
        begin = lap
        
    num_pruned = 0  ### need to define num_pruned regardless of verbosity since it is passed in ticket.

    max_mult_init = {}
    #current_mult = {}
    for A in initial_restriction_list:
        
        CA = [A,]

        cap = min([f[x] for x in A])#used to be min(min([f[x] for x in A]),len(A)), but I can't remember why -- I think we thought the connected small pot lemma was true.  Since f-values are at most 7 for wegner configurations anyway, this should be okay.  It might be worth considering putting in a ceiling (less than cap) of how far we check for large subgraphs with high f-values since this could take a lot of time that might not be necessary.  I mean, what's the practical difference between a max multiplicity of 5 and a max multiplicity of 10?
        while len(CA) < cap:

            if restriction_list_single_check(G, ordering, f, CA, A):
                break

            CA.append(A)
        max_mult_init[set2tup(A)] = len(CA)
        
        if verbosity>1:
            #current_mult[set2tup(A)] = 0
            if counter%tenpercent==0:
                print "       ",counter,"multiplicities checked.  Time on batch:  "+timestring(time.clock()-lap)
                lap = time.clock()
            counter += 1


    if verbosity>1:
        print "        Time for all multiplicity checks:  "+timestring(time.clock()-begin)
        #print max_mult_init.values()
        print "    >> Multiplicities breakdown:"
        for i in range(1,max(max_mult_init.values())+1):
            print "        Number of subgraphs with multiplicity %d:  %d"%(i,len([x for x in max_mult_init.values() if x==i]))
            #if i==max(max_mult_init.values()):
                #print [x for x in max_mult_init.keys() if max_mult_init[x]==i]

        print "    >> Vertex containment breakdown for these leftover subgraphs:"
        for v in range(G.order()):
            print "        Number containing vertex %d: "%(v),len([x for x in initial_restriction_list if v in x])

    #interim_cascade = []
    subgraph_pool = {frozenset(x) for x in initial_restriction_list}
    perminv = []
    ma = sum(max_mult_init.values())+1
    while len(perminv) < G.order():
        containments = [sum([max_mult_init[set2tup(x)] for x in subgraph_pool if v in x]) for v in G.vertices()]
        for j in perminv:
            containments[j] = ma
        mi = min(containments)
        i = containments.index(mi)
        perminv.append(i)
        batch = [x for x in subgraph_pool if i in x]
        #interim_cascade.append(batch)
        subgraph_pool -= set(batch)
    
    if verbosity>1:  print "Relabeling everything using perminv:",perminv
    #Want perminv[i] to be relabeled as i###, and want j to be relabeled as perm[j].
    perm = [perminv.index(x) for x in range(len(perminv))]
    G.relabel(perm)
    f = [f[x] for x in perminv]
    relabeled_subgraph_pool = {frozenset([perm[i] for i in x]) for x in initial_restriction_list}

    max_mult = {}
    current_mult = {}
    for A in relabeled_subgraph_pool:  ### A is new, but max_mult_init is for old labelings -- go back using perminv.
        max_mult[set2tup(A)] = max_mult_init[set2tup({perminv[a] for a in A})]
        current_mult[set2tup(A)] = 0

    ##No order
    #cascade = interim_cascade
    
    subgraphs_list = []
    mult = 1
    while relabeled_subgraph_pool:
        subpool = {x for x in relabeled_subgraph_pool if max_mult[set2tup(x)]==mult}
        relabeled_subgraph_pool -= subpool
        size = G.order()
        while subpool:
            batch = [x for x in subpool if len(x)==size]
            subpool -= set(batch)
            #print batch
            ###
            #batch_rooted = []
            #for i in range(G.order()):
                #batch_rooted.append([])
            #for x in batch:
                #batch_rooted[min(x)].append(set(x))
            #for sub_batch in batch_rooted:
                #subgraphs_list.extend(sub_batch)
            ###
            for x in batch[::-1]:  ### For some reason, this was slightly faster in my small battery of tests -- especially for 7:4x5x5 with {6,8} 2-stem.
                subgraphs_list.append(set(x))
            ###
            #for x in batch:
                #subgraphs_list.append(set(x))
            size -= 1
        mult += 1
    #return  ### future testing?  order vertices to have smallest f*#supersets?  Order subsets to be lex (or rev_lex, or whatever the middle thing is?) within the sub-batches?
    if rev:
        subgraphs_list = subgraphs_list[::-1]
    
    ###Multiplicities greatest to least, then lengths longest to shortest
    #if rev:
        #cascade = []
        #for Y in interim_cascade:
            #X = [(y,max_mult[set2tup(y)]) for y in Y]
            #newX = []
            #if X:
                #mi=max([x[1] for x in X])
            #while len(newX) < len(X):
                #batch = [(i,len(X[i][0])) for i in range(len(X)) if X[i][1]==mi]
                #if not batch:
                    #mi -= 1
                    #continue
                #ma = max([x[1] for x in batch])
                ##batch_sorted = []#sorts the indices so that sizes of subgraphs (with the same multiplicity) are decreasing
                #while ma > 0:#<= max([x[1] for x in batch]):
                    #mini_batch = [x[0] for x in batch if x[1]==ma]
                    ##print mini_batch
                    ##batch_sorted.extend(mini_batch)
                    #ma -= 1
                    #newX.extend([X[i][0] for i in mini_batch])
                #mi -= 1
            #cascade.append(newX)
    
    ##Multiplicities least to greatest, then lengths shortest to longest
    #else:
        #cascade = []
        #for Y in interim_cascade:
            #X = [(y,max_mult[set2tup(y)]) for y in Y]
            #newX = []
            #mi=1
            #while len(newX) < len(X):
                #batch = [(i,len(X[i][0])) for i in range(len(X)) if X[i][1]==mi]
                #if not batch:
                    #mi += 1
                    #continue
                #ma = min([x[1] for x in batch])
                ##batch_sorted = []#sorts the indices so that sizes of subgraphs (with the same multiplicity) are decreasing
                #while ma <= max([x[1] for x in batch]):
                    #mini_batch = [x[0] for x in batch if x[1]==ma]
                    ##print mini_batch
                    ##batch_sorted.extend(mini_batch)
                    #ma += 1
                    #newX.extend([X[i][0] for i in mini_batch])
                #mi += 1
            #cascade.append(newX)
    
    
    #subgraphs_list = []
    #for casc in cascade:
        #subgraphs_list.extend(casc)
    
    
    it = SubgraphsDealer(subgraphs_list,LAscheme=LAscheme)
    
    generator_stack.append(it)                         #Stack of subgraph iterators; cth element generates the c-colorability class.


    it.choose_required_vertex(L,E_stack[0])  ### aka it.required_vertex = ordering[0]
    if it.next_with_required_vertex():
        A = generator_stack[col].subgraph  ### generator_stack[col] == it
    
        colClasses.append(A)
        current_mult[set2tup(A)] += 1
        for v in colClasses[col]:
            L[v]+=1
        E_stack.append([x for x in ordering if L[x]<f[x]])  ### for x in E
    else:
        return

    back = False

    #Begin generation!
    while generator_stack:

        #print colClasses
        #print "E_stack:",E_stack
        #for x in generator_stack[-1].generator_cascade:
            #print x
        #print
        #for x in generator_stack[-1].containment_dict.values():
            #print x
        #print "\n\n"
        
        #Yield the current assignment; receive information from caller.
        #ticket = yield (colClasses,E,L,num_subsets,back,num_pruned,perminv)
        ticket = yield (colClasses,E_stack,generator_stack,back,num_pruned)  ### maybe rename ticket to something like returned_coloring?
        #print "ticket:",ticket
        #ticket is a list, either a coloring for colClasses or else [] (uncolorable assignment).
        back = False
        
        #print "ticket:",ticket

        #Case 1:  Move forward.  (Build the current partial list assignment up one more level.)
        #Since it is a new color, we start/copy a new iterator and add it to the stack (rather than call one we already have).
        #(There is no use for D in this case.)
        if not ticket:#Case:  The partial list assignment was not colorable.
            #Add the new colorability class iterator to the iterator stack, then call the first subgraph, add it to the stack, and update.
            generator_stack.append(generator_stack[col].copy())  #Add new subgraph iterator (for the color col+1) to the stack.
            col+=1                                               #Increment col -- now the class we are generating is for the color col.

            #Because we've just made a copy of the top generator, we need to check the current colorability class for repetition:
            A = generator_stack[col].subgraph
            if current_mult[set2tup(A)] < max_mult[set2tup(A)] and A <= set(E_stack[-1]):
                if len(colClasses) < inprocess_depth:
                    if not restriction_list_single_check(G, ordering, f, colClasses, A):
                        colClasses.append(A)
                        current_mult[set2tup(A)] += 1
                        for v in colClasses[col]:                            #Update L.
                            L[v]+=1
                        E_stack.append([x for x in ordering if L[x]<f[x]])         #Update E.
                        #num_subsets[col] = sum([len(x) for x in generator_stack[col].generator_cascade])+1
                        continue
                else:#If we have more than a few colors and we know A fits, let's just try it!
                    colClasses.append(A)
                    current_mult[set2tup(A)] += 1
                    for v in colClasses[col]:                            #Update L.
                        L[v]+=1
                    E_stack.append([x for x in ordering if L[x]<f[x]])         #Update E.
                    #num_subsets[col] = sum([len(x) for x in generator_stack[col].generator_cascade])+1
                    continue
            #If the current subgraph doesn't fit into E, then we check for the next eligible subgraph after updating restriction lists.
            #We don't update the restriction lists when backtracking since they've already been updated here first.
            
            
            generator_stack[col].choose_required_vertex(L,E_stack[-1])
            generator_stack[col].prune_subgraphs_list_with_E(E_stack[-1])

            #(Update the restrictionList from the current colClasses if we are in the early stages.)
            if len(colClasses) < inprocess_depth:  ### change the condition
                generator_stack[col].prune_subgraphs_list_with_CNS(G, ordering, f, colClasses)
            
            if generator_stack[col].next_with_required_vertex():
                A = generator_stack[col].subgraph
                colClasses.append(A)
                current_mult[set2tup(A)] += 1
                for v in colClasses[col]:                            #Update L.
                    L[v]+=1
                E_stack.append([x for x in ordering if L[x]<f[x]])         #Update E.
                #num_subsets[col] = sum([len(x) for x in generator_stack[col].generator_cascade])+1
                continue
            #If there is no such subgraph, then we are as full as possible and need to backtrack.
            back = True

        #Case 2a:  Backtrack because we tried to fill the current list assignment, but there were no more valid subsets.
        #Skip over the following else portion down to the main backtracking.


        #Case 2b:  Backtrack because the current list assignment colClasses is colorable (on the whole graph).
        #Start by moving the top color forward with respect to not containing the current colorability class.
        else:
            back = True
            #Delete the top colorability class (but not the generator for the class quite yet).
            #print col, colClasses[col], set2tup(colClasses[col]), colClasses
            for v in colClasses[col]:
                L[v] -=1
            current_mult[set2tup(colClasses[col])] -= 1
            del colClasses[col]
            del E_stack[-1]#E = [x for x in ordering if L[x]<f[x]]

            #If the provided coloring (or lack of a provided coloring) is not on all the vertices, then we can't skip ahead in any clever way.
            if len(ticket) < G.order():
                generator_stack[col].choose_required_vertex(L,E_stack[-1])
                generator_stack[col].prune_subgraphs_list_with_E(E_stack[-1])  ### pretty sure this line does nothing
                if generator_stack[col].next_with_required_vertex():
                    A = generator_stack[col].subgraph
                    colClasses.append(A)
                    current_mult[set2tup(A)] += 1
                    for v in colClasses[col]:
                        L[v] +=1
                    E_stack.append([x for x in ordering if L[x]<f[x]])
                    continue
                #If not, then we need to peel off more generators/layers.

            #But if the provided coloring is on all vertices, then we can use the top color class to potentially skip past some colorability class candidates.
            else:
                D = [x for x in colClasses[-1] if ticket[x]==col]  #D is the set of vertices in the coloring's highest color class.

                #First, see if the generator can give us something new:
                generator_stack[col].choose_required_vertex(L,E_stack[-1])
                generator_stack[col].prune_subgraphs_list_with_E(E_stack[-1])  ### pretty sure this line does nothing
                if generator_stack[col].next_with_required_vertex_not_containing_set(D):
                    A = generator_stack[col].subgraph
                    colClasses.append(A)
                    current_mult[set2tup(A)] += 1
                    for v in colClasses[col]:
                        L[v] +=1
                    E_stack.append([x for x in ordering if L[x]<f[x]])
                    continue
                #If not, then we need to peel off more generators/layers.  Note: D no longer serves any purpose past here.


        #At the moment, generator_stack[col] is an iterator that has no valid subgraphs left (either because it was a false start in Case 2a or because we backtracked through it in Case 2b).  But colClasses[col] is undefined (either because it didn't get there in Case 2a or because we deleted it in Case 2b).
        #Case 2:  Backtrack.  Step back as many colors as necessary and then obtain a next subgraph -- add that to colClasses.

        #At this point, we are guaranteed to take at least one step back through our stacks.
        del generator_stack[col]
        col -= 1
        if col<0:
            return
        for v in colClasses[col]:
            L[v] -=1
        current_mult[set2tup(colClasses[col])] -= 1
        del colClasses[col]
        del E_stack[-1]#E = [x for x in ordering if L[x]<f[x]]

        #Now step back as many times as necessary.
        generator_stack[col].choose_required_vertex(L,E_stack[-1])
        while not generator_stack[col].next_with_required_vertex():

            del generator_stack[col]
            col -= 1
            if col<0:#If generator_stack is empty, then we have backtracked entirely.  A StopIteration error will be thrown.
                return
            #Prepare for bumping the now-top generator on the stack.
            for v in colClasses[col]:
                L[v] -=1
            current_mult[set2tup(colClasses[col])] -= 1
            del colClasses[col]
            del E_stack[-1]#E = [x for x in ordering if L[x]<f[x]]
            generator_stack[col].choose_required_vertex(L,E_stack[-1])

        #Now we have successfully bumped the generator to another subgraph.  Push the subgraph onto the colClass stack.
        A = generator_stack[col].subgraph
        colClasses.append(A)
        current_mult[set2tup(A)] += 1
        for v in colClasses[col]:
            L[v]+=1
        E_stack.append([x for x in ordering if L[x]<f[x]])
        continue















def assignmentCheck(G,CAtocopy):
    #Convert the colorability classes assignment to a list assignment.
    CA = []
    for a in CAtocopy:
        CA.append(copy.deepcopy(a))
    n = G.order()
    LA = []                     #LA is the list assignment, where L[v] is the set of colors assignable to vertex v.
    for j in range(n):          #It will change according to our partial colorings during the algorithm.
        LA.append(set())
    for col in range(len(CA)):
        for v in CA[col]:
            LA[v].add(col)

    #Initialize!
    i = 0          #Which vertex we're coloring at the moment.
    c = [-1]*n     #Our coloring function.  c[i] is the current color of vertex i (value -1 means no color assigned).
    stack = []     #Our stack.  Item i is the set of vtxs we deleted color c[i] from (that is, the forward-nbrs which had c[i] in their lists).

    #At the top of this while loop, we have a partial coloring c.  We have either just colored vertex i-1 (in which case c has colored vertices 0 through i-1 and all of i's remaining colors are available to try), or we have backtracked from a partial coloring that didn't work (in which case we had colored i with something before, and now we need to try another of the remaining colors).
    while i < n:   #Go as long as there's something to be colored; failure of all colorings will throw an exception (empty stack) inside.

        #What colors have we not yet tried on i?  We've already tried all permissible colors of label at most c[i].
        remainingColors = {x for x in LA[i] if x>c[i]}  #(The comparison x>c[i] is why we use -1 for vertices with no assigned color.)

        backtrack = False

        #Check if any of the other uncolored vertices have no remaining colors.
        for v in range(i+1,n):
            if not LA[v]:
                backtrack = True

        #Case 1:  We can assign a new color to vertex i and we haven't discovered that farther vertices are empty.  So we do that.
        if remainingColors and not backtrack:
            c[i] = min(remainingColors)
            N = {x for x in G.neighbors(i) if x>i and x in CA[c[i]]}   #No need to bother with back-neighbors, since they have a color assigned.
            for v in N:
                LA[v].remove(c[i])    #Remove the new color from all forward-nbrs of i.
                CA[c[i]].remove(v)
            stack.append(N)
            i+=1

        #Case 2:  There are no colors left for some uncolored vertex, so we pop the stack and backtrack.
        else:
            try:
                oldN = stack.pop()
                for v in oldN:
                    LA[v].add(c[i-1])      #Add the attempted color back to the appropriate vertex lists.
                    CA[c[i-1]].add(v)
                c[i] = -1                #Mark i as uncolored.
                i-=1
            except IndexError:       #If the stack was empty, we've exhausted all possibilities.  Return an empty coloring.
                return []

    return c     #At this point, i == n.







