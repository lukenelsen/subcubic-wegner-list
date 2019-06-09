##Goals:
##    () Separate start-up.
##    () Bring parallelization up to speed for non-Minimal stuff.
##    () Clean up.  (This includes consolidating inprocess CNS-pruning into one routine, which checks just the next subset with the current PLA.)


from sage.misc.misc import powerset
import copy
from sage.graphs.graph import Graph
import time
from sage.rings.polynomial.polynomial_ring_constructor import *
from sage.rings.finite_rings.finite_field_constructor import *
ZZ = sage.rings.integer_ring.IntegerRing()
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.graphs.independent_sets import IndependentSets
from sage.combinat.subset import Subsets
from itertools import combinations

from fractions import Fraction
from numpy import prod


#Takes a graph G, a set of 'eligible' neighbor vertices E of G, and a set S of vertices of G as input.  Returns the neighbors of S in G which are in E.  The neighborhood is open in that we don't require a vertex in S to be a neighbor of itself, but some vertices of S might be in the neighborhood.
def nbrhood(G,E,S):
    N = set()
    for v in S:
        for w in G.neighbors(v):
            if w in E:
                N.add(w)
    return N

#Takes a list/tuple A and a list/tuple/set sublist.  Returns the indicator vector of the intersection of sublist and A (as a list of 0s and 1s).  Example:  convert_to_vector([1,2,3,3,4,6,4,'a'],{343,2,4,5,17,'iamyourfather'}) >> [0,1,0,0,1,0,1,0].
def convert_to_vector(A,sublist): #sublist doesn't actually have to be a subset of A; this effectively operates on (sublist & A).
    vector = []
    for i in range(len(A)):
        if A[i] in sublist:
            vector.append(1)
        else:
            vector.append(0)
    return vector



#Subset generator class.  Takes a list given_list as input, and (after first is called) proceeds via next methods through subsets/siblists of given_list.  The ordering of the generation is the reverse of binary counting by the indicator vector.  Example:  given_list = [3,4,3] --> [3,4,3], [4,3], [3,3], [3], [3,4], [4], [3], [].  Example on the indicator vectors:  [1,1,1], [0,1,1], [1,0,1], [0,0,1], [1,1,0], [0,1,0], [1,0,0], [0,0,0].  The attribute x is the current indicator vector for the subset.  The method x_converted is the current subset.
#Note:  This and SubsetNotEmpty can be adjusted to use a sentinel-type digit instead of a "second" attribute.
class Subset:



    def __init__(self,given_list):
        self.given_list = given_list
        self.n=len(given_list)
        self.second = False         #A boolean indicating whether the first subset (all 1s) has been considered.
        self.x = None



    def copy(self):                         #Only the list x needs copied/deepcopied.
        clone = Subset(self.given_list)
        clone.n = self.n
        clone.second = self.second
        clone.x = copy.deepcopy(self.x)
        return clone


    #Perhaps this would just as well go in the constructor.
    def first(self):
        if self.n>0:
            self.x=[1]*self.n  # indicator vector for the subset; starts at all 1s.
            return True
        else:
            return False



    def x_converted(self):
        return [self.given_list[w] for w in range(self.n) if self.x[w]==1]


    #The next method moves x to the very next subset.
    def next(self):
        if not self.second:
            self.second = True
            return True
        for i in range(self.n):
            if self.x[i]==1:  # i is the leftmost 1
                self.x[i]=0
                i-=1
                while i>=0:
                    self.x[i]=1
                    i-=1
                return True
        else:  # there are no 1s left, so we are done
            return False



    #The next_completely_avoiding method moves x directly to the next subset which doesn't intersect the input y.  y itself must be given as an indicator vector.
    def next_completely_avoiding(self,y):  # we find the next subset whose intersection with y is empty
        if not self.second:
            self.second = True
            if sum(y)==0:         #Uses the fact that if self.second is False, then x is all 1s.
                return True

        # check if x completely avoids y -- if it does not, move to the appropriate next one.
        for i in reversed(range(self.n)):
            if y[i]==1 and self.x[i]==1:  # x does not avoid y     #i is the rightmost 1 in both x and y
                self.x[i]=0
                i-=1
                while i>=0:
                    if y[i] == 0:
                        self.x[i] = 1
                    else:
                        self.x[i] = 0
                    i-=1
                return True

        #This x completely avoids y -- use this to move to the next appropriate one.
        else:
            for i in range(self.n):
                if self.x[i]==1:  # i is the leftmost 1.  (Note:  In this case, y[i] != 1.)
                    self.x[i]=0
                    i-=1
                    while i>=0:
                        #All the lesser bits which are not in y, change these to 1s.  (The rest will stay 0s.)
                        if y[i]==0:
                            self.x[i]=1
                        i-=1
                    return True
            else:  # there are no 1s left, so we are done
                return False


    #The next_next_not_containing method moves x directly to the next subset which doesn't contain the input z.  z itself must be given as an indicator vector.
    def next_not_containing(self,z):  # we find the next subset which does not contain z
        if not self.second:
            self.second = True    #Uses the fact that if self.second is False, then x is all 1s.


        #Preprocessing, Part 1:  Find i, the leftmost 1 of x.
        for a in range(self.n):
            if self.x[a]==1:
                i = a
                break
        else:
            return False     #If x is all 0s, then we're finished here.


        #Case 1:  If z[i]==1, then we know that the very next thing will work.
        if z[i]==1:
            return self.next()


        #Preprocessing, Part 2:  Find j, the leftmost 1 of z.
        for a in range(self.n):
            if z[a]==1:
                j = a
                break
        else:
            return False     #If z is all 0s, then we can't avoid containing it.  We're finished.


        #Case 2:  Search for the leftmost occurrence of x[k]==0 AND z[k]==1 conditioned on k>i.  If it exists, just take the very next thing.
        for a in range(i+1,self.n):
            if self.x[a]==0 and z[a]==1:
                return self.next()


        #Case 3:  There is no such k -- x[max(i,j)]<-0, and overwrite everything to the left to be 1, except for at j.
        #Example:  z = 00110010.  If x = 10110010, then x <- 11010010.  If x = 00001110, then x <- 11010110
        else:
            a = max([i,j])
            self.x[a] = 0
            a-=1
            while a>=0:
                self.x[a] = 1
                a-=1
            self.x[j] = 0         #(Redundant if j==a.)
            return True



    #The next_both_avoiding method moves x directly to the next subset which doesn't intersect the input y and does not contain the input z.  y and z themselves must be given as an indicator vector.
    def next_both_avoiding(self,y,z):
        if not self.second:
            self.second = True    #Uses the fact that if self.second is False, then x is all 1s.

        #Check if z intersects y -- if it does, this is simply a matter of avoiding y.
        for i in range(self.n):
            if y[i]==1 and z[i==1]:
                return self.next_completely_avoiding(y)
        else:
            for i in reversed(range(self.n)):
                if y[i]==1:
                    y_R = i         #y_R is the rightmost 1 of y
                    break
            else:
                return self.next_not_containing(z)        #If y is all 0s, we can't intersect it.

            for i in range(self.n):
                if z[i]==1:
                    z_L = i         #z_L is the leftmost 1 of z
                    break
            else:
                return False        #If z is all 0s, we can't avoid containing it.

        #We've preprocessed, and the case we are in is that y and z don't intersect, and y and z are nonempty.  Either z_L < y_R or z_L > y_R.

        #Case 1:  z_L < y_R.  Go to the next subset avoiding y, and remove z_L (but only if necessary!).
        if z_L < y_R:
            if self.next_completely_avoiding(y):
                for i in range(self.n):
                    if z[i]==1 and self.x[i]==0:     #If x doesn't contain z, we can stop right here.
                        return True
                else:
                    self.x[z_L] = 0                  #Otherwise, we need to pluck z_L.
                    return True
            else:
                return False

        #Case 2:  z_L > y_R.  Go to the next subset not containing z, and from there go to the next subset completely avoiding y.
        else:
            if self.next_not_containing(z):
                if self.next_completely_avoiding(y):
                    return True
                else:
                    return False
            else:
                return False



















#The generator class for connected subgraphs which contain a given root.  Takes a graph G as input along with a list "Eligible" of vertices which may be considered for inclusion in connected subgraphs and a single "root" vertex which must be contained by any such connected subgraph.  It proceeds through subgraphs by building up "levels," deciding to include or not include neighbors of the current level via the subset generator class.
#Like the Subset classes, RootedConnectedSubgraphs has multiple next methods.  One (next) for simply moving forward, one (next_within) for moving to the next subgraph which is contained in a given set, and one (next_within_not_containing) for moving to the next subgraph which is contained in a given set but does not contain another given set.
class RootedConnectedSubgraphs:

    def __init__(self,G,Eligible,root):
        self.G = G
        self.Eligible = Eligible
        self.root = root
        #Store the inputs.

        self.stack = []
        #stack is the stack of levels.  Each element is a list with information about both the current level and the next level.
        #The [0] item is "N_union", the set of vertices which have been considered thus far (either as "yes" or "no" vertices).
        #The [1] item is "S_union", the set of vertices which are currently included (as "yes" vertices).
        #The [2] item is "next_N", the set of vertices which we will consider on the next level (either as "yes" or "no" vertices).
        #The [3] item is "next_subset", the generator of subsets of next_N.  This determines the new S for the next level.

        self.subgraph = set([self.root])
        #subgraph is the set of "yes" vertices at any given time.  (We could call self.stack[-1][1] instead, but this is more explicit.)

        self.second = False
        #second is a boolean attribute used for keeping track of whether we have yet moved past the first subgraph (entire component).
        #second starts as False, is immediately switched to True upon a next method being called, and then stays True.



    def copy(self):
        clone = RootedConnectedSubgraphs(self.G,self.Eligible,self.root)
        if self.stack:
            clone.stack = [[copy.deepcopy(self.stack[0][0]), copy.deepcopy(self.stack[0][1]), copy.deepcopy(self.stack[0][2]), self.stack[0][3].copy()],]
            for i in range(1,len(self.stack)):
                clone.stack.append([copy.deepcopy(self.stack[i][0]), copy.deepcopy(self.stack[i][1]), copy.deepcopy(self.stack[i][2]), self.stack[i][3].copy()])
        else:
            clone.stack = []
        clone.subgraph = copy.deepcopy(self.subgraph)
        clone.second = self.second
        return clone



    #The first method must be called to build up the subgraph and the stack of levels.
    def first(self):
        N_union = set([self.root])            #The union of all vertices we have seen, marked either "yes" or "no."  (Formerly U.)
        S_union = set([self.root])            #The union of all vertices we have marked "yes".
        next_N = {x for x in nbrhood(self.G,self.Eligible,[self.root]) if x not in N_union}  #Vertices we can choose from for the next level.
        next_subset = Subset(list(next_N))       #The generator for the next "yes" vertices.
        next_subset.first()
        self.stack.append([N_union,S_union,next_N,next_subset])

        while True:  #Build up until we return.

            N_union_copy = copy.deepcopy(N_union)
            N_union = N_union | next_N

            #Get next S ("yes" vertices for next level).
            next_S = set(next_subset.x_converted())
            next_subset.second = True

            #Add this S to the chosen vertices.  If next_N is empty, then this is our next connected subgraph.
            S_union = S_union | next_S

            # we compute the N for the next level
            next_N={x for x in nbrhood(self.G,self.Eligible,next_S) if x not in N_union}

            if next_N:  # next_N is not empty
                next_subset = Subset(list(next_N))
                next_subset.first()
                self.stack.append([N_union,S_union,next_N,next_subset])  # the next level

            else:
                # we have completed this subgraph, so stop here for it to be called
                self.subgraph = S_union
                return True



    def next(self):
        if not self.second:
            self.second = True
            return True
        else:
            while True:
                if self.stack:
                    level=self.stack[-1]
                    N_union=level[0]
                    S_union = level[1]
                    next_N = level[2]
                    next_subset=level[3]

                    #Get next S ("yes" vertices for next level).
                    if not next_subset.next():
                        self.stack.pop()
                        continue
                    next_S = set(next_subset.x_converted())

                    #Get next N_union.  (Need to do this before computing next next_N.)
                    N_union_copy = copy.deepcopy(N_union)
                    N_union = N_union | next_N

                    #Add this S to the chosen vertices.  If next_N is empty, then this is our next connected subgraph.
                    S_union = S_union | next_S

                    # we compute the N for the next level
                    next_N={x for x in nbrhood(self.G,self.Eligible,next_S) if x not in N_union}

                    if next_N:  # next_N is not empty
                        next_subset = Subset(list(next_N))
                        next_subset.first()
                        self.stack.append([N_union,S_union,next_N,next_subset])  # the next level
                        continue

                    else:
                        # we have completed this subgraph, so stop here for it to be called
                        self.subgraph = S_union
                        return True

                else:
                    return False



    def next_within(self,E):     #E = "eligible" vertices = complement of "full" vertices.
        if not self.second:
            self.second = True
            if set(self.subgraph) <= set(E):
                return True

        this_E = set(E) & set(self.Eligible)

        I = set(self.G.vertices()) - this_E
        if self.stack and self.stack[-1][1] & I:
            for i in range(len(self.stack)):
                if self.stack[i][1] & I:            #i is the top level which doesn't avoid I -- clear all levels here and below.
                    for j in reversed(range(i,len(self.stack))):
                        del self.stack[j]
                    break

        #At this point, we've backtracked to the appropriate level (possibly not at all).
        #Now we move forward in the next level (using the generator store by the current top level) to the appropriate subset.
        while True:
            if self.stack:
                level=self.stack[-1]
                N_union=level[0]
                S_union = level[1]
                next_N = level[2]
                next_subset=level[3]

                #Get next S ("yes" vertices for next level).
                if not next_subset.next_completely_avoiding(convert_to_vector(next_subset.given_list,I)):
                    self.stack.pop()
                    continue
                next_S = set(next_subset.x_converted())

                #Get next N_union.  (Need to do this before computing next next_N.)
                N_union_copy = copy.deepcopy(N_union)     #For computing next_N.  We can use level[0] in that line instead of N_union_copy.
                N_union = N_union | next_N

                #Add this S to the chosen vertices.  If next_N is empty, then this is our next connected subgraph.
                S_union = S_union | next_S

                # we compute the N for the next level
                next_N={x for x in nbrhood(self.G,this_E,next_S) if x not in N_union}

                if next_N:  # next_N is not empty
                    next_subset = Subset(list(next_N))
                    next_subset.first()
                    self.stack.append([N_union,S_union,next_N,next_subset])  # the next level
                    continue

                else:
                    # we have completed this subgraph, so stop here for it to be called
                    self.subgraph = S_union
                    return True

            else:
                return False



    def next_within_not_containing(self,E,D):     #E = "eligible" vertices = complement of "full" vertices; D = "don't contain" vertices.
        if not self.second:
            self.second = True
            if (not set(D) <= set(self.subgraph)) and set(self.subgraph) <= set(E):
                return True

        this_E = set(E) & set(self.Eligible)

        I = set(self.G.vertices()) - this_E
        if self.stack and (self.stack[-1][1] & I or set(D) <= self.stack[-1][1]):
            for i in range(len(self.stack)):
                if self.stack[i][1] & I or set(D) <= self.stack[i][1]:      #i is the top level which doesn't avoid I or D properly.
                    for j in reversed(range(i,len(self.stack))):
                        del self.stack[j]
                    break

        #At this point, we've backtracked to the appropriate level (possibly not at all).
        #Now we move forward in the next level (using the generator stored by the current top level) to the appropriate subset.
        while True:
            if self.stack:
                level=self.stack[-1]
                N_union=level[0]
                S_union = level[1]
                next_N = level[2]
                next_subset=level[3]

                #Get next S ("yes" vertices for next level).
                y = convert_to_vector(next_subset.given_list,I)
                z_test = (N_union - S_union) & set(D) or set(D) - (N_union | next_N)
                if z_test:         #If we have already said "no" to a vertex in D, or if we can't possibly contain D on the next level.
                    if not next_subset.next_completely_avoiding(y):
                        self.stack.pop()
                        continue
                else:
                    z = convert_to_vector(next_subset.given_list,set(D)-S_union)
                    if not next_subset.next_both_avoiding(y,z):
                        self.stack.pop()
                        continue
                next_S = set(next_subset.x_converted())

                #Get next N_union.  (Need to do this before computing next next_N.)
                N_union_copy = copy.deepcopy(N_union)     #For computing next_N.  We can use level[0] in that line instead of N_union_copy.
                N_union = N_union | next_N

                #Add this S to the chosen vertices.  If next_N is empty, then this is our next connected subgraph.
                S_union = S_union | next_S

                # we compute the N for the next level
                next_N={x for x in nbrhood(self.G,this_E,next_S) if x not in N_union}

                if next_N:  # next_N is not empty
                    next_subset = Subset(list(next_N))
                    next_subset.first()
                    self.stack.append([N_union,S_union,next_N,next_subset])  # the next level
                    continue

                else:
                    # we have completed this subgraph, so stop here for it to be called
                    self.subgraph = S_union
                    return True

            else:
                return False















#A copy of ConnectedSubgraphs, but a list of lists of rooted subgraphs is provided as the cascade.  restrictionList is a list of lists of subsets; restrictionList[vtx] is the list of subgraphs rooted at vtx which are not CNS-reducible.
class RestrictedConnectedSubgraphs_L:

    #ordering is a list of initially eligible vertices, the order of which determines the columns of the generator cascade.
    def __init__(self,G,ordering,restrictionList):
        self.G = G
        self.ordering = ordering
        self.generator_cascade = restrictionList[:]
        
        #containment_dict is a dictionary holding a trimmed version of the cascade for each vertex.
        self.containment_dict = {}
        for v in self.ordering:
            self.containment_dict[v] = [[x for x in column if v in x] for column in self.generator_cascade]

        #finished_roots is the list of currently exhausted columns.
        self.finished_roots = [x for x in range(len(ordering)) if not restrictionList[x]]

        #subgraph is the current connected subgraph (set of vertices).
        self.subgraph = None



    def copy(self):
        leftover_restrictionList = [] #initialize
        for i in range(len(self.generator_cascade)): #overwrite
            leftover_restrictionList.append(self.generator_cascade[i][:])
        clone = RestrictedConnectedSubgraphs_L(self.G,self.ordering,leftover_restrictionList)
        clone.finished_roots = copy.deepcopy(self.finished_roots)
        clone.subgraph = copy.deepcopy(self.subgraph)
        ##clone.containment_dict = copy.deepcopy(self.containment_dict)
        #for v in clone.ordering:
            #clone.containment_dict[v] = [[x for x in column if v in x] for column in clone.generator_cascade]
        return clone



    def update_restriction_list(self,f,CA,maxIS=True):
        self.generator_cascade,num_pruned = restriction_list_updater_inprocess(self.G,self.ordering,f,CA,self.generator_cascade,singletons=True,exclusions=[],maxIS=maxIS)
        for v in self.ordering:
            self.containment_dict[v] = [[x for x in column if v in x] for column in self.generator_cascade]
        self.finished_roots = [x for x in range(len(self.ordering)) if not self.generator_cascade[x]]
        return num_pruned




    #def intersect_restriction_list(self,restrictionList):
        #for i in range(len(self.generator_cascade)):
            #for j in reversed(range(len(self.generator_cascade[i]))):
                #if self.generator_cascade[i][j] in restrictionList[i]:
                    #continue
                #else:
                    #del self.generator_cascade[i][j]



    #As of 2018.06.04, .next and .next_with_root are not used in this class.



    ##L is the list of list sizes
    ##E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]])
    ##The plain next method has a built-in root-chooser and does not skip any subsets.  (L and E are for root selection.)
    #def next(self,L,E):
        #while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.
            ##Check to see whether all the generators in the cascade have finished.
            #if set(E) <= set(self.finished_roots):
                #return False

            ##Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            ##H is ordering with ineligible or finished roots removed.
            #H = [x for x in self.ordering if x in E and x not in self.finished_roots]
##This if statement is redundant.
            #if H:
                #m = min([L[y] for y in H]) #min needs to be over just H.
                #for x in H:
                    #if L[x]==m:
                        #root=x
                        #break
            #else:
                #return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.

            ##Call the next subgraph from the selected root.
            #column = self.generator_cascade[root]
            #if column:
                #self.subgraph = set(column[0])
                #del column[0]
                ##.next() is never called
                ##u = min(self.subgraph)
                ##for v in self.subgraph:
                    ##self.containment_dict[v][u].remove(self.subgraph)
                #return True
            #else:        #If this column is exhausted, make note and try again.
                #self.finished_roots.append(root)
                #continue





    ##This method takes a pre-chosen root, but skips no subgraphs.
    #def next_with_root(self,root):
        #while True:
            #if root in self.finished_roots:
                #return False

            ##Call the next subgraph from the selected root.
            #column = self.generator_cascade[root]
            #if column:
                #self.subgraph = set(column[0])
                #del column[0]
                ##.next_with_root() is never called
                ##u = min(self.subgraph)
                ##for v in self.subgraph:
                    ##self.containment_dict[v][u].remove(self.subgraph)
                #return True
            #else:
                #self.finished_roots.append(root)
                #continue



    #L is the list of list sizes
    #E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]]).  E is always nonempty when this is called.
    #This method has a built-in root-chooser and skips subgraphs which are not contained in E.
    def next_within(self,L,E):

        while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.

            #Check whether it is in the realm of possibility to fill up the remaining vertices.  There is a vertex in E such that it and every earlier-ordered vertex in E are finished as roots if and only if the earliest ordered vertex in E is a finished root.
            if [x for x in self.ordering if x in E][0] in self.finished_roots:
                return False


            #Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            #H is ordering with ineligible or finished roots removed.
            H = [x for x in self.ordering if x in E and x not in self.finished_roots]
            if H:
                m = min([L[y] for y in H]) #min needs to be over just H.
                for x in H:
                    if L[x]==m:
                        root=x
                        break
            else:
                return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.

            #Call the next subgraph from the selected root.
            column = self.generator_cascade[root]
            while column:
                if set(column[0]) <= set(E):
                    self.subgraph = set(column[0])
                    u = min(column[0])
                    for v in self.subgraph:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
                    return True
                else:
                    u = min(column[0])
                    for v in column[0]:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
            #If this column is exhausted, make note and try again.
            self.finished_roots.append(root)
            continue



    #L is the list of list sizes
    #E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]])
    #D is the list of vertices receiving the top color in the coloring certificate -- we can skip colorability classes which contain D.
    #This method has a built-in root-chooser and skips subgraphs which are not contained in E or do contain D.
    def next_within_not_containing(self,L,E,D):

        while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.

            #Check whether it is in the realm of possibility to fill up the remaining vertices.  There is a vertex in E such that it and every earlier-ordered vertex in E are finished as roots if and only if the earliest ordered vertex in E is a finished root.
            if [x for x in self.ordering if x in E][0] in self.finished_roots:
                return False

            #Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            #H is ordering with ineligible or finished roots removed.
            H = [x for x in self.ordering if x in E and x not in self.finished_roots]
            if H:
                m = min([L[y] for y in H]) #min needs to be over just H.
                for x in H:
                    if L[x]==m:
                        root=x
                        break
            else:
                return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.

            #Call the next subgraph from the selected root.
            column = self.generator_cascade[root]
            while column:
                if set(column[0]) <= set(E) and not set(D) <= set(column[0]):
                    self.subgraph = set(column[0])
                    u = min(column[0])
                    for v in self.subgraph:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
                    return True
                else:
                    u = min(column[0])
                    for v in column[0]:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
            #If this column is exhausted, make note and try again.
            self.finished_roots.append(root)
            continue















#A copy of ConnectedSubgraphs, but a list of lists of rooted subgraphs is provided as the cascade.  restrictionList is a list of lists of subsets; restrictionList[vtx] is the list of subgraphs rooted at vtx which are not CNS-reducible.
class RestrictedConnectedSubgraphs_LE:

    #ordering is a list of initially eligible vertices, the order of which determines the columns of the generator cascade.
    def __init__(self,G,ordering,restrictionList):
        self.G = G
        self.ordering = ordering
        self.generator_cascade = restrictionList[:]
        
        #containment_dict is a dictionary holding a trimmed version of the cascade for each vertex.
        self.containment_dict = {}
        for v in self.ordering:
            self.containment_dict[v] = [[x for x in column if v in x] for column in self.generator_cascade]

        #finished_roots is the list of currently exhausted columns.
        self.finished_roots = [x for x in range(len(ordering)) if not restrictionList[x]]

        #subgraph is the current connected subgraph (set of vertices).
        self.subgraph = None



    def copy(self):
        leftover_restrictionList = [] #initialize
        for i in range(len(self.generator_cascade)): #overwrite
            leftover_restrictionList.append(self.generator_cascade[i][:])
        clone = RestrictedConnectedSubgraphs_LE(self.G,self.ordering,leftover_restrictionList)
        clone.finished_roots = copy.deepcopy(self.finished_roots)
        clone.subgraph = copy.deepcopy(self.subgraph)
        ##clone.containment_dict = copy.deepcopy(self.containment_dict)
        #for v in clone.ordering:
            #clone.containment_dict[v] = [[x for x in column if v in x] for column in clone.generator_cascade]
        return clone



    def update_restriction_list(self,f,CA,maxIS=True):
        self.generator_cascade,num_pruned = restriction_list_updater_inprocess(self.G,self.ordering,f,CA,self.generator_cascade,singletons=True,exclusions=[],maxIS=maxIS)
        for v in self.ordering:
            self.containment_dict[v] = [[x for x in column if v in x] for column in self.generator_cascade]
        self.finished_roots = [x for x in range(len(self.ordering)) if not self.generator_cascade[x]]
        return num_pruned




    #def intersect_restriction_list(self,restrictionList):
        #for i in range(len(self.generator_cascade)):
            #for j in reversed(range(len(self.generator_cascade[i]))):
                #if self.generator_cascade[i][j] in restrictionList[i]:
                    #continue
                #else:
                    #del self.generator_cascade[i][j]



    #As of 2018.06.04, .next and .next_with_root are not used in this class.



    ##L is the list of list sizes
    ##E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]])
    ##The plain next method has a built-in root-chooser and does not skip any subsets.  (L and E are for root selection.)
    #def next(self,L,E):
        #while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.
            ##Check to see whether all the generators in the cascade have finished.
            #if set(E) <= set(self.finished_roots):
                #return False

            ##Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            ##H is ordering with ineligible or finished roots removed.
            #H = [x for x in self.ordering if x in E and x not in self.finished_roots]
##This if statement is redundant.
            #if H:
                #m = min([L[y] for y in H]) #min needs to be over just H.
                #for x in H:
                    #if L[x]==m:
                        #root=x
                        #break
            #else:
                #return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.

            ##Call the next subgraph from the selected root.
            #column = self.generator_cascade[root]
            #if column:
                #self.subgraph = set(column[0])
                #del column[0]
                ##.next() is never called
                ##u = min(self.subgraph)
                ##for v in self.subgraph:
                    ##self.containment_dict[v][u].remove(self.subgraph)
                #return True
            #else:        #If this column is exhausted, make note and try again.
                #self.finished_roots.append(root)
                #continue





    ##This method takes a pre-chosen root, but skips no subgraphs.
    #def next_with_root(self,root):
        #while True:
            #if root in self.finished_roots:
                #return False

            ##Call the next subgraph from the selected root.
            #column = self.generator_cascade[root]
            #if column:
                #self.subgraph = set(column[0])
                #del column[0]
                ##.next_with_root() is never called
                ##u = min(self.subgraph)
                ##for v in self.subgraph:
                    ##self.containment_dict[v][u].remove(self.subgraph)
                #return True
            #else:
                #self.finished_roots.append(root)
                #continue



    #L is the list of list sizes
    #E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]]).  E is always nonempty when this is called.
    #This method has a built-in root-chooser and skips subgraphs which are not contained in E.
    def next_within(self,L,E):

        while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.

            #Check whether it is in the realm of possibility to fill up the remaining vertices.  There is a vertex in E such that it and every earlier-ordered vertex in E are finished as roots if and only if the earliest ordered vertex in E is a finished root.
            if [x for x in self.ordering if x in E][0] in self.finished_roots:
                return False


            #Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            #H is ordering with ineligible or finished roots removed.
            H = [x for x in self.ordering if x in E and x not in self.finished_roots]
            if H:
                m = min([L[y] for y in H]) #min needs to be over just H.
                for x in H:
                    if L[x]==m:
                        root=x
                        break
            else:
                return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.
            
            flag = True
            for i in range(len(self.containment_dict[root])):
                if self.containment_dict[root][i]:
                    root = i
                    flag = False
                    break
                
            if flag:
                return False   #In this case, there are no remaining subgraphs which contain root (which needs more colors).

            #Call the next subgraph from the selected root.
            column = self.generator_cascade[root]
            while column:
                if set(column[0]) <= set(E):
                    self.subgraph = set(column[0])
                    u = min(column[0])
                    for v in self.subgraph:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
                    return True
                else:
                    u = min(column[0])
                    for v in column[0]:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
            #If this column is exhausted, make note and try again.
            self.finished_roots.append(root)
            continue



    #L is the list of list sizes
    #E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]])
    #D is the list of vertices receiving the top color in the coloring certificate -- we can skip colorability classes which contain D.
    #This method has a built-in root-chooser and skips subgraphs which are not contained in E or do contain D.
    def next_within_not_containing(self,L,E,D):

        while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.

            #Check whether it is in the realm of possibility to fill up the remaining vertices.  There is a vertex in E such that it and every earlier-ordered vertex in E are finished as roots if and only if the earliest ordered vertex in E is a finished root.
            if [x for x in self.ordering if x in E][0] in self.finished_roots:
                return False

            #Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            #H is ordering with ineligible or finished roots removed.
            H = [x for x in self.ordering if x in E and x not in self.finished_roots]
            if H:
                m = min([L[y] for y in H]) #min needs to be over just H.
                for x in H:
                    if L[x]==m:
                        root=x
                        break
            else:
                return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.
            
            flag = True
            for i in range(len(self.containment_dict[root])):
                if self.containment_dict[root][i]:
                    root = i
                    flag = False
                    break
                
            if flag:
                return False   #In this case, there are no remaining subgraphs which contain root (which needs more colors).

            #Call the next subgraph from the selected root.
            column = self.generator_cascade[root]
            while column:
                if set(column[0]) <= set(E) and not set(D) <= set(column[0]):
                    self.subgraph = set(column[0])
                    u = min(column[0])
                    for v in self.subgraph:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
                    return True
                else:
                    u = min(column[0])
                    for v in column[0]:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
            #If this column is exhausted, make note and try again.
            self.finished_roots.append(root)
            continue















#A copy of ConnectedSubgraphs, but a list of lists of rooted subgraphs is provided as the cascade.  restrictionList is a list of lists of subsets; restrictionList[vtx] is the list of subgraphs rooted at vtx which are not CNS-reducible.
class RestrictedConnectedSubgraphs_LD:

    #ordering is a list of initially eligible vertices, the order of which determines the columns of the generator cascade.
    def __init__(self,G,ordering,restrictionList):
        self.G = G
        self.ordering = ordering
        self.generator_cascade = restrictionList[:]
        
        #containment_dict is a dictionary holding a trimmed version of the cascade for each vertex.
        self.containment_dict = {}
        for v in self.ordering:
            self.containment_dict[v] = [[x for x in column if v in x] for column in self.generator_cascade]

        #finished_roots is the list of currently exhausted columns.
        self.finished_roots = [x for x in range(len(ordering)) if not restrictionList[x]]

        #subgraph is the current connected subgraph (set of vertices).
        self.subgraph = None



    def copy(self):
        leftover_restrictionList = [] #initialize
        for i in range(len(self.generator_cascade)): #overwrite
            leftover_restrictionList.append(self.generator_cascade[i][:])
        clone = RestrictedConnectedSubgraphs_LD(self.G,self.ordering,leftover_restrictionList)
        clone.finished_roots = copy.deepcopy(self.finished_roots)
        clone.subgraph = copy.deepcopy(self.subgraph)
        ##clone.containment_dict = copy.deepcopy(self.containment_dict)
        #for v in clone.ordering:
            #clone.containment_dict[v] = [[x for x in column if v in x] for column in clone.generator_cascade]
        return clone



    def update_restriction_list(self,f,CA,maxIS=True):
        self.generator_cascade,num_pruned = restriction_list_updater_inprocess(self.G,self.ordering,f,CA,self.generator_cascade,singletons=True,exclusions=[],maxIS=maxIS)
        for v in self.ordering:
            self.containment_dict[v] = [[x for x in column if v in x] for column in self.generator_cascade]
        self.finished_roots = [x for x in range(len(self.ordering)) if not self.generator_cascade[x]]
        return num_pruned




    #def intersect_restriction_list(self,restrictionList):
        #for i in range(len(self.generator_cascade)):
            #for j in reversed(range(len(self.generator_cascade[i]))):
                #if self.generator_cascade[i][j] in restrictionList[i]:
                    #continue
                #else:
                    #del self.generator_cascade[i][j]



    #As of 2018.06.04, .next and .next_with_root are not used in this class.



    ##L is the list of list sizes
    ##E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]])
    ##The plain next method has a built-in root-chooser and does not skip any subsets.  (L and E are for root selection.)
    #def next(self,L,E):
        #while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.
            ##Check to see whether all the generators in the cascade have finished.
            #if set(E) <= set(self.finished_roots):
                #return False

            ##Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            ##H is ordering with ineligible or finished roots removed.
            #H = [x for x in self.ordering if x in E and x not in self.finished_roots]
##This if statement is redundant.
            #if H:
                #m = min([L[y] for y in H]) #min needs to be over just H.
                #for x in H:
                    #if L[x]==m:
                        #root=x
                        #break
            #else:
                #return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.

            ##Call the next subgraph from the selected root.
            #column = self.generator_cascade[root]
            #if column:
                #self.subgraph = set(column[0])
                #del column[0]
                ##.next() is never called
                ##u = min(self.subgraph)
                ##for v in self.subgraph:
                    ##self.containment_dict[v][u].remove(self.subgraph)
                #return True
            #else:        #If this column is exhausted, make note and try again.
                #self.finished_roots.append(root)
                #continue





    ##This method takes a pre-chosen root, but skips no subgraphs.
    #def next_with_root(self,root):
        #while True:
            #if root in self.finished_roots:
                #return False

            ##Call the next subgraph from the selected root.
            #column = self.generator_cascade[root]
            #if column:
                #self.subgraph = set(column[0])
                #del column[0]
                ##.next_with_root() is never called
                ##u = min(self.subgraph)
                ##for v in self.subgraph:
                    ##self.containment_dict[v][u].remove(self.subgraph)
                #return True
            #else:
                #self.finished_roots.append(root)
                #continue



    #L is the list of list sizes
    #E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]]).  E is always nonempty when this is called.
    #This method has a built-in root-chooser and skips subgraphs which are not contained in E.
    def next_within(self,L,E):

        while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.

            #Check whether it is in the realm of possibility to fill up the remaining vertices.  There is a vertex in E such that it and every earlier-ordered vertex in E are finished as roots if and only if the earliest ordered vertex in E is a finished root.
            if [x for x in self.ordering if x in E][0] in self.finished_roots:
                return False


            #Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            #H is ordering with ineligible or finished roots removed.
            H = [x for x in self.ordering if x in E]# and x not in self.finished_roots]
            if H:
                m = min([L[y] for y in H]) #min needs to be over just H.
                for x in H:
                    if L[x]==m:
                        root=x
                        break
            else:
                return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.
            
            for column in self.containment_dict[root]:
                while column:
                    if set(column[0]) <= set(E):
                        self.subgraph = set(column[0])
                        u = min(column[0])
                        for v in self.subgraph:
                            self.containment_dict[v][u].remove(self.subgraph)
                        #del column[0]#Don't do this!  column[0] was already deleted via the above line at some point!
                        self.generator_cascade[u].remove(self.subgraph)
                        #print "\nH:",H
                        #print "L:",[L[x] for x in H]
                        #print "root:",root
                        #print "subgraph:",self.subgraph
                        #if not root in self.subgraph:
                            #print "! root is",root,"but not contained in subgraph",self.subgraph
                        return True
                    else:
                        subgraph = set(column[0])
                        u = min(column[0])
                        for v in subgraph:
                            self.containment_dict[v][u].remove(subgraph)
                        self.generator_cascade[u].remove(subgraph)
                        #del column[0]
            #If there are no subgraphs containing root that fit inside E, then there are no full super list assignments.
            return False



    #L is the list of list sizes
    #E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]])
    #D is the list of vertices receiving the top color in the coloring certificate -- we can skip colorability classes which contain D.
    #This method has a built-in root-chooser and skips subgraphs which are not contained in E or do contain D.
    def next_within_not_containing(self,L,E,D):

        while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.

            #Check whether it is in the realm of possibility to fill up the remaining vertices.  There is a vertex in E such that it and every earlier-ordered vertex in E are finished as roots if and only if the earliest ordered vertex in E is a finished root.
            if [x for x in self.ordering if x in E][0] in self.finished_roots:
                return False

            #Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            #H is ordering with ineligible or finished roots removed.
            H = [x for x in self.ordering if x in E]# and x not in self.finished_roots]
            if H:
                m = min([L[y] for y in H]) #min needs to be over just H.
                for x in H:
                    if L[x]==m:
                        root=x
                        break
            else:
                return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.
            
            for column in self.containment_dict[root]:
                while column:
                    if set(column[0]) <= set(E) and not set(D) <= set(column[0]):
                        self.subgraph = set(column[0])
                        u = min(column[0])
                        for v in self.subgraph:
                            self.containment_dict[v][u].remove(self.subgraph)
                        #del column[0]#Don't do this!  column[0] was already deleted via the above line at some point!
                        self.generator_cascade[u].remove(self.subgraph)
                        return True
                    else:
                        subgraph = set(column[0])
                        u = min(column[0])
                        for v in subgraph:
                            self.containment_dict[v][u].remove(subgraph)
                        self.generator_cascade[u].remove(subgraph)
                        #del column[0]
            #If there are no subgraphs containing root that fit inside E, then there are no full super list assignments.
            return False














#A copy of ConnectedSubgraphs, but a list of lists of rooted subgraphs is provided as the cascade.  restrictionList is a list of lists of subsets; restrictionList[vtx] is the list of subgraphs rooted at vtx which are not CNS-reducible.
class RestrictedConnectedSubgraphs_E:

    #ordering is a list of initially eligible vertices, the order of which determines the columns of the generator cascade.
    def __init__(self,G,ordering,restrictionList):
        self.G = G
        self.ordering = ordering
        self.generator_cascade = restrictionList[:]
        
        #containment_dict is a dictionary holding a trimmed version of the cascade for each vertex.
        self.containment_dict = {}
        for v in self.ordering:
            self.containment_dict[v] = [[x for x in column if v in x] for column in self.generator_cascade]

        #finished_roots is the list of currently exhausted columns.
        self.finished_roots = [x for x in range(len(ordering)) if not restrictionList[x]]

        #subgraph is the current connected subgraph (set of vertices).
        self.subgraph = None



    def copy(self):
        leftover_restrictionList = [] #initialize
        for i in range(len(self.generator_cascade)): #overwrite
            leftover_restrictionList.append(self.generator_cascade[i][:])
        clone = RestrictedConnectedSubgraphs_E(self.G,self.ordering,leftover_restrictionList)
        clone.finished_roots = copy.deepcopy(self.finished_roots)
        clone.subgraph = copy.deepcopy(self.subgraph)
        ##clone.containment_dict = copy.deepcopy(self.containment_dict)
        #for v in clone.ordering:
            #clone.containment_dict[v] = [[x for x in column if v in x] for column in clone.generator_cascade]
        return clone



    def update_restriction_list(self,f,CA,maxIS=True):
        self.generator_cascade,num_pruned = restriction_list_updater_inprocess(self.G,self.ordering,f,CA,self.generator_cascade,singletons=True,exclusions=[],maxIS=maxIS)
        for v in self.ordering:
            self.containment_dict[v] = [[x for x in column if v in x] for column in self.generator_cascade]
        self.finished_roots = [x for x in range(len(self.ordering)) if not self.generator_cascade[x]]
        return num_pruned




    #def intersect_restriction_list(self,restrictionList):
        #for i in range(len(self.generator_cascade)):
            #for j in reversed(range(len(self.generator_cascade[i]))):
                #if self.generator_cascade[i][j] in restrictionList[i]:
                    #continue
                #else:
                    #del self.generator_cascade[i][j]



    #As of 2018.06.04, .next and .next_with_root are not used in this class.



    ##L is the list of list sizes
    ##E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]])
    ##The plain next method has a built-in root-chooser and does not skip any subsets.  (L and E are for root selection.)
    #def next(self,L,E):
        #while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.
            ##Check to see whether all the generators in the cascade have finished.
            #if set(E) <= set(self.finished_roots):
                #return False

            ##Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            ##H is ordering with ineligible or finished roots removed.
            #H = [x for x in self.ordering if x in E and x not in self.finished_roots]
##This if statement is redundant.
            #if H:
                #m = min([L[y] for y in H]) #min needs to be over just H.
                #for x in H:
                    #if L[x]==m:
                        #root=x
                        #break
            #else:
                #return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.

            ##Call the next subgraph from the selected root.
            #column = self.generator_cascade[root]
            #if column:
                #self.subgraph = set(column[0])
                #del column[0]
                ##.next() is never called
                ##u = min(self.subgraph)
                ##for v in self.subgraph:
                    ##self.containment_dict[v][u].remove(self.subgraph)
                #return True
            #else:        #If this column is exhausted, make note and try again.
                #self.finished_roots.append(root)
                #continue





    ##This method takes a pre-chosen root, but skips no subgraphs.
    #def next_with_root(self,root):
        #while True:
            #if root in self.finished_roots:
                #return False

            ##Call the next subgraph from the selected root.
            #column = self.generator_cascade[root]
            #if column:
                #self.subgraph = set(column[0])
                #del column[0]
                ##.next_with_root() is never called
                ##u = min(self.subgraph)
                ##for v in self.subgraph:
                    ##self.containment_dict[v][u].remove(self.subgraph)
                #return True
            #else:
                #self.finished_roots.append(root)
                #continue



    #L is the list of list sizes
    #E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]]).  E is always nonempty when this is called.
    #This method has a built-in root-chooser and skips subgraphs which are not contained in E.
    def next_within(self,L,E):

        while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.

            #Check whether it is in the realm of possibility to fill up the remaining vertices.  There is a vertex in E such that it and every earlier-ordered vertex in E are finished as roots if and only if the earliest ordered vertex in E is a finished root.
            if [x for x in self.ordering if x in E][0] in self.finished_roots:
                return False


            #Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            #H is ordering with ineligible or finished roots removed.
            H = [x for x in self.ordering if x in E and x not in self.finished_roots]
            if H:
                root = H[0]
            else:
                return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.

            #Call the next subgraph from the selected root.
            column = self.generator_cascade[root]
            while column:
                if set(column[0]) <= set(E):
                    self.subgraph = set(column[0])
                    u = min(column[0])
                    for v in self.subgraph:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
                    return True
                else:
                    u = min(column[0])
                    for v in column[0]:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
            #If this column is exhausted, make note and try again.
            self.finished_roots.append(root)
            continue



    #L is the list of list sizes
    #E is the list of eligible vertices ([x for x in G.vertices() if L[x]<f[x]])
    #D is the list of vertices receiving the top color in the coloring certificate -- we can skip colorability classes which contain D.
    #This method has a built-in root-chooser and skips subgraphs which are not contained in E or do contain D.
    def next_within_not_containing(self,L,E,D):

        while True: #If the root we choose turns out to have an exhausted generator, then we'll need to try again from the top of the loop.

            #Check whether it is in the realm of possibility to fill up the remaining vertices.  There is a vertex in E such that it and every earlier-ordered vertex in E are finished as roots if and only if the earliest ordered vertex in E is a finished root.
            if [x for x in self.ordering if x in E][0] in self.finished_roots:
                return False

            #Choose the next root.  Currently:  Least L, Least Index selection (according to vertex_ordering)
            #H is ordering with ineligible or finished roots removed.
            H = [x for x in self.ordering if x in E and x not in self.finished_roots]
            if H:
                root = H[0]
            else:
                return False    #In this case, E is a subset of finished_roots.  There are no available roots for subgraphs not contained in E.

            #Call the next subgraph from the selected root.
            column = self.generator_cascade[root]
            while column:
                if set(column[0]) <= set(E) and not set(D) <= set(column[0]):
                    self.subgraph = set(column[0])
                    u = min(column[0])
                    for v in self.subgraph:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
                    return True
                else:
                    u = min(column[0])
                    for v in column[0]:
                        self.containment_dict[v][u].remove(column[0])
                    del column[0]
            #If this column is exhausted, make note and try again.
            self.finished_roots.append(root)
            continue




















#G_copy is the input graph (nonempty) on vertices range(n) and tvec is the list of powers for the corresponding vertex variables.
#WARNING:  .coefficient(x**0) will return the entire polynomial (including all the monomials containing x).  Use .coefficient({x:0}) instead.
def coeff_extract(G,tvec):
    stringlist = ['x'+str(v) for v in G.vertices()]
    P=PolynomialRing(ZZ, G.order(), stringlist)
    try:
        varlist = P.objgens()[1]
        partialpoly = varlist[0]**0
    except IndexError:
        G_copy.show()
        print tvec
        print "Uh-oh!"
    for v in G.vertices():
        prod = varlist[0]**0
        for factor in [varlist[v]-varlist[u] for u in G.neighbors(v) if v < u]:
            prod *= factor
        nextpartialpoly = 0
        for x in range(0,tvec[v]+1):
            nextpartialpoly += partialpoly.coefficient({varlist[v]:x})*prod.coefficient({varlist[v]:(tvec[v]-x)})
        partialpoly = nextpartialpoly
    return partialpoly











def coeff_extract_GF7(G,tvec):
    stringlist = ['x'+str(v) for v in G.vertices()]
    P=PolynomialRing(GF(7), G.order(), stringlist)
    try:
        varlist = P.objgens()[1]
        partialpoly = varlist[0]**0
    except IndexError:
        G_copy.show()
        print tvec
        print "Uh-oh!"
    for v in G.vertices():
        prod = varlist[0]**0
        for factor in [varlist[v]-varlist[u] for u in G.neighbors(v) if v < u]:
            prod *= factor
        nextpartialpoly = 0
        for x in range(0,tvec[v]+1):
            nextpartialpoly += partialpoly.coefficient({varlist[v]:x})*prod.coefficient({varlist[v]:(tvec[v]-x)})
        partialpoly = nextpartialpoly
    return partialpoly











#G_copy is the input graph (nonempty) on vertices range(n) and tvec is the list of powers for the corresponding vertex variables.
def coeff_extract_old(G_copy,tvec):
    G = Graph(G_copy)
    stringlist = ['x'+str(v) for v in G.vertices()]
    P=PolynomialRing(GF(7), G.order(), stringlist)
    try:
        varlist = P.objgens()[1]
        partialpoly = varlist[0]**0
    except IndexError:
        G_copy.show()
        print tvec
        print "Uh-oh!"
    for v in G.vertices():
        prod = 1
        for factor in [varlist[v]-varlist[u] for u in G.neighbors(v)]:
            prod *= factor
        G.delete_vertex(v)
        nextpartialpoly = partialpoly*prod
        partialpoly = nextpartialpoly.coefficient(varlist[v]**tvec[v])
    return partialpoly











#G_copy is the input graph (nonempty) on vertices range(n) and tvec is the list of powers for the corresponding vertex variables.
def coeff_extract_oldZZ(G_copy,tvec):
    G = Graph(G_copy)
    stringlist = ['x'+str(v) for v in G.vertices()]
    P=PolynomialRing(ZZ, G.order(), stringlist)
    try:
        varlist = P.objgens()[1]
        partialpoly = varlist[0]**0
    except IndexError:
        G_copy.show()
        print tvec
        print "Uh-oh!"
    for v in G.vertices():
        prod = 1
        for factor in [varlist[v]-varlist[u] for u in G.neighbors(v)]:
            prod *= factor
        G.delete_vertex(v)
        nextpartialpoly = partialpoly*prod
        partialpoly = nextpartialpoly.coefficient(varlist[v]**tvec[v])
    return partialpoly





def CNS_blind(G_copy, f):
    disparity = sum(f)-len(f) - G_copy.size()
    if disparity<0:
        #print "Can't apply CNS to G,f."
        return False
    else:
        red = False
        for fillvec in IntegerListsLex(disparity,length=len(f),ceiling=[f[x]-1 for x in range(len(f))]):
            tvec = [f[x]-1-fillvec[x] for x in range(len(f))]
            if coeff_extract(G_copy,tvec):
                red = True
                print "Certificate:",tvec
                break
        if red:
            #print "G is f-choosable."
            return True
        else:
            #print "CNS is applicable but inconclusive."
            return False





def CNS_blind_print(G_copy, f, print_mod = 1000):
    begin = time.clock()
    disparity = sum(f)-len(f) - G_copy.size()
    if disparity<0:
        #print "Can't apply CNS to G,f."
        return False
    else:
        count = 0
        red = False
        for fillvec in IntegerListsLex(disparity,length=len(f),ceiling=[f[x]-1 for x in range(len(f))]):
            tvec = [f[x]-1-fillvec[x] for x in range(len(f))]
            count += 1
            if count % print_mod == 0:
                print "%d tvecs computed so far in %s."%(count,timestring(time.clock()-begin))
            if coeff_extract(G_copy,tvec):
                red = True
                print "Certificate:",tvec
                break
        print "%d tvecs computed in %s."%(count,timestring(time.clock()-begin))
        if red:
            #print "G is f-choosable."
            return True
        else:
            #print "CNS is applicable but inconclusive."
            return False





def CNS_blind_QNS(G_copy, f):
    disparity = sum(f)-len(f) - G_copy.size()
    if disparity<0:
        #print "Can't apply CNS to G,f."
        return False
    else:
        red = False
        for fillvec in IntegerListsLex(disparity,length=len(f),ceiling=[f[x]-1 for x in range(len(f))]):
            tvec = [f[x]-1-fillvec[x] for x in range(len(f))]
            if coeff_QNS(G_copy,tvec):
                red = True
                print "Certificate:",tvec
                break
        if red:
            #print "G is f-choosable."
            return True
        else:
            #print "CNS is applicable but inconclusive."
            return False






def QNS_blind_test(G_copy, f, version):
    disparity = sum(f)-len(f) - G_copy.size()
    if disparity<0:
        print "Can't apply CNS to G,f."
    else:
        count = 0
        total = 0
        if version == 0:
            for fillvec in IntegerListsLex(disparity,length=len(f),ceiling=[f[x]-1 for x in range(len(f))]):
                total += 1
                tvec = [f[x]-1-fillvec[x] for x in range(len(f))]
                if coeff_extract(G_copy,tvec):
                    count += 1
        if version == 1:
            for fillvec in IntegerListsLex(disparity,length=len(f),ceiling=[f[x]-1 for x in range(len(f))]):
                total += 1
                tvec = [f[x]-1-fillvec[x] for x in range(len(f))]
                if coeff_QNS1(G_copy,tvec):
                    count += 1
        if version == 2:
            for fillvec in IntegerListsLex(disparity,length=len(f),ceiling=[f[x]-1 for x in range(len(f))]):
                total += 1
                tvec = [f[x]-1-fillvec[x] for x in range(len(f))]
                if coeff_QNS2(G_copy,tvec):
                    count += 1
        if version == 3:
            for fillvec in IntegerListsLex(disparity,length=len(f),ceiling=[f[x]-1 for x in range(len(f))]):
                total += 1
                tvec = [f[x]-1-fillvec[x] for x in range(len(f))]
                if coeff_QNS3(G_copy,tvec):
                    count += 1
        print "Done!  %d out of %d monomials had nonzero coefficient."%(count,total)















#Generates proper colorings of G using the list assignment L(v) := range(tvec[v]+1).
def colorings_and_values(G,tvec):
    n = G.order()
    LA = [set(range(t+1)) for t in tvec]
    CA = [{v for v in range(n) if tvec[v] >= h} for h in range(max(tvec)+1)]

    #Initialize!
    i = 0          #Which vertex we're coloring at the moment.
    c = [-1]*n     #Our coloring function.  c[i] is the current color of vertex i (value -1 means no color assigned).
    stack = []     #Our stack.  Item i is the set of vtxs we deleted color c[i] from (that is, the forward-nbrs which had c[i] in their lists).
    
    partial_vals = [1,]  #Stack of the products of the factors of the graph polynomial which have been induced by the colored vertices.
    #partial_subs = []  #Stack of the subproducts of the normalizing factors of the coefficient formula corresponding the colored vertices.

    remainingColors = None
    
    #At the top of this while loop, we have a partial coloring c.  We have either just colored vertex i-1 (in which case c has colored vertices 0 through i-1 and all of i's remaining colors are available to try), or we have backtracked from a partial coloring that didn't work (in which case we had colored i with something before, and now we need to try another of the remaining colors).
    while True:
        
        if i > n-1: #Full coloring!  Yield, then keep going.
            yield c,partial_vals[-1]
            remainingColors = {x for x in LA[-1] if x>c[-1]}
            while remainingColors: #If we have another color for the last vertex, use it for another coloring.
                c[-1] = min(remainingColors)
                partial_vals[-1] = partial_vals[-2]*int(prod([c[x]-c[-1] for x in G.neighbors(n-1)]))
                yield c,partial_vals[-1]
                remainingColors = {x for x in LA[-1] if x>c[-1]}
            stack.pop()
            backtrack = True
            i -= 1
        else:
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
            partial_vals.append(partial_vals[-1]*int(prod([c[x]-c[-1] for x in G.neighbors(i) if x<i])))
            i+=1

        #Case 2:  There are no colors left for some uncolored vertex, so we pop the stack and backtrack.
        else:
            try:
                oldN = stack.pop()
                for v in oldN:
                    LA[v].add(c[i-1])      #Add the attempted color back to the appropriate vertex lists.
                    CA[c[i-1]].add(v)
                c[i] = -1                #Mark i as uncolored.
                del partial_vals[-1]
                i-=1
            except IndexError:       #If the stack was empty, we've exhausted all possibilities.  Return an empty coloring.
                return




#def eval_poly(G,point):
    #prod = 1
    #for edge in G.edges():
        #prod *= point[edge[0]] - point[edge[1]]
    #return prod














def CNS_initial_power(G_copy, f):
    if sum(f)-len(f) < G_copy.size():
        print "Can't apply CNS to G,f."
        return False
    
    counter = 0
    base = [0]*len(f)
    h = 0
    r = [i for i,x in enumerate(f) if x-1 > h]
    s = G_copy.size()
    while s > len(r):
        print "\ns:",s
        print "base:",base
        print "r:",r
        for i in r:
            base[i] += 1
        h += 1
        s -= len(r)
        r = [i for i,x in enumerate(f) if x-1 > h]
    
    print "\ns:",s
    print "base:",base
    print "r:",r
    print "ceiling:",[f[x]-1-base[x] for x in range(len(f))]
    for fillvec in IntegerListsLex(s,length=len(f),ceiling=[f[x]-1-base[x] for x in range(len(f))]):
        tvec = [base[x]+fillvec[x] for x in range(len(f))]
        print coeff_extract(G_copy,tvec),tvec
        counter += 1
        if counter > 1000:
            return
        if coeff_extract(G_copy,tvec):
            print "G is f-choosable!  Certificate:",tvec
            return True
    print "CNS is inconclusive."
    return False






def CNS_initial_fgap(G_copy, f):
    if sum(f)-len(f) < G_copy.size():
        print "Can't apply CNS to G,f."
        return False
    
    counter = 0
    base = [0]*len(f)
    h = max(f)-2
    r = [i for i,x in enumerate(f) if x-1 > h]
    s = G_copy.size()
    while s > len(r):
        print "\ns:",s
        print "base:",base
        print "r:",r
        for i in r:
            base[i] += 1
        h -= 1
        s -= len(r)
        r = [i for i,x in enumerate(f) if x-1 > h]
    
    print "\ns:",s
    print "base:",base
    print "r:",r
    print "ceiling:",[f[x]-1-base[x] for x in range(len(f))]
    for fillvec in IntegerListsLex(s,length=len(f),ceiling=[f[x]-1-base[x] for x in range(len(f))]):
        tvec = [base[x]+fillvec[x] for x in range(len(f))]
        print coeff_extract(G_copy,tvec),tvec
        counter += 1
        if counter > 1000:
            return
        if coeff_extract(G_copy,tvec):
            print "G is f-choosable!  Certificate:",tvec
            return True
    print "CNS is inconclusive."
    return False






def CNS_initial_degreegap(G_copy, f):
    if sum(f)-len(f) < G_copy.size():
        print "Can't apply CNS to G,f."
        return False
    
    counter = 0
    deg = [G_copy.degree(x) for x in G_copy.vertices()]
    base = [0]*len(f)
    h = max(deg)-1
    r = [i for i,x in enumerate(deg) if x > h and f[i]-1 > base[i]]
    s = G_copy.size()
    while s > len(r):
        print "\ns:",s
        print "base:",base
        print "r:",r
        for i in r:
            base[i] += 1
        h -= 1
        s -= len(r)
        r = [i for i,x in enumerate(deg) if x > h and f[i]-1 > base[i]]
    
    print "\ns:",s
    print "base:",base
    print "r:",r
    print "ceiling:",[f[x]-1-base[x] for x in range(len(f))]
    for fillvec in IntegerListsLex(s,length=len(f),ceiling=[f[x]-1-base[x] for x in range(len(f))]):
        tvec = [base[x]+fillvec[x] for x in range(len(f))]
        print coeff_extract(G_copy,tvec),tvec
        counter += 1
        if counter > 1000:
            return
        if coeff_extract(G_copy,tvec):
            print "G is f-choosable!  Certificate:",tvec
            return True
    print "CNS is inconclusive."
    return False








def CNS_smart(G, f):
    
    #INITIALIZE
    
    #n is the number of vertices.
    n = G.order()
    if n != len(f):
        print "Error:  Length of f does not match order of G."
        return

    if [x for x in f if x<1]:
        #print "Error: f-vector contains non-positive entries."
        return False
    
    #gapleft is the remaining number of steps below the sum of (f-1)-values for our monomial.
    #gapleft begins as the disparity and should drop down to 0.
    #In general, it is given by disparity - sum([f[x]-1-deg[x] for x in range(v)]).
    gapleft = sum(f) - len(f) - G.size()
    if gapleft < 0:
        return False
    
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
    
    G,f = Gf_relabeling(G,f,perminv)
    
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
                    return False
            
            #Drop down the degree (and note this in diff).
            deg[-1] -= 1
            gapleft -= 1
            #Now that we've backtracked, we try moving forward again by finding the next coefficient piece.
        
        
        else:  #If there is hope for expanding the current monomial.
            
            #FORWARD
            
            #If we have built the entire monomial, then stop.
            if v > n-1:
                #print "G is f-choosable.  Certificate:",deg
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
        
        
        #print "time:",time.clock()-begin
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
        








def CNS_smart_all(G, f):
    
    #INITIALIZE
    
    #n is the number of vertices.
    n = G.order()
    if n != len(f):
        print "Error:  Length of f does not match order of G."
        return

    if [x for x in f if x<1]:
        print "Error: f-vector contains non-positive entries."
        return
    
    #gapleft is the remaining number of steps below the sum of (f-1)-values for our monomial.
    #gapleft begins as the disparity and should drop down to 0.
    #In general, it is given by disparity - sum([f[x]-1-deg[x] for x in range(v)]).
    gapleft = sum(f) - len(f) - G.size()
    if gapleft < 0:
        print "Error: negative disparity."
        return
    
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
    
    G,f = Gf_relabeling(G,f,perminv)
    
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
    
    
    count = 0
    total = 0
    
    #START LOOP
    
    while True:
        
        if coeff[-1] == 0 or gapleft < 0 or v > n-1:  #If we have already encountered an issue with expansion.
            
            #If we have built the entire monomial, then backtrack.
            if v > n-1:
                #print "G is f-choosable.  Certificate:",deg
                count += 1
            
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
                    print "Done!  %d monomials had nonzero coefficient out of %d total expansion steps."%(count,total)
                    return
            
            #Drop down the degree (and note this in diff).
            deg[-1] -= 1
            gapleft -= 1
            #Now that we've backtracked, we try moving forward again by finding the next coefficient piece.
        
        
        else:  #If there is hope for expanding the current monomial.
            
            #FORWARD
            
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
        
        
        #print "time:",time.clock()-begin
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
        
        total += 1
        
        #Increment and go back to the top of the loop!
        v += 1
        








def CNS_smart_print_ordering(G, f, ordering):
    begin = time.clock()
    #ma = 0
    #mi = 0
    #coes = set([])
    #count = 0
    
    #zerogapcount = [0,0,0]#coeff=0, gap<0, both
    #abcount = {}
    #for i in range(len(f)):
        #abcount[i] = 0
    
    #INITIALIZE
    
    #n is the number of vertices.
    n = G.order()
    if n != len(f):
        print "Error:  Length of f does not match order of G."
        return

    if [x for x in f if x<1]:
        print "Error: f-vector contains non-positive entries."
        return False
    
    #gapleft is the remaining number of steps below the sum of (f-1)-values for our monomial.
    #gapleft begins as the disparity and should drop down to 0.
    #In general, it is given by disparity - sum([f[x]-1-deg[x] for x in range(v)]).
    gapleft = sum(f) - len(f) - G.size()
    if gapleft < 0:
        print "    >> Disparity:  %d.  Cannot apply Nullstellensatz."%(gapleft)
        return False
    else:
        print "    >> Disparity:  %d.  Proceeding."%(gapleft)
    
    print "    >> Relabeling via",ordering
    G,f = Gf_relabeling(G,f,ordering)
    
    #Polynomial stuff.
    stringlist = ['x'+str(x) for x in range(n)]
    P=PolynomialRing(ZZ, n, stringlist)
    varlist = P.objgens()[1]
    
    #deg, mins, and coeff all act as stacks.
    deg = []    #deg is a list holding the degrees we're currently using.
    mins = []   #mins is a list holding the lower bounds for the degrees we can use in deg.
    coeff = []  #coeff is a list holding the coefficients we've found for varlist[i]**deg[i].
    #Entries in deg start high, but go down by 1 to the corresponding mins value as we backtrack.
    
    #The coeff stack is always one longer than deg and mins; the first entry should be the polynomial 1.
    coeff.append(varlist[0]**0)
    
    #v is the vertex/variable we are inspecting.
    v = 0
    
    
    #START LOOP
    
    while True:
        #count += 1
        #if count %1000 == 0:
            #print count
            #print "\ntime:",time.clock()-begin
            #print "v:",v
            #print "deg:",deg
            #print "mins:",mins
            #print "coeff monomial counts:",[len(x.monomials()) for x in coeff]
            #print "coes:",coes
            #print "gapleft:",gapleft
        
        if coeff[-1] == 0 or gapleft < 0:  #If we have already encountered an issue with expansion.
            #if coeff[-1] ==0:
                #if gapleft < 0:
                    #zerogapcount[2] += 1
                #else:
                    #zerogapcount[0] += 1
            #else:
                #zerogapcount[1] += 1
            
            #print "Backtrack."
            #print "deg:",deg
            #print "coeff[-1]:",coeff[-1]
            #for x in coeff:
                #print x
            #print
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
                #print "gapleft:",gapleft
                del deg[-1]
                del mins[-1]
                
                #If we've marched back through everything, then we've exhausted our options.
                if v < 1:
                    print "    >> The CNS is inconclusive; no relevant monomials have nonzero coefficient."
                    #print "mi: %d    ma: %d"%(mi,ma)
                    #print "Backtracks from zero coefficient: %d.  Backtracks from negative gap: %d.  Backtracks from both: %d."%(zerogapcount[0],zerogapcount[1],zerogapcount[2])
                    #print "Skips from ab:",sum(abcount.values()), [abcount[x] for x in range(len(f))]
                    print "Time:",timestring(time.clock()-begin)
                    return False
            
            #Drop down the degree (and note this in diff).
            deg[-1] -= 1
            gapleft -= 1
            #Now that we've backtracked, we try moving forward again by finding the next coefficient piece.
        
        
        else:  #If there is hope for expanding the current monomial.
            
            #coes = coes | set(coeff[-1].coefficients())
            
            #y = min(coeff[-1].coefficients())
            #z = max(coeff[-1].coefficients())
            #if y < mi:
                #mi = y
            #if z > ma:
                #ma = z
            
            #print "Forward."
            #FORWARD
            
            #If we have built the entire monomial, then stop.
            if v > n-1:
                print "    >> The graph is choosable!  \n    >>Certificate:  The monomial with power list "+str(deg)+" has coefficient "+str(coeff[-1])+"."
                #print "mi: %d    ma: %d"%(mi,ma)
                #print "Backtracks from zero coefficient: %d.  Backtracks from negative gap: %d.  Backtracks from both: %d."%(zerogapcount[0],zerogapcount[1],zerogapcount[2])
                #print "Skips from ab:",sum(abcount.values()), [abcount[x] for x in range(len(f))]
                print "Time:",timestring(time.clock()-begin)
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
                #abcount[v] += f[v]-1-a-b
                deg.append(a+b)
                gapleft -= f[v]-1-a-b
            else:
                deg.append(f[v]-1)
            #Now find the next coefficient piece and try to keep going.
        
        
        #print "time:",time.clock()-begin
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
        








def CNS_smart_ordering(G, f, ordering):
    begin = time.clock()
    #ma = 0
    #mi = 0
    #coes = set([])
    #count = 0
    
    #zerogapcount = [0,0,0]#coeff=0, gap<0, both
    #abcount = {}
    #for i in range(len(f)):
        #abcount[i] = 0
    
    #INITIALIZE
    
    #n is the number of vertices.
    n = G.order()
    if n != len(f):
        print "Error:  Length of f does not match order of G."
        return

    if [x for x in f if x<1]:
        print "Error: f-vector contains non-positive entries."
        return False
    
    #gapleft is the remaining number of steps below the sum of (f-1)-values for our monomial.
    #gapleft begins as the disparity and should drop down to 0.
    #In general, it is given by disparity - sum([f[x]-1-deg[x] for x in range(v)]).
    gapleft = sum(f) - len(f) - G.size()
    if gapleft < 0:
        #print "    >> Disparity:  %d.  Cannot apply Nullstellensatz."%(gapleft)
        return False
    #else:
        #print "    >> Disparity:  %d.  Proceeding."%(gapleft)
    
    #print "    >> Relabeling via",ordering
    G,f = Gf_relabeling(G,f,ordering)
    
    #Polynomial stuff.
    stringlist = ['x'+str(x) for x in range(n)]
    P=PolynomialRing(ZZ, n, stringlist)
    varlist = P.objgens()[1]
    
    #deg, mins, and coeff all act as stacks.
    deg = []    #deg is a list holding the degrees we're currently using.
    mins = []   #mins is a list holding the lower bounds for the degrees we can use in deg.
    coeff = []  #coeff is a list holding the coefficients we've found for varlist[i]**deg[i].
    #Entries in deg start high, but go down by 1 to the corresponding mins value as we backtrack.
    
    #The coeff stack is always one longer than deg and mins; the first entry should be the polynomial 1.
    coeff.append(varlist[0]**0)
    
    #v is the vertex/variable we are inspecting.
    v = 0
    
    
    #START LOOP
    
    while True:
        #count += 1
        #if count %1000 == 0:
            #print count
            #print "\ntime:",time.clock()-begin
            #print "v:",v
            #print "deg:",deg
            #print "mins:",mins
            #print "coeff monomial counts:",[len(x.monomials()) for x in coeff]
            #print "coes:",coes
            #print "gapleft:",gapleft
        
        if coeff[-1] == 0 or gapleft < 0:  #If we have already encountered an issue with expansion.
            #if coeff[-1] ==0:
                #if gapleft < 0:
                    #zerogapcount[2] += 1
                #else:
                    #zerogapcount[0] += 1
            #else:
                #zerogapcount[1] += 1
            
            #print "Backtrack."
            #print "deg:",deg
            #print "coeff[-1]:",coeff[-1]
            #for x in coeff:
                #print x
            #print
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
                #print "gapleft:",gapleft
                del deg[-1]
                del mins[-1]
                
                #If we've marched back through everything, then we've exhausted our options.
                if v < 1:
                    #print "    >> The CNS is inconclusive; no relevant monomials have nonzero coefficient."
                    #print "mi: %d    ma: %d"%(mi,ma)
                    #print "Backtracks from zero coefficient: %d.  Backtracks from negative gap: %d.  Backtracks from both: %d."%(zerogapcount[0],zerogapcount[1],zerogapcount[2])
                    #print "Skips from ab:",sum(abcount.values()), [abcount[x] for x in range(len(f))]
                    #print "Time:",time.clock()-begin
                    #return False
                    return time.clock()-begin
            
            #Drop down the degree (and note this in diff).
            deg[-1] -= 1
            gapleft -= 1
            #Now that we've backtracked, we try moving forward again by finding the next coefficient piece.
        
        
        else:  #If there is hope for expanding the current monomial.
            
            #coes = coes | set(coeff[-1].coefficients())
            
            #y = min(coeff[-1].coefficients())
            #z = max(coeff[-1].coefficients())
            #if y < mi:
                #mi = y
            #if z > ma:
                #ma = z
            
            #print "Forward."
            #FORWARD
            
            #If we have built the entire monomial, then stop.
            if v > n-1:
                #print "    >> The graph is choosable!  \n    >>Certificate:  The monomial with power list "+str(deg)+" has coefficient "+str(coeff[-1])+"."
                #print "mi: %d    ma: %d"%(mi,ma)
                #print "Backtracks from zero coefficient: %d.  Backtracks from negative gap: %d.  Backtracks from both: %d."%(zerogapcount[0],zerogapcount[1],zerogapcount[2])
                #print "Skips from ab:",sum(abcount.values()), [abcount[x] for x in range(len(f))]
                #print "Time:",time.clock()-begin
                #return False
                return time.clock()-begin
            
            #Otherwise, append next values to deg and mins.
            
            #f[v]-1-gapleft is a lower bound for deg[v].  (Any power less than this number 
            #will cause the monomial to already have too low of total degree, even if the 
            #remaining powers are all the f-1 values.)
            mins.append(max(0,f[v]-1-gapleft))
            
            a = coeff[-1].degree(varlist[v])  #a is the degree of xv in the residual coefficient.
            b = len([x for x in G.neighbors(v) if v < x])  #b is the number of forward neighbors of v.
            
            #The monomial we have built so far cannot have a power of xv higher than a+b.
            if a+b < f[v]-1:
                #abcount[v] += f[v]-1-a-b
                deg.append(a+b)
                gapleft -= f[v]-1-a-b
            else:
                deg.append(f[v]-1)
            #Now find the next coefficient piece and try to keep going.
        
        
        #print "time:",time.clock()-begin
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
        








def CNS_smart_print(G, f):
    #begin = time.clock()
    #ma = 0
    #mi = 0
    #coes = set([])
    #count = 0
    
    #zerogapcount = [0,0,0]#coeff=0, gap<0, both
    #abcount = {}
    #for i in range(len(f)):
        #abcount[i] = 0
    
    #INITIALIZE
    
    #n is the number of vertices.
    n = G.order()
    if n != len(f):
        print "Error:  Length of f does not match order of G."
        return

    if [x for x in f if x<1]:
        print "Error: f-vector contains non-positive entries."
        return False
    
    #gapleft is the remaining number of steps below the sum of (f-1)-values for our monomial.
    #gapleft begins as the disparity and should drop down to 0.
    #In general, it is given by disparity - sum([f[x]-1-deg[x] for x in range(v)]).
    gapleft = sum(f) - len(f) - G.size()
    if gapleft < 0:
        print "    >> Disparity:  %d.  Cannot apply Nullstellensatz."%(gapleft)
        return False
    else:
        print "    >> Disparity:  %d.  Proceeding."%(gapleft)
    
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
    
    print "    >> Relabeling via",perminv
    G,f = Gf_relabeling(G,f,perminv)
    
    #Polynomial stuff.
    stringlist = ['x'+str(x) for x in range(n)]
    P=PolynomialRing(ZZ, n, stringlist)
    varlist = P.objgens()[1]
    
    #deg, mins, and coeff all act as stacks.
    deg = []    #deg is a list holding the degrees we're currently using.
    mins = []   #mins is a list holding the lower bounds for the degrees we can use in deg.
    coeff = []  #coeff is a list holding the coefficients we've found for varlist[i]**deg[i].
    #Entries in deg start high, but go down by 1 to the corresponding mins value as we backtrack.
    
    #The coeff stack is always one longer than deg and mins; the first entry should be the polynomial 1.
    coeff.append(varlist[0]**0)
    
    #v is the vertex/variable we are inspecting.
    v = 0
    
    
    #START LOOP
    
    while True:
        #count += 1
        #if count %1000 == 0:
            #print count
            #print "\ntime:",time.clock()-begin
            #print "v:",v
            #print "deg:",deg
            #print "mins:",mins
            #print "coeff monomial counts:",[len(x.monomials()) for x in coeff]
            #print "coes:",coes
            #print "gapleft:",gapleft
        
        if coeff[-1] == 0 or gapleft < 0:  #If we have already encountered an issue with expansion.
            #if coeff[-1] ==0:
                #if gapleft < 0:
                    #zerogapcount[2] += 1
                #else:
                    #zerogapcount[0] += 1
            #else:
                #zerogapcount[1] += 1
            
            #print "Backtrack."
            #print "deg:",deg
            #print "coeff[-1]:",coeff[-1]
            #for x in coeff:
                #print x
            #print
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
                #print "gapleft:",gapleft
                del deg[-1]
                del mins[-1]
                
                #If we've marched back through everything, then we've exhausted our options.
                if v < 1:
                    print "    >> The CNS is inconclusive; no relevant monomials have nonzero coefficient."
                    #print "mi: %d    ma: %d"%(mi,ma)
                    #print "Backtracks from zero coefficient: %d.  Backtracks from negative gap: %d.  Backtracks from both: %d."%(zerogapcount[0],zerogapcount[1],zerogapcount[2])
                    #print "Skips from ab:",sum(abcount.values()), [abcount[x] for x in range(len(f))]
                    return False
            
            #Drop down the degree (and note this in diff).
            deg[-1] -= 1
            gapleft -= 1
            #Now that we've backtracked, we try moving forward again by finding the next coefficient piece.
        
        
        else:  #If there is hope for expanding the current monomial.
            
            #coes = coes | set(coeff[-1].coefficients())
            
            #y = min(coeff[-1].coefficients())
            #z = max(coeff[-1].coefficients())
            #if y < mi:
                #mi = y
            #if z > ma:
                #ma = z
            
            #print "Forward."
            #FORWARD
            
            #If we have built the entire monomial, then stop.
            if v > n-1:
                print "    >> The graph is choosable!  \n    >>Certificate:  The monomial with power list "+str(deg)+" has coefficient "+str(coeff[-1])+"."
                #print "mi: %d    ma: %d"%(mi,ma)
                #print "Backtracks from zero coefficient: %d.  Backtracks from negative gap: %d.  Backtracks from both: %d."%(zerogapcount[0],zerogapcount[1],zerogapcount[2])
                #print "Skips from ab:",sum(abcount.values()), [abcount[x] for x in range(len(f))]
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
                #abcount[v] += f[v]-1-a-b
                deg.append(a+b)
                gapleft -= f[v]-1-a-b
            else:
                deg.append(f[v]-1)
            #Now find the next coefficient piece and try to keep going.
        
        
        #print "time:",time.clock()-begin
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

















def graph_poly(G_copy,override=False):
    if not override and G_copy.size() > 20:
        print "Error:  Graph size arbitrarily limited to 20 edges to prevent computing too large of polynomials."
        return
    
    G = Graph(G_copy.adjacency_matrix())
    stringlist = ['x'+str(v) for v in G.vertices()]
    P=PolynomialRing(ZZ, G.order(), stringlist)

    varlist = P.objgens()[1]
    poly = varlist[0]**0

    for v in G.vertices():
        prod = 1
        for factor in [varlist[v]-varlist[u] for u in G.neighbors(v)]:
            prod *= factor
        G.delete_vertex(v)#change this?
        poly *= prod
    return poly






def eval_poly(G,point):
    prod = 1
    for edge in G.edges():
        prod *= point[edge[0]] - point[edge[1]]
    return prod






















def CNS_print(G_copy, f):
    G = Graph(G_copy)
    s = ''
    for v in G.vertices():
        i=0
        for u in G.neighbors(v):
            i += 1
        s += str(i)+' '
        G.delete_vertex(v)
    
    disparity = sum(f)-len(f) - G_copy.size()
    print 'disparity:',disparity
    if disparity<0:
        #print "Can't apply CNS to G,f."
        return False
    else:
        begin = time.clock()
        count = 0
        red = False
        first = None
        for fillvec in IntegerListsLex(disparity,length=len(f),ceiling=[f[x]-1 for x in range(len(f))]):
            count += 1
            tvec = [f[x]-1-fillvec[x] for x in range(len(f))]
            boo = coeff_extract(G_copy,tvec)
            if boo:
                if not red:
                    lap = time.clock()
                    red = True
                    first = count
                    first_tvec = tvec
                #break
        print 'degree descent (for coeff_extract):',s
        print 'time for all (%d):'%(count),time.clock()-begin
        print
        if red:
            #print "G is f-choosable."
            print 'time until first (%d):'%(first),lap-begin
            print 'first tvec:',first_tvec
            return True
        else:
            #print "CNS is applicable but inconclusive."
            return False



#Generates proper colorings of G using the list assignment L(v) := range(tvec[v]+1).
def colorings(G,tvec):
    n = G.order()
    LA = [set(range(t+1)) for t in tvec]
    CA = [{v for v in range(n) if tvec[v] >= h} for h in range(max(tvec)+1)]

    #Initialize!
    i = 0          #Which vertex we're coloring at the moment.
    c = [-1]*n     #Our coloring function.  c[i] is the current color of vertex i (value -1 means no color assigned).
    stack = []     #Our stack.  Item i is the set of vtxs we deleted color c[i] from (that is, the forward-nbrs which had c[i] in their lists).

    remainingColors = None
    
    #At the top of this while loop, we have a partial coloring c.  We have either just colored vertex i-1 (in which case c has colored vertices 0 through i-1 and all of i's remaining colors are available to try), or we have backtracked from a partial coloring that didn't work (in which case we had colored i with something before, and now we need to try another of the remaining colors).
    while True:
        
        if i > n-1: #Full coloring!  Yield, then keep going.
            yield c
            remainingColors = {x for x in LA[-1] if x>c[-1]}
            while remainingColors: #If we have another color for the last vertex, use it for another coloring.
                c[-1] = min(remainingColors)
                yield c
                remainingColors = {x for x in LA[-1] if x>c[-1]}
            stack.pop()
            backtrack = True
            i -= 1
        else:
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
                return




#Naive calculation
def coeff_QNS1(G,tvec):
    n = G.order()
    total = 0
    for c in colorings(G,tvec):#psi = 1
        prod = Fraction(eval_poly(G,c))
        for j in range(n):
            for x in range(c[j]):
                prod /= c[j] - x
            for x in range(c[j]+1,tvec[j]+1):
                prod /= c[j] - x
        total += prod
    return int(total)




#Multiplies and then divides instead of always divides.
def coeff_QNS2(G,tvec):
    n = G.order()
    total = 0
    for c in colorings(G,tvec):#psi = 1
        prod = Fraction(eval_poly(G,c))
        div = 1
        for j in range(n):
            for x in range(c[j]):
                div *= c[j] - x
            for x in range(c[j]+1,tvec[j]+1):
                div *= c[j] - x
        prod /= div
        total += prod
    return int(total)




#Multiplies and then divides instead of always divides.
def coeff_QNS3(G,tvec):
    subproductsArray = QNSmakeSubproducts(G,tvec)
    n = G.order()
    total = 0
    for c,prod in colorings_and_values(G,tvec):
        div = 1
        for j in range(n):
            div *= subproductsArray[j][c[j]]
        prod = Fraction(prod,div)
        total += prod
    return int(total)




#Multiplies and then divides instead of always divides.
def coeff_QNS4(G,tvec):
    subproductsArray = QNSmakeSubproducts(G,tvec)
    n = G.order()
    total = 0
    for c in colorings(G,tvec):#psi = 1
        prod = Fraction(eval_poly(G,c))
        div = 1
        for j in range(n):
            div *= subproductsArray[j][c[j]]
        prod = Fraction(prod,div)
        total += prod
    return int(total)




def QNSmakeSubproducts(G,tvec):
    array = []
    for j in tvec:
        array.append([])
        for i in range(j+1):
            prod = 1
            for x in range(i):
                prod *= i - x
            for x in range(i+1,j+1):
                prod *= i - x
            array[-1].append(prod)
    return array






def IndSets(graph,maximal=False,complement=False):
    for I in IndependentSets(graph,maximal=maximal,complement=complement):
        yield I








#Returns a colorability class cascade (from a given cascade) of only the subgraphs which are not CNS-reducible when adding them to the current list assignment CA.
#Keeps track of independent sets I that have shown subsets A to be CNS-reducible, and for each new A' first tests to see if I \subseteq A' \subseteq A.
#exclusions is a list of sets, none of which can be contained by any set added to restrictionList.  (This is used for reducers.)
def restriction_list_updater_inprocess(G_copy, ordering, f, CA, given_restrictionList, singletons=True, exclusions=[],maxIS=True):
    
    num_pruned = 0

    restrictionList=[]    #we will consider each item in given_restrictionList and copy it into restrictionList if we can't reduce with it
    for i in range(G_copy.order()):
        restrictionList.append([])

    reducibleList = []

    for vtx in ordering:
        for A in given_restrictionList[vtx]:

            if not singletons and len(A)<2:
                continue

            up = False
            for B in exclusions:
                if B <= A:
                    up = True
            if up:
                continue

            red = False

            #if A is reducible for CA, continue/break to move on to next subset

            #first see if A is reducible by previous observations
            for x in reducibleList:
                if x[0] <= A and A <= x[1]:
                    red = True
                    break
            if red:
                continue

            #now see if A can be shown reducible by CNS
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
                    for I in IndSets(H,maximal=maxIS):#check if the coloring given by 'coloring' and I shows via CNS that A is reducible
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
                            
                            if len(newf) < 1 or CNS_smart(subG,newf):
                                red = True
                                reducibleList.append((set(I),A))
                                break

                            #for fillvec in IntegerListsLex(sum(newf)-len(newf) - subG.size(),length=len(newf),ceiling=[newf[x]-1 for x in range(len(newf))]):
                                #tvec = [newf[x]-1-fillvec[x] for x in range(len(newf))]
                                #if not tvec or coeff_extract(subG,tvec):#If tvec is [], then we've actually just colored all of G.
                                    #red = True
                                    #reducibleList.append((set(I),A))
                                    #break
                            #if red:
                                #break
                    if red:
                        break
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

            if not red:
                restrictionList[vtx].append(A)
            else:
                num_pruned += 1

    return restrictionList,num_pruned








#Printing version of above, without keeping tabs on pruned subgraphs.
def restriction_list_updater_initial(G_copy, ordering, f, CA, given_restrictionList, singletons=True, exclusions=[],maxIS=True):

    restrictionList=[]    #we will consider each item in given_restrictionList and copy it into restrictionList if we can't reduce with it
    for i in range(G_copy.order()):
        restrictionList.append([])

    reducibleList = []
    
    tenpercent = int(sum([len(x) for x in given_restrictionList])/100)+1
    counter = 0
    lap = time.clock()
    begin = lap

    for vtx in ordering:
        for A in given_restrictionList[vtx]:
            if counter %tenpercent==tenpercent-1:
                print "        %d vetted.  Time on batch:  "%(counter)+timestring(time.clock()-lap)
                lap = time.clock()
            counter += 1

            if not singletons and len(A)<2:
                continue

            up = False
            for B in exclusions:
                if B <= A:
                    up = True
            if up:
                continue

            red = False

            #if A is reducible for CA, continue/break to move on to next subset

            #first see if A is reducible by previous observations
            for x in reducibleList:
                if x[0] <= A and A <= x[1]:
                    red = True
                    break
            if red:
                continue

            #now see if A can be shown reducible by CNS
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
                    for I in IndSets(H,maximal=maxIS):#check if the coloring given by 'coloring' and I shows via CNS that A is reducible
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
                            
                            if len(newf) < 1 or CNS_smart(subG,newf):
                                red = True
                                reducibleList.append((set(I),A))
                                break

                            #for fillvec in IntegerListsLex(sum(newf)-len(newf) - subG.size(),length=len(newf),ceiling=[newf[x]-1 for x in range(len(newf))]):
                                #tvec = [newf[x]-1-fillvec[x] for x in range(len(newf))]
                                #if not tvec or coeff_extract(subG,tvec):#If tvec is [], then we've actually just colored all of G.
                                    #red = True
                                    #reducibleList.append((set(I),A))
                                    #break
                            #if red:
                                #break
                    if red:
                        break
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

            if not red:
                restrictionList[vtx].append(A)

    print "        Time vetting all:  "+timestring(time.clock()-begin)
    #for x in restrictionList:
        #for y in x:
            #print y
        #print
    return restrictionList





#Checks whether the subset A is CNS-reducible when adding it to the current list assignment CA.  Returns True if the subset is reducible for the given CA; returns False otherwise (which means it must stay in the restriction list).
def restriction_list_single_check(G_copy, ordering, f, CA, A, maxIS=True):
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
            for I in IndSets(H,maximal=maxIS):#check if the coloring given by 'coloring' and I shows via CNS that A is reducible
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

                    if len(newf) < 1 or CNS_smart(subG,newf):
                        return True

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
















#Hijacked for non-Minimal, but not currently leveraging any Minial properties.
#Please have the vertices be named 0,...,(n-1).
def fChoosableMinimalAll(G_copy,f_copy,inprocess=3,print_mod=1000,maxIS=True,LAscheme="L",rev=False):
    ordering_copy = range(G_copy.order())
    
    begin = time.clock()

    it = Subsets1(G_copy.vertices())
    for i in range(1+G_copy.order()+G_copy.size()):
        it.next()
    for a in it:
        subG = G_copy.subgraph(a)
        if not subG.is_connected():
            continue
        subordering = [ordering_copy[x] for x in range(len(ordering_copy)) if ordering_copy[x] in a]
        subf = [f_copy[ordering_copy.index(x)] for x in subordering]
        skip = False
        for v in a:
            if subf[subordering.index(v)] > subG.degree(v):
                skip=True
                break
        if skip:
            continue
        subG.relabel()
        print "\na:",a
        z = fChoosableReturnBad(subG,subf,inprocess=inprocess,print_mod=print_mod,maxIS=maxIS,LAscheme=LAscheme,rev=rev)
        print "Total time so far:",timestring(time.clock()-begin)
        if z:
            return False
        #z = fChoosableMinimalSubgraph(subG,subf,inprocess)
        #if not z[0]:
            #return False,z[1]
    return True,[]













def fChoosableMinimalSubgraph(G,f,inprocess=3):

    #Trying to catch input error.
    if len(f) != G.order():
        print "Error: f-vector length does not match graph order."
        return


    it = restrictedListAssignmentsMinimal(G,f,inprocess_depth=inprocess)

    try:
        x = it.next()           #The first call to a Python iterator must be .next().
        #^^Recall: x is a 2-tuple (CA,E,L), where CA is a colorability class assignment, E is the list of vertices whose lists aren't full, and L is the list of the sizes of the current lists of the vertices.

        #At the top of this loop, we have obtained an x from the listAssignment iterator.  We check its representability -- unless it is representable, uncolorable, and full, we take the next x from the iterator.
        while True:

            #If any full super-list assignment requires at least as many colors as vertices, we can turn back right now.
            if len(f) - len(x[0]) - max([f[v]-x[2][v] for v in range(len(f))]) < 1:
                b = []
                c = []

            else:

                b = assignmentRepresentationCheck(x[0])

                #If the list assignment is representable, we need to check with certainty whether it is colorable.
                if b:
                    c = assignmentCheckInduced(G,x[0])

                    #If the list assignment is not colorable and full, then stop -- G is not f-choosable.
                    if not c and not x[1]:
                        return False,x[0]

                #If the list assignment is not representable, we know it is colorable but don't want to take the time to find the coloring.
                else:
                    c = []

            #Tell the listAssignment iterator to go to the next relevant list assignment (either moving forward or backtracking).
            send_value = (b,c)
            x = it.send(send_value)

    except StopIteration:           #If we never find an uncolorable full list, then any f-list-assignment is colorable.
        return True,None









def Subsets1(given_list):
    for x in Subsets(given_list):
        yield x
    return









#Takes a colorability class assignment as a list of sets, returns either a list of distinct representatives for the colors or [] if none exist.
#Here, a representative for a color must have that color in its list.
def assignmentRepresentationCheck(CA):
    d = []              #d is the list of representatives.  d[j] is the vertex representing color j.
    col = 0             #col is the color we are about to choose a representative for.
    U = set([])         #U is the set of representatives we are currently using.
    remaining = copy.deepcopy(CA)   #remaining is a list of sets of the vertices which have yet to be tried as a representative in this branch of the space.

    while col < len(CA):
        up = False

        #Move forward and add a representative, if possible.
        for v in remaining[col]-U:#This loop will executive once if the set remaining[col]-U is nonempty, not at all if it is empty.
            d.append(v)
            U.add(v)
            col +=1
            up = True
            break

        #If we moved forward, go back to the top of the loop.
        if up:
            continue

        #Otherwise, backtrack.
        remaining[col:] = copy.deepcopy(CA[col:])
        col -=1
        if col < 0:#In this case, d is empty.
            break
        U.remove(d[col])
        remaining[col].remove(d[col])
        del d[col]

    return d






def assignmentCheckInduced(G_big,CAtocopy):
    #Convert the colorability classes assignment to a list assignment.
    A = set([])
    for a in CAtocopy:
        A = A | a
    phi = dict()
    for v in A:
        phi[v] = len([x for x in A if x<v])
    G = G_big.subgraph(A)
    G.relabel()
    CA = []
    for a in CAtocopy:
        CA.append({phi[x] for x in a})
    #^^this is the prep for the induced stuff.

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













#inprocess_depth controls how many partial list assignments are checked for the reducibility when extending during the normal routine.
def restrictedListAssignmentsMinimal(G,f,inprocess_depth=3):
    #moved from parameters
    ordering = range(G.order())
    reducers=[]

    #Initialize the things!
    generator_stack = []
    colClasses = []                                    #colClasses[i] is {vertices with color i in their list}.  We yield this.
    L = [0]*G.order()                                  #Keeps track of how many colors each vertex is currently assigned.
    D = None
    E = [x for x in ordering if L[x]<f[x]]             #Eligible vertices for colorability classes.
    col = 0                                            #Keeps track of next colorability class('s color).
    ticket = None

    reducers = [set(x) for x in reducers]

    initial_restriction_list_dict = tuple_restriction_list_maker(G, ordering, f, max_depth=0, exclusions=reducers)

    max_mult = {}
    current_mult = {}
    for vtxList in initial_restriction_list_dict[()]:
        for A in vtxList:
            max_mult[set2tup(A)] = multiplicity_checker(G,ordering,f,A)
            current_mult[set2tup(A)] = 0

    #First round, to make generator_stack nonempty.
    it = RestrictedConnectedSubgraphs(G,E,initial_restriction_list_dict[()])
    generator_stack.append(it)                         #Stack of subgraph iterators; cth element generates the c-colorability class.

    if it.next_within(L,E):
        A = generator_stack[col].subgraph
        colClasses.append(A)
        current_mult[set2tup(A)] += 1
        for v in colClasses[col]:
            L[v]+=1
        E = [x for x in ordering if L[x]<f[x]]
    else:
        return


    #Begin generation!
    while generator_stack:
        #Yield the current assignment; receive information from caller.
        ticket = yield (colClasses,E,L)
        #ticket is a tuple of the form (list, list).
        #ticket[0] is a list of distinct representatives for colClasses, or [] if colClasses is not representable.
        #ticket[1] will be either a coloring for colClasses or else [] (uncolorable assignment or no coloring attempted).

        #Case 1:  Move forward.  (Build the current partial list assignment up one more level.)
        #Since it is a new color, we start/copy a new iterator and add it to the stack (rather than call one we already have).
        #(There is no use for D in this case.)
        if ticket[0] and not ticket[1]:#Case:  The partial list assignment was representable, then we found it was not colorable.
            #Add the new colorability class iterator to the iterator stack, then call the first subgraph, add it to the stack, and update.
            generator_stack.append(generator_stack[col].copy())  #Add new subgraph iterator (for the color col+1) to the stack.
            col+=1                                               #Increment col -- now the class we are generating is for the color col.

            #Because we've just made a copy of the top generator, we need to check the current colorability class for repetition:
            A = generator_stack[col].subgraph
            if current_mult[set2tup(A)] < max_mult[set2tup(A)] and A <= set(E):
                if len(colClasses) < inprocess_depth:
                    if not restriction_list_single_check(G, ordering, f, colClasses, A):
                        colClasses.append(A)
                        current_mult[set2tup(A)] += 1
                        for v in colClasses[col]:                            #Update L.
                            L[v]+=1
                        E = [x for x in ordering if L[x]<f[x]]               #Update E.
                        continue
                else:#If we have more than a few colors and we know A fits, let's just try it!
                    colClasses.append(A)
                    current_mult[set2tup(A)] += 1
                    for v in colClasses[col]:                            #Update L.
                        L[v]+=1
                    E = [x for x in ordering if L[x]<f[x]]               #Update E.
                    continue
            #If the current subgraph doesn't fit into E, then we check for the next eligible subgraph after updating restrIction lists.
            #We don't update the restriction lists when backtracking since they've already been updated here first.

            #(Update the restrictionList from the current colClasses if we are in the early stages.)
            if len(colClasses) < inprocess_depth:
                generator_stack[col].update_restriction_list(f,colClasses)
            
            #Now go and get the next subgraph.
            if generator_stack[col].next_within(L,E):
                A = generator_stack[col].subgraph
                colClasses.append(A)
                current_mult[set2tup(A)] += 1
                for v in colClasses[col]:                            #Update L.
                    L[v]+=1
                E = [x for x in ordering if L[x]<f[x]]               #Update E.
                continue
            #If there is no such subgraph, then we are as full as possible and need to backtrack.

        #Case 2a:  Backtrack because we tried to fill the current list assignment, but there were no more valid subsets.
        #Skip over the following else portion down to the main backtracking.


        #Case 2b:  Backtrack because the current list assignment colClasses is colorable (on the vertices colClasses sees).
        #This could be because colClasses is not representable (in which case ticket[0] is []) or because we found a coloring (ticket[1]).
        #Start by moving the top color forward with respect to not containing the current colorability class.
        else:
            #Delete the top colorability class (but not the generator for the class quite yet).
            for v in colClasses[col]:
                L[v] -=1
            current_mult[set2tup(colClasses[col])] -= 1
            del colClasses[col]
            E = [x for x in ordering if L[x]<f[x]]

            #If the provided coloring (or lack of a provided coloring) is not on all the vertices, then we can't skip ahead in any clever way.
            if len(ticket[1]) < G.order():
                #First, see if the generator can give us something new:
                if generator_stack[col].next_within(L,E):
                    A = generator_stack[col].subgraph
                    colClasses.append(A)
                    current_mult[set2tup(A)] += 1
                    for v in colClasses[col]:
                        L[v] +=1
                    E = [x for x in ordering if L[x]<f[x]]
                    continue
                #If not, then we need to peel off more generators/layers.

            #But if the provided coloring is on all vertices, then we can use the top color class to potentially skip past some colorability class candidates.
            else:
                D = [x for x in colClasses[-1] if ticket[1][x]==col]  #D is the set of vertices in the coloring's highest color class.

                #First, see if the generator can give us something new:
                if generator_stack[col].next_within_not_containing(L,E,D):
                    A = generator_stack[col].subgraph
                    colClasses.append(A)
                    current_mult[set2tup(A)] += 1
                    for v in colClasses[col]:
                        L[v] +=1
                    E = [x for x in ordering if L[x]<f[x]]
                    continue
                #If not, then we need to peel off more generators/layers.  Note: D no longer serves any purpose past here.


        #At the moment, generator_stack[col] is an iterator that has no valid subgraphs left (either because it was a false start in Case 2a or because we backtracked through it in Case 2b).  But colClasses[col] is undefined (either because it didn't get there in Case 2a or because we deleted it in Case 2b).
        #Case 2:  Backtrack.  Step back as many colors as necessary and then obtain a next subgraph -- add that to colClasses.

#Is this code block redundant?  Do we need to avoid calling next_within in the first step?
        #At this point, we are guaranteed to take at least one step back through our stacks.
        del generator_stack[col]
        col -= 1
        if col<0:
            return
        for v in colClasses[col]:
            L[v] -=1
        current_mult[set2tup(colClasses[col])] -= 1
        del colClasses[col]
        E = [x for x in ordering if L[x]<f[x]]

        #Now step back as many times as necessary.
        while not generator_stack[col].next_within(L,E):

            del generator_stack[col]
            col -= 1
            if col<0:#If generator_stack is empty, then we have backtracked entirely.  A StopIteration error will be thrown.
                return
            #Prepare for bumping the now-top generator on the stack.
            for v in colClasses[col]:
                L[v] -=1
            current_mult[set2tup(colClasses[col])] -= 1
            del colClasses[col]
            E = [x for x in ordering if L[x]<f[x]]

        #Now we have successfully bumped the generator to another subgraph.  Push the subgraph onto the colClass stack.
        A = generator_stack[col].subgraph
        colClasses.append(A)
        current_mult[set2tup(A)] += 1
        for v in colClasses[col]:
            L[v]+=1
        E = [x for x in ordering if L[x]<f[x]]
        continue

















#Takes a subset A and returns the maximum multiplicity for which it is not CNS-reducible.  (A is itself assumed to not be CNS-reducible, and we know that the multiplicity must be at most |A|.  When G-copy is unchoosable-minimal, we know that the multiplicity will be less than |A|.)
def multiplicity_checker(G_copy, ordering, f, A, minimal=True,maxIS=True):

    CA = [A,]

    if minimal:
        cap = min(min([f[x] for x in A]),len(A)-1)
    else:
        cap = min([f[x] for x in A])#used to be min(min([f[x] for x in A]),len(A)), but I can't remember why -- I think we thought the connected small pot lemma was true.  Since f-values are at most 7 for wegner configurations anyway, this should be okay.  It might be worth considering putting in a ceiling (less than cap) of how far we check for large subgraphs with high f-values since this could take a lot of time that might not be necessary.  I mean, what's the practical difference between a max multiplicity of 5 and a max multiplicity of 10?
    while len(CA) < cap:

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
                for I in IndSets(H,maximal=maxIS):#check if the coloring given by 'coloring' and I shows via CNS that A is reducible
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
                        
                        if len(newf) < 1 or CNS_smart(subG,newf):
                            return len(CA)

                        #for fillvec in IntegerListsLex(sum(newf)-len(newf) - subG.size(),length=len(newf),ceiling=[newf[x]-1 for x in range(len(newf))]):
                            #tvec = [newf[x]-1-fillvec[x] for x in range(len(newf))]
                            #if not tvec or coeff_extract(subG,tvec):#If tvec is [], then we've actually just colored all of G.
                                #return len(CA)

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
                if c==0:#If we've searched all possibilities for CA.append(A)
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
            ISgens.append(IndSets(H,maximal=True))#since we want to properly color as many vertices as we can, we use maximal
            new = set(ISgens[c].next())
            if c==0:
                coloring.append((new,new))
            else:
                coloring.append((new,new|coloring[c-1][1]))
            c+=1

        CA.append(A)

    return cap










#Returns a dictionary which has partial list assignments as keys (in the form of colorability class assignments as tuples of tuples) with each corresponding value being a colorability class cascade (list of lists of sets) containing the subsets which are CNS-reducible for the partial list assignment.  The partial list assignments are only those which use non-reducible subsets, and partial list assignments with repeated subsets are considered.  (Even though multiplicities are currently handled separately in the iterator, it is possible to backtrack to a partial list assignment which has repreated subsets.)  max_depth is the maximum length of the partial list assignments we'll consider.
def initial_restriction_list_maker(G_copy, ordering, f, maxIS=True):
    
    CA = []#CA is from the inprocess version of this routine.  Here, it can be removed entirely.
    
    begin = time.clock()
    
    count = 0

    for vtx in ordering:
        subsets = RootedConnectedSubgraphs(G_copy,ordering[ordering.index(vtx):],vtx)
        subsets.first()
        while subsets.next():
            count += 1
    
    lap = time.clock()
    
    print "    >> %d connected subgraphs of our graph.  (Took %s to count them.)  Now vetting them all."%(count,timestring(lap-begin))
    

    restrictionList=[]    #we will consider each item in given_restrictionList and copy it into restrictionList if we can't reduce with it
    for i in range(G_copy.order()):
        restrictionList.append([])

    reducibleList = []
    
    tenpercent = int(count/100)+1
    counter = 0
    lap = time.clock()
    begin = lap

    for vtx in ordering:
        subsets = RootedConnectedSubgraphs(G_copy,ordering[ordering.index(vtx):],vtx)
        subsets.first()
        while subsets.next():
            A = subsets.subgraph
            
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
                continue

            #now see if A can be shown reducible by CNS
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
                    for I in IndSets(H,maximal=maxIS):#check if the coloring given by 'coloring' and I shows via CNS that A is reducible
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
                            
                            if len(newf) < 1 or CNS_smart(subG,newf):
                                red = True
                                reducibleList.append((set(I),A))
                                break

                            #for fillvec in IntegerListsLex(sum(newf)-len(newf) - subG.size(),length=len(newf),ceiling=[newf[x]-1 for x in range(len(newf))]):
                                #tvec = [newf[x]-1-fillvec[x] for x in range(len(newf))]
                                #if not tvec or coeff_extract(subG,tvec):#If tvec is [], then we've actually just colored all of G.
                                    #red = True
                                    #reducibleList.append((set(I),A))
                                    #break
                            #if red:
                                #break
                    if red:
                        break
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

            if not red:
                restrictionList[vtx].append(A)

    print "        Time vetting all:  "+timestring(time.clock()-begin)
    #for x in restrictionList:
        #for y in x:
            #print y
        #print
    #print len(reducibleList)
    return restrictionList










#Returns a dictionary which has partial list assignments as keys (in the form of colorability class assignments as tuples of tuples) with each corresponding value being a colorability class cascade (list of lists of sets) containing the subsets which are CNS-reducible for the partial list assignment.  The partial list assignments are only those which use non-reducible subsets, and partial list assignments with repeated subsets are considered.  (Even though multiplicities are currently handled separately in the iterator, it is possible to backtrack to a partial list assignment which has repreated subsets.)  max_depth is the maximum length of the partial list assignments we'll consider.
def initial_restriction_list_makerNoPrint(G_copy, ordering, f, maxIS=True):
    
    CA = []#CA is from the inprocess version of this routine.  Here, it can be removed entirely.
    
    #begin = time.clock()
    
    #count = 0

    #for vtx in ordering:
        #subsets = RootedConnectedSubgraphs(G_copy,ordering[ordering.index(vtx):],vtx)
        #subsets.first()
        #while subsets.next():
            #count += 1
    
    #lap = time.clock()
    
    #print "    >> %d connected subgraphs of our graph.  (Took %s to count them.)  Now vetting them all."%(count,timestring(lap-begin))
    

    restrictionList=[]    #we will consider each item in given_restrictionList and copy it into restrictionList if we can't reduce with it
    for i in range(G_copy.order()):
        restrictionList.append([])

    reducibleList = []
    
    #tenpercent = int(count/100)+1
    #counter = 0
    #lap = time.clock()
    #begin = lap

    for vtx in ordering:
        subsets = RootedConnectedSubgraphs(G_copy,ordering[ordering.index(vtx):],vtx)
        subsets.first()
        while subsets.next():
            A = subsets.subgraph
            
            #if counter %tenpercent==tenpercent-1:
                #print "        %d vetted.  Time on batch:  "%(counter)+timestring(time.clock()-lap)
                #lap = time.clock()
            #counter += 1

            red = False

            #if A is reducible for CA, continue/break to move on to next subset

            #first see if A is reducible by previous observations
            for x in reducibleList:
                if x[0] <= A and A <= x[1]:
                    red = True
                    break
            if red:
                continue

            #now see if A can be shown reducible by CNS
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
                    for I in IndSets(H,maximal=maxIS):#check if the coloring given by 'coloring' and I shows via CNS that A is reducible
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
                            
                            if len(newf) < 1 or CNS_smart(subG,newf):
                                red = True
                                reducibleList.append((set(I),A))
                                break

                            #for fillvec in IntegerListsLex(sum(newf)-len(newf) - subG.size(),length=len(newf),ceiling=[newf[x]-1 for x in range(len(newf))]):
                                #tvec = [newf[x]-1-fillvec[x] for x in range(len(newf))]
                                #if not tvec or coeff_extract(subG,tvec):#If tvec is [], then we've actually just colored all of G.
                                    #red = True
                                    #reducibleList.append((set(I),A))
                                    #break
                            #if red:
                                #break
                    if red:
                        break
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

            if not red:
                restrictionList[vtx].append(A)

    #print "        Time vetting all:  "+timestring(time.clock()-begin)
    #for x in restrictionList:
        #for y in x:
            #print y
        #print
    #print len(reducibleList)
    return restrictionList











#Returns a dictionary which has partial list assignments as keys (in the form of colorability class assignments as tuples of tuples) with each corresponding value being a colorability class cascade (list of lists of sets) containing the subsets which are CNS-reducible for the partial list assignment.  The partial list assignments are only those which use non-reducible subsets, and partial list assignments with repeated subsets are considered.  (Even though multiplicities are currently handled separately in the iterator, it is possible to backtrack to a partial list assignment which has repreated subsets.)  max_depth is the maximum length of the partial list assignments we'll consider.
def tuple_restriction_list_maker(G, ordering, f, minimal=True, max_depth=1, exclusions=[],maxIS=True):
    
    begin = time.clock()

    restriction_dict = {}

    fullCascade=[]
    for i in range(G.order()):
        fullCascade.append([])
    for vtx in ordering:
        subsets = RootedConnectedSubgraphs(G,ordering[ordering.index(vtx):],vtx)
        subsets.first()
        while subsets.next():
            A = subsets.subgraph
            fullCascade[vtx].append(A)
            
    print "    >>",sum([len(x) for x in fullCascade]),"connected subgraphs of our graph.  Now vetting them all."

    if minimal:#If G is unchoosable-minimal, then we know that singletons cannot be part of a bad list assignment.
        restriction_dict[()] = restriction_list_updater_initial(G, ordering, f, [], fullCascade, singletons=False,exclusions=exclusions,maxIS=maxIS)
    else:
        restriction_dict[()] = restriction_list_updater_initial(G, ordering, f, [], fullCascade,exclusions=exclusions,maxIS=maxIS)

    #PLAs = [[[],],] #PLAs keeps track of all the partial list assignments we might use.

    #while len(PLAs) < max_depth:

        #PLAs.append([])

        #for CA in PLAs[-2]:
            #CA_tuple = tuple(set2tup(x) for x in CA)

            #if CA:
                #a = ordering.index([x for x in ordering if x in CA[-1]][0])
            #else:
                #a = 0
            #for v in ordering[a:]:
                #new_rList = copy.deepcopy(restriction_dict[CA_tuple])
                #for A in restriction_dict[CA_tuple][v]:
                    #if CA or ordering[0] in A:
                        #CA.append(A)
                        #cascade = restriction_list_updater(G, ordering, f, CA, new_rList)
                        #restriction_dict[CA_tuple+(set2tup(A),)] = cascade
                        #for L in cascade: #If the cascade is empty, then we don't need to add CA to PLAs.
                            #if L:
                                #PLAs[-1].append(CA[:])
                                #break
                        #del CA[-1]
                    ##The only trim we can make on the restriction list passed to restriction_list_updater is in each A's column.
                    #del new_rList[v][0]

    return restriction_dict











#-----------------------------------------------------------------------
#Stuff for "Par" versions






def fChoosableMinimalAllPar(G_copy,f_copy,inprocess=3,m=100,r=0):
    begin = time.clock()
    lap = time.clock()
    
    ordering_copy = range(G_copy.order())

    it = Subsets1(G_copy.vertices())
    for i in range(1+G_copy.order()+G_copy.size()):#This should be fixed sometime to reflect a more intelligent means of skipping singletons and edges, etc., if we want to do so.  Right now it does not skip all edges and any graphs we work with do not have adjacent 1s, so we're fine.
        it.next()
    counter=0
    for a in it:
        if counter%m==r:
            print "\nResidue Class %d, Task #%d:"%(r,int(counter/m)),a
            subG = G_copy.subgraph(a)
            subordering = [ordering_copy[x] for x in range(len(ordering_copy)) if ordering_copy[x] in a]
            subf = [f_copy[ordering_copy.index(x)] for x in subordering]
            skip = False
            for v in a:
                if subf[subordering.index(v)] > subG.degree(v):
                    skip=True
                    break
            if skip:
                print "Skip!  (Degree exceeded.)"
                counter += 1
                continue
            subG.relabel()
            if not subG.is_connected():
                print "Skip!  (Disconnected.)"
                counter += 1
                continue
            print "Running subroutine."
            z = fChoosableMinimalSubgraphReturnBad(subG,subf,inprocess)
            if z[0]:#For the ReturnBad subroutine, z[0] is the list of bad LAs, if there are any.
                print "Bad LA!  Stopping."
                print "Lap time:  "+timestring(time.clock()-lap)
                print "Total time:  "+timestring(time.clock()-begin)
                return False
            else:
                print "Good!"
            print "Lap time:  "+timestring(time.clock()-lap)
            print "Total time:  "+timestring(time.clock()-begin)
            lap = time.clock()
        counter += 1

    print "\nDone!  No bad LAs in this residue class."
    print "Lap time:  "+timestring(time.clock()-lap)
    print "Total time:  "+timestring(time.clock()-begin)
    return True










#d is depth of the tree, m is the modulus, r is the residue for the job.
#d should be at least 2 -- if it is 1, we need to put in a residue check for the first it.next() call.
def fChoosableMinimalSubgraphParTreeCount(G,f,inprocess=3,d=3,count_modulus=1000):
    
    begin = time.clock()
    print "Starting up!"
    
    top_node = -1  #Indexes which node at depth d we're currently at.
    explore = False  #boolean which indicates to search all the descendants of the current node.

    #Trying to catch input error.
    if len(f) != G.order():
        print "Error: f-vector length does not match graph order."
        return


    it = restrictedListAssignmentsMinimal(G,f,inprocess_depth=inprocess)

    try:
        x = it.next()           #The first call to a Python iterator must be .next().
        #^^Recall: x is a 2-tuple (CA,E,L), where CA is a colorability class assignment, E is the list of vertices whose lists aren't full, and L is the list of the sizes of the current lists of the vertices.

        #At the top of this loop, we have obtained an x from the listAssignment iterator.  We check its representability -- unless it is representable, uncolorable, and full, we take the next x from the iterator.
        while True:

            #If any full super-list assignment requires at least as many colors as vertices, we can turn back right now.
            if len(f) - len(x[0]) - max([f[v]-x[2][v] for v in range(len(f))]) < 1:
                b = []
                c = []

            else:

                b = assignmentRepresentationCheck(x[0])

                #If the list assignment is representable, we need to check with certainty whether it is colorable.
                if b:
                    c = assignmentCheckInduced(G,x[0])

                    #If the list assignment is not colorable and full, then stop -- G is not f-choosable.
                    if not c and not x[1]:
                        print "\n Oh, no! (Bad LA found)"
                        print " "*10+"*"*10+" BAD LIST ASSIGNMENT:"+"*"*10
                        for i in range(len(x[0])):
                            print " "*10+subsetstring(x[0][i],G.order())
                        print "\nStopped after "+timestring(time.clock()-begin)
                        return False,x[0]

                #If the list assignment is not representable, we know it is colorable but don't want to take the time to find the coloring.
                else:
                    c = []

            #Here's the new stuff for parallelization:
            #The generator moves forward iff b = True and c = False.
            #We want to move forward only if the LA was not colorable AND we fall into the correct residue class.
            if b and not c and len(x[0])==d:  #Nodes at depth d and which want to move forward are precisely the tasks we are splitting up.
                top_node += 1
                if top_node % count_modulus == 0:
                    print "Counted %d nodes so far - it has taken %s."%(top_node,timestring(time.clock()-begin))
                b = []
            
            #Tell the listAssignment iterator to go to the next relevant list assignment (either moving forward or backtracking).
            send_value = (b,c)
            
            x = it.send(send_value)

    except StopIteration:           #If we never find an uncolorable full list, then any f-list-assignment is colorable.
        print "Total of %d nodes to explore (having %d colorability classes).  Time to go through those layers of the tree:  %s"%(top_node+1,d,timestring(time.clock()-begin))










#d is depth of the tree, m is the modulus, r is the residue for the job.
#d should be at least 2 -- if it is 1, we need to put in a residue check for the first it.next() call.
def fChoosableMinimalSubgraphPar(G,f,inprocess=3,d=3,m=100,r=1):
    
    begin = time.clock()    
    lap = time.clock()
    print "Starting up!"
    
    top_node = -1  #Indexes which node at depth d we're currently at.
    explore = False  #boolean which indicates to search all the descendants of the current node.

    #Trying to catch input error.
    if len(f) != G.order():
        print "Error: f-vector length does not match graph order."
        return


    it = restrictedListAssignmentsMinimal(G,f,inprocess_depth=inprocess)

    try:
        x = it.next()           #The first call to a Python iterator must be .next().
        #^^Recall: x is a 2-tuple (CA,E,L), where CA is a colorability class assignment, E is the list of vertices whose lists aren't full, and L is the list of the sizes of the current lists of the vertices.

        #At the top of this loop, we have obtained an x from the listAssignment iterator.  We check its representability -- unless it is representable, uncolorable, and full, we take the next x from the iterator.
        while True:

            #If any full super-list assignment requires at least as many colors as vertices, we can turn back right now.
            if len(f) - len(x[0]) - max([f[v]-x[2][v] for v in range(len(f))]) < 1:
                b = []
                c = []

            else:

                b = assignmentRepresentationCheck(x[0])

                #If the list assignment is representable, we need to check with certainty whether it is colorable.
                if b:
                    c = assignmentCheckInduced(G,x[0])

                    #If the list assignment is not colorable and full, then stop -- G is not f-choosable.
                    if not c and not x[1]:
                        print "\n Oh, no! (Bad LA found)"
                        print " "*10+"*"*10+" BAD LIST ASSIGNMENT:"+"*"*10
                        for i in range(len(x[0])):
                            print " "*10+subsetstring(x[0][i],G.order())
                        print "\nStopped with this residue class after a total time of "+timestring(time.clock()-begin)
                        return False,x[0]

                #If the list assignment is not representable, we know it is colorable but don't want to take the time to find the coloring.
                else:
                    c = []

            #Here's the new stuff for parallelization:
            #The generator moves forward iff b = True and c = False.
            #We want to move forward only if the LA was not colorable AND we fall into the correct residue class.
            if b and not c and len(x[0])==d:  #Nodes at depth d and which want to move forward are precisely the tasks we are splitting up.
                top_node += 1
                if top_node % m != r:  #If we aren't in the correct residue class, backtrack instead of explore.
                    b = []
                else:
                    print "Lap Time (including in-between nodes):  %s"%(timestring(time.clock()-lap))
                    lap = time.clock()
                    print "\nResidue class:  %d.  Current task #:  %d."%(r,int(top_node/m))
            
            #Tell the listAssignment iterator to go to the next relevant list assignment (either moving forward or backtracking).
            send_value = (b,c)
            
            x = it.send(send_value)

    except StopIteration:           #If we never find an uncolorable full list, then any f-list-assignment is colorable.
        print "Lap Time (including in-between nodes):  %s"%(timestring(time.clock()-lap))
        print "\n No bad LAs found!"
        print "Finished with this residue class in a total time of "+timestring(time.clock()-begin)
        return True,None


















#-----------------------------------------------------------------------
#Stuff for "ReturnBad" versions



#Please have the vertices be named 0,...,(n-1).
def fChoosableMinimalAllReturnBad(G_copy,f_copy,inprocess=3,modulus=1000,name="G"):
    begin = time.clock()
#     if not ordering_copy:
#         n = G_copy.order()
    ordering_copy = range(G_copy.order())
    
    
    #
    
    
    
    masterboo = True

    print "Subgraphs not listed are ruled out immediately via degree reduction or disconnectedness."

    counter = 0

    it = Subsets1(G_copy.vertices())
    for i in range(1+G_copy.order()+G_copy.size()):
        it.next()
        counter += 1
    for a in it:
        counter += 1
        subG = G_copy.subgraph(a)
        subordering = [ordering_copy[x] for x in range(len(ordering_copy)) if ordering_copy[x] in a]
        subf = [f_copy[ordering_copy.index(x)] for x in subordering]
        skip = False
        for v in a:
            if subf[subordering.index(v)] > subG.degree(v):
                skip=True
                break
        if skip:
            continue
        if not subG.is_connected():
            continue
        subG.relabel()
        print "#%d"%(counter),a
        if CNS(subG,subf):
            print "> CNS applies!  Total time so far: "+timestring(time.clock()-begin)+"\n"
            continue
        else:
            print "> Initiating brute force.  Total time so far: "+timestring(time.clock()-begin)
        z,N = fChoosableMinimalSubgraphReturnBad(subG,subf,inprocess,modulus)
        if z:
            masterboo = False
            print a,"\n",N,"total list assignments directly examined.",len(z),"are bad list assignments.\n"
            count = 1
            fi = open(name+"-"+str(preprocessing)+str(inprocess)+str(extra)+"-store.txt","w+")
            fi.write(str(a)+"\n\n")
            for x in z:
                print "Bad LA #"+str(count)+" ("+str(x[1])+") ("+str(x[2])+"):\n",x[0],"\n"
                fi.write(str(x[0])+"\n\n")
                count += 1
            fi.close()
        else:
            print "> No bad list assignments!"
        print "> Total time so far: "+timestring(time.clock()-begin)+"\n\n"
    print "> Done!  Took",timestring(time.clock()-begin),"total."
    if masterboo:
        print "The graph is f-choosable."
    else:
        print "The graph is NOT f-choosable."











def fChoosableMinimalSubgraphReturnBad(G,f,inprocess=3,modulus=1000):

    start = time.clock()
    lap = start

    if len(f) != G.order():
        print "Error: f-vector length does not match graph order."
        return


    it = restrictedListAssignmentsMinimal(G,f,inprocess_depth=inprocess)

    z = []

    nodes_count = 0
    partial_count = 0
    first_count = 0
    unrep_count = 0
    bad_nums = []
    node_code = [0,]

    try:
        x = it.next()           #The first call to a Python iterator must be .next().
        #^^Recall: x is a 2-tuple (CA,E,L), where CA is a colorability class assignment, E is the list of vertices whose lists aren't full, and L is the list of the sizes of the current lists of the vertices.

        print "(Start-up took",timestring(time.clock()-lap),".)"
        lap = time.clock()

        nodes_count += 1
        partial_count +=1
        node_code[0] += 1

        #At the top of this loop, we have obtained an x from the listAssignment iterator.  We check its representability -- unless it is representable, uncolorable, and full, we take the next x from the iterator.
        while True:

            if x[1]:
                partial_count += 1

            #If any full super-list assignment requires at least as many colors as vertices, we can turn back right now.
            if len(f) - len(x[0]) - max([f[v]-x[2][v] for v in range(len(f))]) < 1:
                b = []
                c = []
                first_count +=1

            else:

                b = assignmentRepresentationCheck(x[0])

                #If the list assignment is representable, we need to check with certainty whether it is colorable.
                if b:
                    c = assignmentCheckInduced(G,x[0])

                    #If the list assignment is not colorable and full, then add it to the list.
                    if not c and not x[1]:
                        z.append((x[0][:],nodes_count,node_code[:]))
                        print " "*10+"*"*10+" BAD LIST ASSIGNMENT:"+"*"*10
                        for i in range(len(x[0])):
                            print " "*10+subsetstring(x[0][i],G_copy.order()),"(",node_code[i],")"
                        print
                        bad_nums.append(nodes_count)
                        #Then switch b to be False/empty so that the iterator backtracks.
                        b = []


                #If the list assignment is not representable, we know it is colorable but don't want to take the time to find the coloring.
                else:
                    c = []
                    unrep_count +=1

            #Tell the listAssignment iterator to go to the next relevant list assignment (either moving forward or backtracking).
            send_value = (b,c)
            x = it.send(send_value)

            if len(x[0]) > len(node_code):
                node_code.append(1)
            else:
                node_code[len(x[0])-1] += 1
                for i in range(len(x[0]),len(node_code)):
                    del node_code[-1]

            nodes_count += 1
            if nodes_count%modulus == 0:
                print "-----> ",nodes_count,"partial list assignments inspected so far.  # bad LAs:",len(bad_nums),".  # full good LAs:",nodes_count-len(bad_nums)-partial_count,".  # caught by first check:",first_count,".  # caught by rep check:",unrep_count,".  # caught by direct coloring:",nodes_count-len(bad_nums)-first_count-unrep_count,"."
#                print "^^^^^>  #Bad list assignment node numbers:",bad_nums
                print "^^^^^>  Time taken on this batch of PLAs:",timestring(time.clock()-lap),"    Time taken on the whole (subgraph) job so far:",timestring(time.clock()-start)
                print "^^^^^>  Current partial list assignment:"
                for i in range(len(x[0])):
                    print subsetstring(x[0][i],G.order()),"(",node_code[i],")"
                print "\n"*3
                lap = time.clock()

    except StopIteration:           #After checking all of the list assignments, return all the bad ones.
        print "Finished.  %d PLAs visited which weren't pruned via CNS methods."%(nodes_count)
        return z,nodes_count





#----------------------------------------------------------------------------------
#Miscellaneous




#The tuple_restriction_list_maker fixes an ordering on all subsets (by least-ordered vertex) because it inherits fullCascade and then turns lists of colorability classes into tuples of tuples (to make them hashable).  However, our algorithm does not respect the vertex ordering but instead pokes at different vertex-rooted generators as it sees fit.  This procedure converts colClasses in the algorithm to the EXACT same tuple (i.e. appropriately orders it like tuple_restriction_list_maker would) so that it can be used as a key.
def colClasses_tupler(ordering,colClasses_copy):
    colClasses = colClasses_copy[:]
    colClasses_tuple = ()
    for vtx in ordering:
        vtx_tuple = ()
        li = []
        for i in range(len(colClasses)):
            if vtx in colClasses[i]:
                vtx_tuple += (set2tup(colClasses[i]),)
                li.append(i)
        if vtx_tuple:
            colClasses_tuple += vtx_tuple
        for i in reversed(li):
            del colClasses[i]
    return colClasses_tuple





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









#----------------------------------------------------------------------------------
#Not Minimal stuff










def fChoosableReturnBad(G_copy,f_copy,inprocess=3,print_mod=1000,maxIS=True,LAscheme="L",rev=False):
    
    if not LAscheme in ["L","E","LE","LD"]:
        print 'Error:  LAscheme must be one of "L", "E", "LE", or "LD".'
        return
    
    num_pruned = 0

    start = time.clock()

    if len(f_copy) != G_copy.order():
        print "Error: f-vector length does not match graph order."
        return

    if [x for x in f_copy if x<1]:
        print "Error: f-vector contains non-positive entries.  (If correct f, graph is not f-choosable.)"
        return

    G,f,same,low = reduce(G_copy,f_copy)
    
    if not same:
        print "The input configuration was reduced for having f-values either higher than the degree or equal to 1."
        if not f:
            print "In fact, it reduced entirely; the graph is f-choosable via a greedy argument."
            return []
        G.show()
        print "f:",f
    if low:
        print "Reduced f-vector contains non-positive entries.  Graph is not f-choosable."
        return False
    
    print "Testing:  Is the graph (%d vertices, %d edges) "%(G.order(),G.size())+str(f)+"-choosable? \n\nFirst, we consider whether the Combinatorial Nullstellensatz applies."
    boo = CNS_smart_print(G,f)
    print "    (Time so far:  "+timestring(time.clock()-start)+")"
    if boo:
        return []

    start2 = time.clock()

    print "\nNow we move on to our exhaustive search.  \nTo start up our process, we find and remove connected subgraphs which cannot be part of a bad list assignment, and then check the multiplicities of those that remain."

    it = restrictedListAssignments(G,f,inprocess_depth=inprocess,maxIS=maxIS,LAscheme=LAscheme,rev=rev)

    z = []

    nodes_count = 0
    full_count = 0
    bottom_count = 0
    bad_nums = []
    node_code = [0,]

    try:
        x = it.next()           #The first call to a Python iterator must be .next().
        #^^Recall: x is a tuple (CA,E,L,num_subsets,back,num_pruned), where CA is a colorability class assignment, E is the list of vertices whose lists aren't full, L is the list of the sizes of the current lists of the vertices, num_subsets is the list of the number of remaining subsets at each level of the subtree, and back is a boolean indicating whether we just got to this node after backtracking.  num_pruned is the sum of the number of subsets which were pruned inprocess at each step.  perminv is the reordering of vertices -- we use it only once to relabel the graph here.
        
        G,f = Gf_relabeling(G,f,x[6])
        
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

            if not x[1][-1]:
                full_count += 1

            #If there is an unfilled vertex which has no supergraph candidates remaining, then turn back!
            #(If the list assignment is full, then this cannot happen.)
            flag = True
            #print "E:",x[1][-1]
            #print "->",[sum([len(y) for y in x[3][-1].containment_dict[v]]) for v in x[1][-1]]
            #The top subgraph has already been removed from the containment_dict.  If it can be used again, then we need not check its vertices.  Otherwise, we check them.
            if x[7][set2tup(x[3][-1].subgraph)] < x[8][set2tup(x[3][-1].subgraph)]:
                for v in x[1][-1]:
                    if not v in x[3][-1].subgraph and sum([len(y) for y in x[3][-1].containment_dict[v]]) < 1:
                        flag = False
                        break
            else:
                for v in x[1][-1]:
                    if sum([len(y) for y in x[3][-1].containment_dict[v]]) < 1:
                        flag = False
                        break
                
            
            if flag:
                c = assignmentCheck(G,x[0])

                #If the list assignment is not colorable and full, then add it to the list.
                if not c and not x[1][-1]:
                    bad_nums.append(nodes_count)
                    z.append((x[0][:],nodes_count,node_code[:]))
                    print " "*10+"*"*10+" BAD LIST ASSIGNMENT #%d:"%(len(bad_nums))+"*"*10
                    for i in range(len(x[0])):
                        print " "*10+subsetstring(x[0][i],G.order()),"(",node_code[i],")"#|",x[3][i],")"
                    print
                    #Then switch c to be nonempty so that the iterator backtracks.
                    c = [0,]
            else:
                c = [0,]


            if x[4]:
                bottom_count += 1
            
            #Tell the listAssignment iterator to go to the next relevant list assignment (either moving forward or backtracking).
            x = it.send(c)

            if len(x[0]) > len(node_code):
                node_code.append(1)
            else:
                node_code[len(x[0])-1] += 1
                for i in range(len(x[0]),len(node_code)):
                    del node_code[-1]

            num_pruned += x[5]

            nodes_count += 1
            if nodes_count%print_mod == 0:
                print "-----> ",nodes_count,"partial list assignments inspected so far.  # bad LAs:",len(bad_nums),".  # full LAs:",full_count,".  # bottom LAs:",bottom_count,".  # inprocess prunes:",num_pruned,"."
                print "^^^^^>  Time taken on this batch of PLAs:",timestring(time.clock()-lap),"    Time taken on the whole (subgraph) job so far:",timestring(time.clock()-start)
                print "^^^^^>  Current partial list assignment:"
                print subsetstring(x[0][0],G.order()),"(",node_code[0],"|",len(x[3][0].generator_cascade[0]),")"
                for i in range(1,len(x[0])):
                    #print subsetstring(x[0][i],G.order()),"(",node_code[i],"|",x[3][i],")"
                    print subsetstring(x[0][i],G.order()),"(",node_code[i],"|",len([y for Y in x[3][i].generator_cascade for y in Y if y <= set(x[1][i])]),")"
                print "\n"*3
                lap = time.clock()
                #print "E_stack:",x[1]
                #for g in x[3]:
                    #print g.generator_cascade
                    #print

    except StopIteration:           #After checking all of the list assignments, return all the bad ones.
        if nodes_count < 1:#If CNS inconclusive but a vertex is not contained by any remaining subgraphs, then the normal start3 gets skipped.
            start3 = time.clock()
        end = time.clock()
        print "Finished.  %d PLAs visited which weren't pruned via CNS methods."%(nodes_count)
        print "# bad LAs:  %d.  # full LAs:  %d.  # bottom LAs:  %d.  # inprocess prunes:  %d"%(len(bad_nums),full_count,bottom_count,num_pruned)
        print "Time taken on initial CNS:       "+timestring(start2-start)
        print "Time taken on subgraph pruning:  "+timestring(start3-start2)
        print "Time taken on exhaustive search: "+timestring(end-start3)
        print "Total time taken on entire job:  "+timestring(end-start)
        return z









#Identical to fChoosableReturnBad, but does not store bad list assignments.  No bad_nums, z.
def fChoosable(G_copy,f_copy,inprocess=3,print_mod=1000,maxIS=True,LAscheme="L",rev=False):
    
    if not LAscheme in ["L","E","LE","LD"]:
        print 'Error:  LAscheme must be one of "L", "E", "LE", or "LD".'
        return
    
    num_pruned = 0

    start = time.clock()

    if len(f_copy) != G_copy.order():
        print "Error: f-vector length does not match graph order."
        return

    if [x for x in f_copy if x<1]:
        print "Error: f-vector contains non-positive entries.  (If correct f, graph is not f-choosable.)"
        return

    G,f,same,low = reduce(G_copy,f_copy)
    
    if not same:
        print "The input configuration was reduced for having f-values either higher than the degree or equal to 1."
        if not f:
            print "In fact, it reduced entirely; the graph is f-choosable via a greedy argument."
            return []
        G.show()
        print "f:",f
    if low:
        print "Reduced f-vector contains non-positive entries.  Graph is not f-choosable."
        return False
    
    print "Testing:  Is the graph (%d vertices, %d edges) "%(G.order(),G.size())+str(f)+"-choosable? \n\nFirst, we consider whether the Combinatorial Nullstellensatz applies."
    boo = CNS_smart_print(G,f)
    print "    (Time so far:  "+timestring(time.clock()-start)+")"
    if boo:
        return []

    start2 = time.clock()

    print "\nNow we move on to our exhaustive search.  \nTo start up our process, we find and remove connected subgraphs which cannot be part of a bad list assignment, and then check the multiplicities of those that remain."

    it = restrictedListAssignments(G,f,inprocess_depth=inprocess,maxIS=maxIS,LAscheme=LAscheme,rev=rev)

    nodes_count = 0
    full_count = 0
    bottom_count = 0
    node_code = [0,]

    try:
        x = it.next()           #The first call to a Python iterator must be .next().
        #^^Recall: x is a tuple (CA,E,L,num_subsets,back,num_pruned), where CA is a colorability class assignment, E is the list of vertices whose lists aren't full, L is the list of the sizes of the current lists of the vertices, num_subsets is the list of the number of remaining subsets at each level of the subtree, and back is a boolean indicating whether we just got to this node after backtracking.  num_pruned is the sum of the number of subsets which were pruned inprocess at each step.  perminv is the reordering of vertices -- we use it only once to relabel the graph here.
        
        G,f = Gf_relabeling(G,f,x[6])
        
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

            if not x[1][-1]:
                full_count += 1

            #If there is an unfilled vertex which has no supergraph candidates remaining, then turn back!
            #(If the list assignment is full, then this cannot happen.)
            flag = True
            #print "E:",x[1][-1]
            #print "->",[sum([len(y) for y in x[3][-1].containment_dict[v]]) for v in x[1][-1]]
            #The top subgraph has already been removed from the containment_dict.  If it can be used again, then we need not check its vertices.  Otherwise, we check them.
            if x[7][set2tup(x[3][-1].subgraph)] < x[8][set2tup(x[3][-1].subgraph)]:
                for v in x[1][-1]:
                    if not v in x[3][-1].subgraph and sum([len(y) for y in x[3][-1].containment_dict[v]]) < 1:
                        flag = False
                        break
            else:
                for v in x[1][-1]:
                    if sum([len(y) for y in x[3][-1].containment_dict[v]]) < 1:
                        flag = False
                        break
                
            
            if flag:
                c = assignmentCheck(G,x[0])

                #If the list assignment is not colorable and full, then report it and stop.
                if not c and not x[1][-1]:
                    print " "*10+"*"*10+" BAD LIST ASSIGNMENT:"+"*"*10
                    for i in range(len(x[0])):
                        print " "*10+subsetstring(x[0][i],G.order()),"(",node_code[i],")"#|",x[3][i],")"
                    print
                    #Then switch c to be nonempty so that the iterator backtracks.
                    break
            else:
                c = [0,]


            if x[4]:
                bottom_count += 1
            
            #Tell the listAssignment iterator to go to the next relevant list assignment (either moving forward or backtracking).
            x = it.send(c)

            if len(x[0]) > len(node_code):
                node_code.append(1)
            else:
                node_code[len(x[0])-1] += 1
                for i in range(len(x[0]),len(node_code)):
                    del node_code[-1]

            num_pruned += x[5]

            nodes_count += 1
            if nodes_count%print_mod == 0:
                print "-----> ",nodes_count,"partial list assignments inspected so far.  # full LAs:",full_count,".  # bottom LAs:",bottom_count,".  # inprocess prunes:",num_pruned,"."
                print "^^^^^>  Time taken on this batch of PLAs:",timestring(time.clock()-lap),"    Time taken on the whole (subgraph) job so far:",timestring(time.clock()-start)
                print "^^^^^>  Current partial list assignment:"
                print subsetstring(x[0][0],G.order()),"(",node_code[0],"|",len(x[3][0].generator_cascade[0]),")"
                for i in range(1,len(x[0])):
                    #print subsetstring(x[0][i],G.order()),"(",node_code[i],"|",x[3][i],")"
                    print subsetstring(x[0][i],G.order()),"(",node_code[i],"|",len([y for Y in x[3][i].generator_cascade for y in Y if y <= set(x[1][i])]),")"
                print "\n"*3
                lap = time.clock()
                #print "E_stack:",x[1]
                #for g in x[3]:
                    #print g.generator_cascade
                    #print

    except StopIteration:
        print "\nNo bad list assignments!\n"

    if nodes_count < 1:#If CNS inconclusive but a vertex is not contained by any remaining subgraphs, then the normal start3 gets skipped.
        start3 = time.clock()
    end = time.clock()
    print "Finished.  %d PLAs visited which weren't pruned via CNS methods."%(nodes_count)
    print "# full LAs:  %d.  # bottom LAs:  %d.  # inprocess prunes:  %d"%(full_count,bottom_count,num_pruned)
    print "Time taken on initial CNS:       "+timestring(start2-start)
    print "Time taken on subgraph pruning:  "+timestring(start3-start2)
    print "Time taken on exhaustive search: "+timestring(end-start3)
    print "Total time taken on entire job:  "+timestring(end-start)
    return









#Identical to fChoosable, except for any printing and returns (True/False,'one-word description') instead of a list.
def fChoosableNoPrint(G_copy,f_copy,inprocess=3,print_mod=1000,maxIS=True,LAscheme="L",rev=False):
    
    if not LAscheme in ["L","E","LE","LD"]:
        return (False,'error')
    
    if len(f_copy) != G_copy.order():
        return (False,'error')

    if [x for x in f_copy if x<1]:
        return (False,'error')

    G,f,same,low = reduce(G_copy,f_copy)
    
    if not same:
        if not f:
            return (True,'greedy')
    if low:
        return (False,'error')
    
    if CNS_smart(G,f):
        return (True,'CNS')


    it = restrictedListAssignmentsNoPrint(G,f,inprocess_depth=inprocess,maxIS=maxIS,LAscheme=LAscheme,rev=rev)

    try:
        x = it.next()           #The first call to a Python iterator must be .next().
        #^^Recall: x is a tuple (CA,E,L,num_subsets,back,num_pruned), where CA is a colorability class assignment, E is the list of vertices whose lists aren't full, L is the list of the sizes of the current lists of the vertices, num_subsets is the list of the number of remaining subsets at each level of the subtree, and back is a boolean indicating whether we just got to this node after backtracking.  num_pruned is the sum of the number of subsets which were pruned inprocess at each step.  perminv is the reordering of vertices -- we use it only once to relabel the graph here.
        
        G,f = Gf_relabeling(G,f,x[6])
        
        #At the top of this loop, we have obtained an x from the listAssignment iterator.  We check its representability -- unless it is representable, uncolorable, and full, we take the next x from the iterator.
        while True:

            #If there is an unfilled vertex which has no supergraph candidates remaining, then turn back!
            #(If the list assignment is full, then this cannot happen.)
            flag = True
            #The top subgraph has already been removed from the containment_dict.  If it can be used again, then we need not check its vertices.  Otherwise, we check them.
            if x[7][set2tup(x[3][-1].subgraph)] < x[8][set2tup(x[3][-1].subgraph)]:
                for v in x[1][-1]:
                    if not v in x[3][-1].subgraph and sum([len(y) for y in x[3][-1].containment_dict[v]]) < 1:
                        flag = False
                        break
            else:
                for v in x[1][-1]:
                    if sum([len(y) for y in x[3][-1].containment_dict[v]]) < 1:
                        flag = False
                        break
                
            
            if flag:
                c = assignmentCheck(G,x[0])

                #If the list assignment is not colorable and full, then report it and stop.
                if not c and not x[1][-1]:
                    return (False,'brute')
            else:
                c = [0,]

            #Tell the listAssignment iterator to go to the next relevant list assignment (either moving forward or backtracking).
            x = it.send(c)

    except StopIteration:
        #No bad list assignments!
        return (True,'brute')






def relabeling(G_copy,f_copy,cascade_copy,perminv):
    #Want perminv[i] to be relabeled as i.
    perm = [perminv.index(x) for x in range(len(perminv))]
    G = Graph(G_copy)
    G.relabel(perm)
    f = [f_copy[x] for x in perminv]
    
    cascade = []
    for x in cascade_copy:
        cascade.append([])
    for x in cascade_copy:
        for y in x:
            s = {perm[z] for z in y}
            cascade[min(s)].append(s)
    
    return G,f,cascade





def Gf_relabeling(G_copy,f_copy,perminv):
    #Want perminv[i] to be relabeled as i.
    perm = [perminv.index(x) for x in range(len(perminv))]
    G = Graph(G_copy)
    G.relabel(perm)
    f = [f_copy[x] for x in perminv]    
    return G,f




#inprocess_depth controls how many partial list assignments are checked for the reducibility when extending during the normal routine.
def restrictedListAssignments(G,f,inprocess_depth=3,maxIS=True,LAscheme="L",rev=False):
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

    initial_restriction_list = initial_restriction_list_maker(G, ordering, f, maxIS=maxIS)

    #num_subsets is going to keep track of how many nodes we have at each level of the subtree.  Entries are overwritten every time we go down to start the level.
    #num_subsets = [sum([len(x) for x in initial_restriction_list_dict[()]]),]
    #for i in range(sum(f)):
        #num_subsets.append(0)

    #print "    >>",num_subsets[0],"remaining subgraphs can potentially be used to build list assignments."
    print "    >>",sum([len(x) for x in initial_restriction_list]),"remaining subgraphs can potentially be used to build list assignments."
    
    if sum([len(x) for x in initial_restriction_list]) == 0:
        print "    >> So we can't build any bad list assignments!"
        return

    tenpercent = int(sum([len(x) for x in initial_restriction_list])/100)+1#int(num_subsets[0]/10)+1
    lap = time.clock()
    counter = 1
    begin = lap
    
    num_pruned = 0

    max_mult_init = {}
    #current_mult = {}
    for vtxList in initial_restriction_list:
        for A in vtxList:
            max_mult_init[set2tup(A)] = multiplicity_checker(G,ordering,f,A,minimal=False,maxIS=maxIS)
            #current_mult[set2tup(A)] = 0
            if counter%tenpercent==0:
                print "       ",counter,"multiplicities checked.  Time on batch:  "+timestring(time.clock()-lap)
                lap = time.clock()
            counter += 1


    print "        Time for all multiplicity checks:  "+timestring(time.clock()-begin)
    #print max_mult_init.values()
    print "    >> Multiplicities breakdown:"
    for i in range(1,max(max_mult_init.values())+1):
        print "        Number of subgraphs with multiplicity %d:  %d"%(i,len([x for x in max_mult_init.values() if x==i]))
        #if i==max(max_mult_init.values()):
            #print [x for x in max_mult_init.keys() if max_mult_init[x]==i]

    print "    >> Vertex containment breakdown for these leftover subgraphs:"
    for v in range(G.order()):
        print "        Number containing vertex %d: "%(v),sum([len([x for x in y if v in x]) for y in initial_restriction_list[:v+1]])

    #mult_sums = [sum([sum([max_mult_init[set2tup(x)] for x in y if v in x]) for y in initial_restriction_list_dict[()][:v+1]]) for v in range(G.order())]
    #print "    >> Vertex multiplicity-containment breakdown for these leftover subgraphs:"
    #for v in range(G.order()):
        #print "        Vertex %d: "%(v),mult_sums[v]

    #mi = min(mult_sums)
    #ma = max(mult_sums)+1
    #perminv = []
    #while mi < ma:
        #i = mult_sums.index(mi)
        #perminv.append(i)
        #mult_sums[i] = ma
        #mi = min(mult_sums)
    #print "Relabeling everything using perminv:",perminv

    #len_sums = [sum([sum([len(x)*max_mult_init[set2tup(x)] for x in y if v in x]) for y in initial_restriction_list_dict[()][:v+1]]) for v in range(G.order())]
    #print "    >> Vertex length-multiplicity-containment breakdown for these leftover subgraphs:"
    #for v in range(G.order()):
        #print "        Vertex %d: "%(v),len_sums[v]

    #mi = min(len_sums)
    #ma = max(len_sums)+1
    #perminv = []
    #while mi < ma:
        #i = len_sums.index(mi)
        #perminv.append(i)
        #len_sums[i] = ma
        #mi = min(len_sums)

    #print "Relabeling everything using perminv:",perminv
    #G,f,interim_cascade = relabeling(G,f,initial_restriction_list_dict[()],perminv)
    
    #interim_cascade = []
    subgraph_pool = {frozenset(x) for y in initial_restriction_list for x in y}
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
    
    print "Relabeling everything using perminv:",perminv
    G,f,interim_cascade = relabeling(G,f,initial_restriction_list,perminv)

    max_mult = {}
    current_mult = {}
    for vtxList in interim_cascade:
        for A in vtxList:
            max_mult[set2tup(A)] = max_mult_init[set2tup({perminv[a] for a in A})]
            current_mult[set2tup(A)] = 0

    ##No order
    #cascade = interim_cascade
    
    ##Multiplicities greatest to least, then lengths longest to shortest
    if rev:
        cascade = []
        for Y in interim_cascade:
            X = [(y,max_mult[set2tup(y)]) for y in Y]
            newX = []
            if X:
                mi=max([x[1] for x in X])
            while len(newX) < len(X):
                batch = [(i,len(X[i][0])) for i in range(len(X)) if X[i][1]==mi]
                if not batch:
                    mi -= 1
                    continue
                ma = max([x[1] for x in batch])
                #batch_sorted = []#sorts the indices so that sizes of subgraphs (with the same multiplicity) are decreasing
                while ma > 0:#<= max([x[1] for x in batch]):
                    mini_batch = [x[0] for x in batch if x[1]==ma]
                    #print mini_batch
                    #batch_sorted.extend(mini_batch)
                    ma -= 1
                    newX.extend([X[i][0] for i in mini_batch])
                mi -= 1
            cascade.append(newX)
    
    #Multiplicities least to greatest, then lengths shortest to longest
    else:
        cascade = []
        for Y in interim_cascade:
            X = [(y,max_mult[set2tup(y)]) for y in Y]
            newX = []
            mi=1
            while len(newX) < len(X):
                batch = [(i,len(X[i][0])) for i in range(len(X)) if X[i][1]==mi]
                if not batch:
                    mi += 1
                    continue
                ma = min([x[1] for x in batch])
                #batch_sorted = []#sorts the indices so that sizes of subgraphs (with the same multiplicity) are decreasing
                while ma <= max([x[1] for x in batch]):
                    mini_batch = [x[0] for x in batch if x[1]==ma]
                    #print mini_batch
                    #batch_sorted.extend(mini_batch)
                    ma += 1
                    newX.extend([X[i][0] for i in mini_batch])
                mi += 1
            cascade.append(newX)
    
    ##Multiplicities greatest to least, then lengths shortest to longest#just such a bad idea.
    #cascade = []
    #for Y in interim_cascade:
        #X = [(y,max_mult[set2tup(y)]) for y in Y]
        #newX = []
        #if not X:
            #cascade.append(newX)
            #continue
        #mi=max([x[1] for x in X])
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
            #mi -= 1
        #cascade.append(newX)
    
    ##Just lengths (currently shortest to longest)
    #cascade = []
    #for Y in interim_cascade:
        #X = list(Y)
        #newX = []
        #mi=1
        #while len(newX) < len(X):
            #batch = [x for x in X if len(x)==mi]
            #newX.extend(batch)
            #mi += 1
        #cascade.append(newX)

    #for Y in cascade:
        #for y in Y:
            #print max_mult[set2tup(y)], len(y), y
        #print
    #print cascade

    #First round, to make generator_stack nonempty.
    if LAscheme == "L":
        it = RestrictedConnectedSubgraphs_L(G,E_stack[0],cascade)#initial_restriction_list_dict[()])
    elif LAscheme == "E":
        it = RestrictedConnectedSubgraphs_E(G,E_stack[0],cascade)#initial_restriction_list_dict[()])
    elif LAscheme == "LE":
        it = RestrictedConnectedSubgraphs_LE(G,E_stack[0],cascade)#initial_restriction_list_dict[()])
    elif LAscheme == "LD":
        it = RestrictedConnectedSubgraphs_LD(G,E_stack[0],cascade)#initial_restriction_list_dict[()])
    generator_stack.append(it)                         #Stack of subgraph iterators; cth element generates the c-colorability class.

    if it.next_within(L,E_stack[0]):
        A = generator_stack[col].subgraph
        colClasses.append(A)
        current_mult[set2tup(A)] += 1
        for v in colClasses[col]:
            L[v]+=1
        E_stack.append([x for x in ordering if L[x]<f[x]])
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
        ticket = yield (colClasses,E_stack,L,generator_stack,back,num_pruned,perminv,current_mult,max_mult)
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

            #(Update the restrictionList from the current colClasses if we are in the early stages.)
            if len(colClasses) < inprocess_depth:
                num_pruned += generator_stack[col].update_restriction_list(f,colClasses,maxIS=maxIS)
            
            #Now go and get the next subgraph.
            if generator_stack[col].next_within(L,E_stack[-1]):
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
                #First, see if the generator can give us something new:
                if generator_stack[col].next_within(L,E_stack[-1]):
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
                if generator_stack[col].next_within_not_containing(L,E_stack[-1],D):
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
        while not generator_stack[col].next_within(L,E_stack[-1]):

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

        #Now we have successfully bumped the generator to another subgraph.  Push the subgraph onto the colClass stack.
        A = generator_stack[col].subgraph
        colClasses.append(A)
        current_mult[set2tup(A)] += 1
        for v in colClasses[col]:
            L[v]+=1
        E_stack.append([x for x in ordering if L[x]<f[x]])
        continue




#Identical to RestrictedListAssignments, but without print statements.
def restrictedListAssignmentsNoPrint(G,f,inprocess_depth=3,maxIS=True,LAscheme="L",rev=False):
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

    initial_restriction_list = initial_restriction_list_makerNoPrint(G, ordering, f, maxIS=maxIS)

    #num_subsets is going to keep track of how many nodes we have at each level of the subtree.  Entries are overwritten every time we go down to start the level.
    #num_subsets = [sum([len(x) for x in initial_restriction_list_dict[()]]),]
    #for i in range(sum(f)):
        #num_subsets.append(0)

    #print "    >>",num_subsets[0],"remaining subgraphs can potentially be used to build list assignments."
    #print "    >>",sum([len(x) for x in initial_restriction_list]),"remaining subgraphs can potentially be used to build list assignments."
    
    if sum([len(x) for x in initial_restriction_list]) == 0:
        #print "    >> So we can't build any bad list assignments!"
        return

    #tenpercent = int(sum([len(x) for x in initial_restriction_list])/100)+1#int(num_subsets[0]/10)+1
    #lap = time.clock()
    #counter = 1
    #begin = lap
    
    num_pruned = 0

    max_mult_init = {}
    #current_mult = {}
    for vtxList in initial_restriction_list:
        for A in vtxList:
            max_mult_init[set2tup(A)] = multiplicity_checker(G,ordering,f,A,minimal=False,maxIS=maxIS)
            #current_mult[set2tup(A)] = 0
            #if counter%tenpercent==0:
                #print "       ",counter,"multiplicities checked.  Time on batch:  "+timestring(time.clock()-lap)
                #lap = time.clock()
            #counter += 1

    #print "        Time for all multiplicity checks:  "+timestring(time.clock()-begin)
    #print max_mult_init.values()
    #print "    >> Multiplicities breakdown:"
    #for i in range(1,max(max_mult_init.values())+1):
        #print "        Number of subgraphs with multiplicity %d:  %d"%(i,len([x for x in max_mult_init.values() if x==i]))
        ##if i==max(max_mult_init.values()):
            ##print [x for x in max_mult_init.keys() if max_mult_init[x]==i]

    #print "    >> Vertex containment breakdown for these leftover subgraphs:"
    #for v in range(G.order()):
        #print "        Number containing vertex %d: "%(v),sum([len([x for x in y if v in x]) for y in initial_restriction_list[:v+1]])

    #mult_sums = [sum([sum([max_mult_init[set2tup(x)] for x in y if v in x]) for y in initial_restriction_list_dict[()][:v+1]]) for v in range(G.order())]
    #print "    >> Vertex multiplicity-containment breakdown for these leftover subgraphs:"
    #for v in range(G.order()):
        #print "        Vertex %d: "%(v),mult_sums[v]

    #mi = min(mult_sums)
    #ma = max(mult_sums)+1
    #perminv = []
    #while mi < ma:
        #i = mult_sums.index(mi)
        #perminv.append(i)
        #mult_sums[i] = ma
        #mi = min(mult_sums)
    #print "Relabeling everything using perminv:",perminv

    #len_sums = [sum([sum([len(x)*max_mult_init[set2tup(x)] for x in y if v in x]) for y in initial_restriction_list_dict[()][:v+1]]) for v in range(G.order())]
    #print "    >> Vertex length-multiplicity-containment breakdown for these leftover subgraphs:"
    #for v in range(G.order()):
        #print "        Vertex %d: "%(v),len_sums[v]

    #mi = min(len_sums)
    #ma = max(len_sums)+1
    #perminv = []
    #while mi < ma:
        #i = len_sums.index(mi)
        #perminv.append(i)
        #len_sums[i] = ma
        #mi = min(len_sums)

    #print "Relabeling everything using perminv:",perminv
    #G,f,interim_cascade = relabeling(G,f,initial_restriction_list_dict[()],perminv)
    
    #interim_cascade = []
    subgraph_pool = {frozenset(x) for y in initial_restriction_list for x in y}
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
    
    #print "Relabeling everything using perminv:",perminv
    G,f,interim_cascade = relabeling(G,f,initial_restriction_list,perminv)

    max_mult = {}
    current_mult = {}
    for vtxList in interim_cascade:
        for A in vtxList:
            max_mult[set2tup(A)] = max_mult_init[set2tup({perminv[a] for a in A})]
            current_mult[set2tup(A)] = 0

    ##No order
    #cascade = interim_cascade
    
    ##Multiplicities greatest to least, then lengths longest to shortest
    if rev:
        cascade = []
        for Y in interim_cascade:
            X = [(y,max_mult[set2tup(y)]) for y in Y]
            newX = []
            if X:
                mi=max([x[1] for x in X])
            while len(newX) < len(X):
                batch = [(i,len(X[i][0])) for i in range(len(X)) if X[i][1]==mi]
                if not batch:
                    mi -= 1
                    continue
                ma = max([x[1] for x in batch])
                #batch_sorted = []#sorts the indices so that sizes of subgraphs (with the same multiplicity) are decreasing
                while ma > 0:#<= max([x[1] for x in batch]):
                    mini_batch = [x[0] for x in batch if x[1]==ma]
                    #print mini_batch
                    #batch_sorted.extend(mini_batch)
                    ma -= 1
                    newX.extend([X[i][0] for i in mini_batch])
                mi -= 1
            cascade.append(newX)
    
    #Multiplicities least to greatest, then lengths shortest to longest
    else:
        cascade = []
        for Y in interim_cascade:
            X = [(y,max_mult[set2tup(y)]) for y in Y]
            newX = []
            mi=1
            while len(newX) < len(X):
                batch = [(i,len(X[i][0])) for i in range(len(X)) if X[i][1]==mi]
                if not batch:
                    mi += 1
                    continue
                ma = min([x[1] for x in batch])
                #batch_sorted = []#sorts the indices so that sizes of subgraphs (with the same multiplicity) are decreasing
                while ma <= max([x[1] for x in batch]):
                    mini_batch = [x[0] for x in batch if x[1]==ma]
                    #print mini_batch
                    #batch_sorted.extend(mini_batch)
                    ma += 1
                    newX.extend([X[i][0] for i in mini_batch])
                mi += 1
            cascade.append(newX)
    
    ##Multiplicities greatest to least, then lengths shortest to longest#just such a bad idea.
    #cascade = []
    #for Y in interim_cascade:
        #X = [(y,max_mult[set2tup(y)]) for y in Y]
        #newX = []
        #if not X:
            #cascade.append(newX)
            #continue
        #mi=max([x[1] for x in X])
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
            #mi -= 1
        #cascade.append(newX)
    
    ##Just lengths (currently shortest to longest)
    #cascade = []
    #for Y in interim_cascade:
        #X = list(Y)
        #newX = []
        #mi=1
        #while len(newX) < len(X):
            #batch = [x for x in X if len(x)==mi]
            #newX.extend(batch)
            #mi += 1
        #cascade.append(newX)

    #for Y in cascade:
        #for y in Y:
            #print max_mult[set2tup(y)], len(y), y
        #print
    #print cascade

    #First round, to make generator_stack nonempty.
    if LAscheme == "L":
        it = RestrictedConnectedSubgraphs_L(G,E_stack[0],cascade)#initial_restriction_list_dict[()])
    elif LAscheme == "E":
        it = RestrictedConnectedSubgraphs_E(G,E_stack[0],cascade)#initial_restriction_list_dict[()])
    elif LAscheme == "LE":
        it = RestrictedConnectedSubgraphs_LE(G,E_stack[0],cascade)#initial_restriction_list_dict[()])
    elif LAscheme == "LD":
        it = RestrictedConnectedSubgraphs_LD(G,E_stack[0],cascade)#initial_restriction_list_dict[()])
    generator_stack.append(it)                         #Stack of subgraph iterators; cth element generates the c-colorability class.

    if it.next_within(L,E_stack[0]):
        A = generator_stack[col].subgraph
        colClasses.append(A)
        current_mult[set2tup(A)] += 1
        for v in colClasses[col]:
            L[v]+=1
        E_stack.append([x for x in ordering if L[x]<f[x]])
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
        ticket = yield (colClasses,E_stack,L,generator_stack,back,num_pruned,perminv,current_mult,max_mult)
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

            #(Update the restrictionList from the current colClasses if we are in the early stages.)
            if len(colClasses) < inprocess_depth:
                num_pruned += generator_stack[col].update_restriction_list(f,colClasses,maxIS=maxIS)
            
            #Now go and get the next subgraph.
            if generator_stack[col].next_within(L,E_stack[-1]):
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
                #First, see if the generator can give us something new:
                if generator_stack[col].next_within(L,E_stack[-1]):
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
                if generator_stack[col].next_within_not_containing(L,E_stack[-1],D):
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
        while not generator_stack[col].next_within(L,E_stack[-1]):

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







