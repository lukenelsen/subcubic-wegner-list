# To-Do:
    # After paper is finished, check references to Lemmas/Observations/Definitions.


#------------------------------------------------------
# Contents
#------------------------------------------------------

#    makeChainsDict
#    check_reducible_last_block
#    facial_length_list_generator
#    final_charge_lower_bound
#    run_lemma_28
#    (main command)



#------------------------------------------------------
# Libraries, Comments, etc.
#------------------------------------------------------

from fractions import Fraction
# In final_charge_lower_bound, we subtract noninteger rationals.

import time
# For displaying runtime in run_lemma_??.

# Instead of importing checking_realizations.py, we redefine timestring here.
def timestring(time):
    # Displays time length:  #m #.#s
    
    t = int(time)
    d = str(int((10*(time-t))%10))
    s = str(t%60)
    m = int(t/60)
    string = str(m)+"m "+s+"."+d+"s"
    return string



# For configuration names, we use notation slightly different from the notation used in the paper.  The configuration b:s_1...s_t is encoded as "cbas_1...s_t" where any "*" is replaced with "x".  (The "c" precedes the (c)entral face and the "a" precedes the facial length list (a)round the central face.)  For example, "9:3**4" is encoded as "c9a3xx4" and "4*555*4" is encoded as "cxa4x555x4".



#------------------------------------------------------
# The Routines
#------------------------------------------------------




def makeChainsDict():
    # This makes a dictionary (called ChainsDict) used for storing the reducible configurations as information which will be useful for our discharging.  For each central face length between 7 and 14, we want to store all the lists of entries which cannot be a sublist of a facial length list around a central face of that given length (in a minimal counterexample).  We call these lists "reducible chains" around the given central face lengths.
    # For example, since cxa355 is a reducible configuration, we know that ['3','5','5'] and ['5','5','3'] are reducible chains around any central face.  This means that neither of these two lists can be equal to three consecutive entries in a facial length list around any central face (including wrap-around).  So for each k between 7 and 14, ChainsDict[k] should contain ['3','5','5'] and ['5','5','3'].
    # For another example, since cxa375 is a reducible configuration, we know that ['3','x','5'] and ['5','x','3'] are reducible chains around a central 7-face.  (However, they are not reducible chains around a central 8-face.)  This means that certain facial length lists around a central 7-face cannot occur in a minimal counterexample, such as [3,4,5,5,5,5,5], [5,5,5,+,+,3,+], or [5,+,3,+,+,+,+].  So ChainsDict[7] should contain ['3','x','5'] and ['5','x','3'].
    # Since we are going to use reducible chains to compare against facial length lists with entries given as 3, 4, 5, or 6+, we are really only interested in chains with entries from {'3','4','5','x'}.  So the configuration cxa3x555 is of interest, but the configuration cxa356 is not.
    # Since we are going to use reducible chains to make observations about central faces of length 7-14, we are also only interested in chains that can possibly occur around such central faces.  So the configurations cxa375 and cxa3x4 are of interest, but the configuration c4a55 is not.
    
    # First, we define a list of the relevant configurations from the target set.  We comment out some of these configurations since they contribute nothing to ChainsDict.  (Although they could have been deleted or left uncommented, perhaps this will help for understanding which configurations are relevant to this stage.)  Next to each configuration we comment out, we also state a brief reason for its exclusion.
    # In addition to the above-stated reasons for excluding configurations, there is another:  some configurations have already been used in Lemma ?? to prevent the generation of any facial length lists containing them, and they serve no further purpose here.
    
    RelevantTargetSet =[
    # open chain:  start with 3 and end with 3
    #'c3a3',  # Lemma ?? only
    'cxa3x3',
    'cxa3xx3',
    'cxa3xxx3',
    'cxa3xx5x3',
    
    # open chain:  start with 3 and end with 4
    #'c3a4',  # Lemma ?? only
    'cxa3x4',
    'cxa35x4',
    'cxa3x54',
    'cxa3xx54',
    'cxa3xxx54',
    
    # open chain:  start with 3 and end with 5
    #'cxa355',  # Lemma ?? only
    'cxa375',
    'cxa3x555',
    
    # open chain:  start with 3 and end with 6
    #'c3a6',  # Contains a specified (6+)-face which would not be a central face of length 7-14
    #'cxa356',  # Contains a specified (6+)-face which would not be a central face of length 7-14
    
    # open chain:  start with 4 and end with 4
    #'c4a4',  # Lemma ?? only
    #'cxa454',  # Lemma ?? only
    #'cxa464',  # Contains a specified (6+)-face which would not be a central face of length 7-14
    'cxa474',
    'cxa4x54',
    'cxa4x55x4',
    'cxa4x555x4',
    'cxa4x4x4x4',
    
    # open chain:  start with 4 and end with 5
    'cxa4x45',
    #'cxa4555',  # Lemma ?? only
    'cxa4x5555',
    
    # open chain:  start with 4 and end with 6
    #'cxa456',  # Contains a specified (6+)-face which would not be a central face of length 7-14
    
    # open chain:  start with 5 and end with 5
    #'cxa535',  # Lemma ?? only
    #'cxa555555',  # Lemma ?? only
    
    # open chain:  start with 5 and end with 6
    #'cxa546',  # Contains a specified (6+)-face which would not be a central face of length 7-14
    
    # central 3-face
    'c3a57',
    
    # central 4-face
    #'c4a55',  # No way to have a central face of length 7-14
    #'c4a56',  # No way to have a central face of length 7-14
    'c4a57',
    #'c4a66',  # No way to have a central face of length 7-14
    'c4a585',
    'c4a676',
    'c4a686',
    
    # central 5-face
    #'c5a555',  # No way to have a central face of length 7-14
    #'c5a556',  # No way to have a central face of length 7-14
    #'c5a565',  # No way to have a central face of length 7-14
    #'c5a566',  # No way to have a central face of length 7-14
    'c5a575',
    #'c5a656',  # No way to have a central face of length 7-14
    #'c5a666',  # No way to have a central face of length 7-14
    #'c5a4x55',  # No way to have a central face of length 7-14
    #'c5a5x65',  # No way to have a central face of length 7-14
    #'c5a5585',  # No way to have a central face of length 7-14
    #'c5a5x66',  # No way to have a central face of length 7-14
    #'c5a55x6',  # No way to have a central face of length 7-14
    #'c5a56x6',  # No way to have a central face of length 7-14
    #'c5a6x66',  # No way to have a central face of length 7-14
    
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
    'c9a35x5'
        ]
    
    
    
    # Now we begin making ChainsDict!
    
    ChainsDict = {}  # Initialize ChainsDict.
    for cent in range(7,15):
        ChainsDict[cent] = []
    
    # Note that the keys for ChainsDict are int type.  However, everything we store will be strings.
    
    # For each reducible configuration given in the target set, we identify all the ways that a face of length 7-14 can be the central face for the configuration.  For each of these ways, we store the resulting chain of faces around the central faces (and the reverse of the chain, if it is distinct).
    for s in RelevantTargetSet:
        if s[1]=='x':
        
            # General case:  cxa s_1 s_2 ... s_t.  Store s_1 s_2 ... s_t (and the reverse) in each central face length's list.
            li = []  # li will be [s_1, s_2, ..., s_t] (where each s_i is a string).
            for a in s[3:]:
                li.append(a)
            for cent in range(7,15):
                ChainsDict[cent].append(li[:])  # Store s_1 s_2 ... s_t.
            rev_li = li[::-1]  # rev_li is the reverse of li.
            if rev_li != li:  # If the reverse is the same, then we don't want to store it again.
                for cent in range(7,15):
                    ChainsDict[cent].append(rev_li[:])  # Store s_t ... s_2 s_1 (reversed).
            
            # Additional special cases:  if t is small enough (at most 3), some of the outer faces could possibly be viewed as the central face.  (Section ?6??)
            # In the cases of t = 1,2, by inspection we have no configurations in the target set with an integer entry at least 7.
            # In the case of t = 3 for cxaABC, the configuration is equivalent to cBaAxC.  So if B > 6, then we should also store AxC and CxA in ChainsDict[B].
            if len(s)==6:  # If t==3.
                A,B,C = s[3],s[4],s[5]
                if B != 'x' and int(B) > 6:
                    ChainsDict[int(B)].append([A,'x',C])  # Store AxC in B's list of reducible chains.
                    if A != C:  # If the reverse is distinct.
                        ChainsDict[int(B)].append([C,'x',A])  # Store CxA in K's list of reducible chains.
        
        else:
            K = s[1]  # 'K' is easier for reference.
            
            # General case:  cKa s_1 s_2 ... s_t.  Store s_1 s_2 ... s_t (and the reverse) in K's list of reducible chains.
            if int(K) > 6:
                li = []  # li will be [s_1, s_2, ..., s_t] (where each s_i is a string).
                for a in s[3:]:
                    li.append(a)
                ChainsDict[int(K)].append(li[:])  # Store s_1 s_2 ... s_t.
                rev_li = li[::-1]  # rev_li is the reverse of li.
                if rev_li != li:  # If the reverse is the same, then we don't want to store it again.
                    ChainsDict[int(K)].append(rev_li[:])
            
            # Additional special cases:  if t is small enough (at most 3), some of the outer faces could possibly be viewed as the central face.  (Section ?6??)
            # In the case of t = 1, by inspection we have no configurations in the target set with an integer entry at least 7.
            # In the case of t = 2 for cKaAB, the configuration is equivalent to cAaKB.  So if A > 6, then we should also store KB and BK in A's list of reducible chains.   (Similarly store KA and AK in B's list if B > 6 and B != A.)
            if len(s)==5:  # If t==2.
                A,B = s[3],s[4]  # (Note that neither A nor B are 'x'.)
                if int(A) > 6 and A != K:  # If A == K, then switching won't give a distinct chain.
                    ChainsDict[int(A)].append([K,B])
                    if K != B:  # If K == B, then reversing the chain won't give a distinct chain.
                        ChainsDict[int(A)].append([B,K])
                if int(B) > 6 and B != K and B != A:  # If B is the same as K or A, then switching won't give a distinct chain.
                    ChainsDict[int(B)].append([K,A])
                    if K != A:  # If K == A, then reversing the chain won't give a distinct chain.
                        ChainsDict[int(B)].append([A,K])
            # In the case of t = 3 for cKaABC, the configuration is equivalent to cBaAKC.  So if B > 6 and B != K, then we should also store AKC and CKA in ChainsDict[B].  This is almost identical to the case when K = 'x' and t = 3.
            # Note:  For all the configurations in the target set with this form, B != 'x'.  Otherwise we would need to add to all central face lists instead of just one.
            if len(s)==6:  # If t==3.
                A,B,C = s[3],s[4],s[5]
                if int(B) > 6 and B != K:
                    ChainsDict[int(B)].append([A,K,C])  # Store AKC in B's list of reducible chains.
                    if A != C:  # If the reverse is distinct.
                        ChainsDict[int(B)].append([C,K,A])  # Store CKA in K's list of reducible chains.
    
    return ChainsDict
















def check_reducible_last_block(stack,block_length,central_face_length,ChainsDict):
    # For each reducible chain in ChainsDict[central_face_length], we check to see if it fits anywhere which overlaps the most recent block in the facial length list.  We do this by filling the rest of the given stack with ('x')s and then extending with another copy of the stack to create a longer list which can handle cyclic wrap-around.  Then, for each reducible chain, we compare every sublist of the chain's length which overlaps the new block to see if the sublist is a match.
    
    # We make a list of three copies of the current facial length list.  The second copy of the facial length list is where our inspections are centered.  The first and third copies are for possible cyclic wrap-around to the left and to the right.
    # The current facial length list is the stack followed by trailing ('x')s.
    wraparound_stack = ( stack + ['x']*(central_face_length-len(stack)) )*3
    
    # a and b are the indices of the first and last entries of the last block (in the second copy of the stack).
    a = central_face_length + len(stack) - (block_length + 1)
    b = central_face_length + len(stack) - 2
    # Since the last block of stack is followed by an 'x', we need to step back block_length + 1 steps from central_face_length + len(stack) to get to the first entry in the last block.  Likewise, we need to take 2 steps backs from central_face_length + len(stack) to get to the last entry.  See the example below.
        # stack = 53x5555x
        # block = 5555
        # block_length = 4
        # central_face_length = 10
        # len(stack) = 8
        #
        # 53x5555xxx53x5555xxx53x5555xxx  <-- wraparound_stack
        #              a  b z             <-- a = 13, b = 16, (z := central_face_length + len(stack) = 18)
        # 000000000011111111112222222222  <-- indexing guide:  tens place
        # 012345678901234567890123456789  <-- indexing guide:  ones place
    
    for reducible_chain in ChainsDict[central_face_length]:  # Inspect each reducible chain individually.
        chain_length = len(reducible_chain)
        for i in range(a-chain_length+1,b+1):  # i will be the starting index for the sublist.
            # Compare to the sublist starting at i.  Since the chain_length is at most 7 by inspection and the central_face_length is at least 7, the sublist won't overlap itself.
            sublist = wraparound_stack[i:i+chain_length]
            
            # We simply check the reducible_chain against the sublist entry by entry.  Since reversed directions of reducible chains are stored separately, we only need to check the reducible chain and the sublist in the same direction.
            match_bool = True
            for i in range(chain_length):
                
                # An 'x' in a reducible configuration matches everything.
                if reducible_chain[i]=='x':
                    continue  # These entries match, so check the next entries.
                
                # A specified number in a reducible configuration matches the same number in the sublist.
                elif reducible_chain[i]==str(sublist[i]):
                    continue  # These entries match, so check the next entries.
                
                # There is no other kind of match.  Different specified numbers do not match.
                # Note:  If we had added chains to ChainsDict which included entries not in {'3','4','5','x'}, there would still have been no other kind of matches.  A specified number in a reducible configuration does not match a '+' in a facial length list, because the '+' represents multiple numbers some of which are not equal to the specified number.  (For example, a 6 in the reducible configuration would not match the '+' in a facial length list because the '+' includes the possibility of a 10-face.)
                else:
                    break  # We found two entries which don't match, so this reducible chain does not match the sublist.
            
            else:
                return True  # If we find a match, then the new block has resulted in a reducible configuration.
            # If we did not match the sublist, then shift the sublist over and try again.
        # If that reducible chain does not show up anywhere overlapping our new block, then check the next reducible chain.
    return False  # If we don't find any match for any reducible chain, then the new block has not resulted in a reducible configuration.

















def facial_length_list_generator(central_face_length,ChainsDict):
    # Generates facial length lists around a specified face length which do not contain any configurations from the target set.
    # Facial length lists have entries from {3,4,5,'+'}, where '+' represents length 6+.
    # Facial length lists are generated via backtrack-and-search by building stacks of nonreducible blocks separated by a positive number of '+' entries.  This is more efficient than directly generating all the words from the alphabet {3,4,5,'+'} of a given length.  We also check at each forward step whether the new block forms any reducible configuration with the rest of the current stack.
    # We do not attempt to remove repeated words up to cyclic shifts or reversals.
    
    # First, we list the possible words from {3,4,5} which do not contain any reducible configurations.  According to Lemma ??, that list consists of the first 17 items listed below.  This can be verified by computer by generating {3,4,5}-words while checking for the occurence of any configuration in the target set.
    # In addition to those 17 blocks, we also include the singleton list ['+'].  We will use this final "block" as a spacer between the other blocks in our generation.
    # The comments at the end of each line give the index for each block.
    blocks = [
    [3],  # 0
    [3,5],  # 1
    [4],  # 2
    [4,5],  # 3
    [4,5,5],  # 4
    [5],  # 5
    [5,3],  # 6
    [5,4],  # 7
    [5,4,5],  # 8
    [5,4,5,5],  # 9
    [5,5],  # 10
    [5,5,4],  # 11
    [5,5,4,5],  # 12
    [5,5,4,5,5],  # 13
    [5,5,5],  # 14
    [5,5,5,5],  # 15
    [5,5,5,5,5],  #  16
    ['+']  # 17
        ]
    # The integer num is defined to be the index of the special ['+'] block, and we will use num in our generation below.
    num = 17
    
    # Our generation process will be a backtrack-and-search with two stacks:  block_stack and actual_stack.  block_stack will store integer values (the indices of the blocks given above) to represent which blocks are on the stack.  actual_stack will store values from {3,4,5,'+'} and will be extended by the blocks given above.  At any given time, we could convert information from one to the other.  However, it is easiest to update both simultaneously.
    block_stack = [0,num]  # Initialize block_stack for the loop.
    actual_stack = blocks[0] + blocks[num]  # Initialize actual_stack for the loop.
    
    # At the top of the loop, we have just put a new block onto our two stacks (if it was one of the 17 {3,4,5}-words, then it has also been automatically followed by the '+' spacer block).
    # The first thing we will do in the loop is check to see if this addition has violated the length requirement or has resulted in a reducible configuration.  If it has, we will backtrack and try the next block.  If it hasn't then we will try to add the least-indexed available block or, if the stack is actually full, then we will yield it and then backtrack.
    # The least-indexed available block will always be given by block_stack[0].  This somewhat (but not entirely) limits our duplication of facial length lists by removing some cyclic shifts.
    while block_stack[0] < num:  # While the first block is still a {3,4,5}-word (and not a '+').
        
        # Proceed with the small battery of tests to see if the current stack is viable.
        # First, check that actual_stack is not too long.  If so, just backtrack.
        # Second, check that actual_stack does not now result in a reducible configuration.  This is only necessary to check if the second-to-last block is a {3,4,5}-word and not the spacer '+'.  If it does result in a reducible configuration (in which case we have block_stack[-2] < num and check_reducible_last_block returns True), then backtrack.
        if ( len(actual_stack) <= central_face_length and 
            (block_stack[-2] == num or not check_reducible_last_block(actual_stack,len(blocks[block_stack[-2]]),central_face_length,ChainsDict) ) ) :
            
            # At this point, we have two options:  keep building or yield.
            if len(actual_stack) < central_face_length:  # If the stack is not full, then build it up more.
                # The first candidate for building is the block currently in the first position.
                block_stack.append(block_stack[0])
                actual_stack.extend(blocks[block_stack[0]])
                # Since the block currently in first position is a {3,4,5}-word (according to the condition on the while loop), we automatically add a spacer.
                block_stack.append(num)
                actual_stack.append('+')
                # Now instead of backtracking, return to the top of the loop.
                continue
            else:  # If the stack is full, yield it and then backtrack.
                yield actual_stack
        
        # Backtrack.  Go back to the last {3,4,5}-word and replace it with the next block.
        while block_stack[-1] == num:  # While the last block is the spacer '+'.
            actual_stack.pop()  # Pop the trailing ('+')s.
            block_stack.pop()  # Pop the trailing (num)s.
        actual_stack = actual_stack[:-len(blocks[block_stack[-1]])]  # Remove the last {3,4,5}-word.
        block_stack[-1] += 1  # Increment the last {3,4,5}=word to the next block.
        actual_stack.extend(blocks[block_stack[-1]])  # Add the actual entries of the new block.
        if block_stack[-1] < num:  # If the new block is a {3,4,5}-word (not '+'), then add '+'.
            block_stack.append(num)
            actual_stack.append('+')
    
    # By the time the first block gets to the spacer, there is only one facial length list left.  This is the list of all '+' symbols.  Since the facial length list contains all (6+)-faces, the central face will certainly have nonnegative charge and there is no need to yield it.
    return









def final_charge_lower_bound(central_face_length,facial_length_list):
    # Applies Lemma 25 for large faces; to obtain a lower bound for the final charge a central face has left after its adjacent small faces pull charge.
    # Assumes length at least 7 and full facial length list with entries from {3,4,5,'+'}.
    
    charge = Fraction(int(central_face_length-6),1)  # Set initial charge.
    
    # Apply Lemma 25.1:  A 3-face pulls 1 charge unless it is next to a 5-face -- then 3/2 charge.
    # threes is the list of indices with entries equal to 3.
    threes = [z for z in range(central_face_length) if facial_length_list[z]==3]
    for i in threes:
        # We check the entries immediately before or after the 3 at index i.  Python automatically interprets -1 as the last entry in the list, so we can simply use (i-1) for the index of the previous entry.  However, we mod out (i+1) by the list length just in case i is the last index for the list.
        if 5 in [facial_length_list[i-1],facial_length_list[(i+1)%central_face_length]]:  
            charge -= Fraction(3,2)
        else:
            charge -= Fraction(1,1)
    
    # Apply Lemma 25.2:  A 4-face pulls at most 2/3 charge if there is another 4-face two spots away or if the central face length is less than 9; otherwise,it pulls at most 1 charge.
    # fours is the list of indices with entries equal to 4.
    fours = [z for z in range(central_face_length) if facial_length_list[z]==4]
    for i in fours:
        # We check the entries two spots before or after the 4 at index i.  Again, we mod out by the list length for the cases that i is small or large.
        if (i-2)%central_face_length in fours or (i+2)%central_face_length in fours or central_face_length < 9:
            charge -= Fraction(2,3)
        else:
            charge -= Fraction(1,1)
    
    # Apply Lemma 25.3:  A 5-face pulls at most 1/4 charge if it has a neighboring 3-face, and it pulls at most 1/2 charge if it has two neighboring 5-faces and the central face length is at least 9; otherwise it pulls at most 1/3 charge.
    # fives is the list of indices with entries equal to 5.
    fives = [z for z in range(central_face_length) if facial_length_list[z]==5]
    for i in fives:
        # nbrs is the list containing the entry on the immediately left of i and the entry on the immediate right of i.  Again, we mod out for the (i+1) case.
        nbrs = [facial_length_list[i-1],facial_length_list[(i+1)%central_face_length]]
        if 3 in nbrs:  # If there is a 3 next to the 5
            charge -= Fraction(1,4)
        elif nbrs == [5,5] and central_face_length>=9:  # If the 5 is next to two 5s and the central face is big enough
            charge -= Fraction(1,2)
        else:
            charge -= Fraction(1,3)
    
    return charge









def run_lemma_28():
    # Generate all facial length lists (around central faces of length 7-14) which do not contain any configurations from the target set.  Check each of these to see if the lower bound on the final charge left with the central face is nonnegative.
    
    time_zero = time.clock()
    
    # Initialize list and dictionaries of reducible configurations.
    ChainsDict = makeChainsDict()
    
    unhappy_total = 0
    for k in range(7,15):
        lap = time.clock()
        print "Now analyzing central face length %d."%(k)
        print "Facial length list classes which pull too much charge:"
        
        # Generate the ways the eligible blocks form full stacks around this central face length.  (Without any reducible configurations.)
        gen_count = 0
        unhappy_count = 0
        for FLL in facial_length_list_generator(k,ChainsDict):
            gen_count += 1
            
            # Print and count any facial length list classes which pull too much charge.
            if final_charge_lower_bound(k,FLL) < 0:
                unhappy_count += 1
                s = ""
                for a in FLL:
                    s += str(a)
                print "        "+s+"  ( Lower Bound =",final_charge_lower_bound(k,FLL),")"
        
        unhappy_total += unhappy_count
        
        if unhappy_count==0:
            print "        (there are none)"
        print "    Of the %d facial length list classes around a %d-face with no reducible\n      configurations, %d have possibly negative charge.)"%(gen_count,k,unhappy_count)
        print "    Time to generate and check: "+timestring(time.clock()-lap)
        print
    
    print "In total, there are %d facial length list classes around central faces\n  of length between 7 and 14 which pull too much charge."%(unhappy_total)
    print "Total time: "+timestring(time.clock()-time_zero)





if __name__ == "__main__":
    run_lemma_28()




