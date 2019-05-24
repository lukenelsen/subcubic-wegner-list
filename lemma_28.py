# To-Do:
    # Revise/finalize comments in makeChainsDict.
    # After paper is finished, rename lemma_28 appropriately.


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

import time
# For displaying runtime in run_lemma_28.

def timestring(time):
    # Displays time length:  ?m ?.?s
    
    t = int(time)
    d = str(int((10*(time-t))%10))
    s = str(t%60)
    m = int(t/60)
    string = str(m)+"m "+s+"."+d+"s"
    return string



#------------------------------------------------------
# The Routines
#------------------------------------------------------

# The comments in makeChainsDict can certainly use some fine-tuning.


def makeChainsDict():
    # This makes a dictionary (called ChainsDict) used for storing the reducible configurations as information which will be useful for our discharging.  For each central face length between 7 and 14, we want to store all the lists of entries which cannot be a sublist of a facial length list around a central face of that given length (in a minimal counterexample).  We call these lists "reducible chains".
    # For example, since cxa375 is a reducible configuration, we know that ['3','x','5'] and ['5','x','3'] are reducible chains around a central 7-face.  However, they are not reducible chains around a central 8-face.
    # Since we are going to use reducible chains to compare against facial length lists with entries given as 3, 4, 5, or 6+, we are really only interested in chains with entries from {'3','4','5','x'}.

    # First, define the reducible configurations in the target set.  We use notation which carried over from initial work rather than the final notation used in the paper.  Here, we use 'x' instead of '*' to indicate the wild card symbol.  A 'cxa' prefix means an unspecified central face length.  (e.g., 'cxa3x54' is the configuration 3*54.)  A 'cKa' prefix for an integer K means a central face of specified length K.  The 'c' stands for "central" and the 'a' stands for "around".
    # We could leave some of these configurations out of the list since they are irrelevant to ChainsDict, but we leave them all for a sense of completeness.  (e.g., the configuration c4a55 has no way of having a central face of length between 7 and 14.)  (e.g., the configuration cxa454 can have a central face of length 7 through 4, but already serves its purpose in Lemma ?? by controlling the {3,4,5}-words.......)
    TargetSet =[
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
    'c9a35x5'
    'c9a455x4',
    'c9a455x5',
    'c9a545x5',
    'c9a35x5'
        ]

    ChainsDict = {}  # Initialize ChainsDict.
    for cent in range(7,15):
        ChainsDict[cent] = []

    # For each reducible configuration given in the target set, we identify all the ways that a face of length 7-14 can be the central face for the configuration.  For each of these ways, we store the resulting chain of faces around the central faces (and the reverse of the chain, if it is distinct).
    for s in TargetSet:
        if s[1]=='x':

            # General case:  cxa s_1 s_2 ... s_t.  Store s_1 s_2 ... s_t (and possibly the reverse) in each central face length's list.
            li = []
            for a in s[3:]:
                li.append(a)
            for cent in range(7,15):
                ChainsDict[cent].append(li[:])  # store s_1 s_2 ... s_t
            rev_li = li[::-1]
            if rev_li != li:
                for cent in range(7,15):
                    ChainsDict[cent].append(rev_li[:])  # store s_t ... s_2 s_1

            # Additional special cases:  if t is small enough, some of the outer faces could possibly be viewed as the central face.  (Observation ???)
            # In the case of t = 2, we have no configurations in the target set with an integer entry at least 7.
            # In the case of t = 3 for cxaJKL if K > 6, then we should also store cKaJxL/cKaLxJ.
            if len(s)==6 and s[4] != 'x' and int(s[4]) > 6:  # If K != 'x' and K > 6
                ChainsDict[int(s[4])].append([s[3],'x',s[5]])  # store JxL in K's list
                if s[3] != s[5]:#If J != L
                    ChainsDict[int(s[4])].append([s[5],'x',s[3]])  # store LxJ in K's list

        else:

            # General case:  cKa s_1 s_2 ... s_t.  Store s_1 s_2 ... s_t (and possibly the reverse) in K's list.
            if int(s[1]) > 6:
                li = []
                for a in s[3:]:
                    li.append(a)
                ChainsDict[int(s[1])].append(li[:])
                rev_li = li[::-1]
                if rev_li != li:
                    ChainsDict[int(s[1])].append(rev_li[:])

            # Additional special cases:  if t is small enough, some of the outer faces could possibly be viewed as the central face.  (Observation ???)
            # In the case of t = 1, we have no configurations in the target set with an integer entry at least 7.
            # In the case of t = 2 for cKaJL if J > 6, then we should also store cJaKL/cJaLK.  (Similarly if L > 6.)
            if len(s)==5:#(Neither J nor K are 'x'.)
                if int(s[3]) > 6 and s[3] != s[1]:  # If J > 6 and J != K
                    ChainsDict[int(s[4])].append([s[1],s[4]])  # store KL in J's list.
                    if s[1] != s[4]:  # If K != L
                        ChainsDict[int(s[3])].append([s[4],s[1]])  # store LK in J's list.
                if int(s[4]) > 6 and s[4] != s[1] and s[4] != s[3]:  # If L > 6 and L != K and L != J
                    ChainsDict[int(s[4])].append([s[1],s[3]])  # store KJ in L's list.
                    if s[1] != s[3]:  # If K != J
                        ChainsDict[int(s[4])].append([s[3],s[1]])  # store JK in L's list.
            # In the case of t = 3 for cKaJLM if L > 6, then we should also store cLaJKM/cLaMKJ.  This is almost identical to the case when s[1] = 'x'.
            # For all the configurations in the target set with this form, K != 'x'.
            if len(s)==6 and int(s[4]) > 6 and s[4] != s[1]:  # If L > 6 and L != K
                ChainsDict[int(s[4])].append([s[3],s[1],s[5]])  # store JKM in L's list
                if s[3] != s[5]:#If M != J
                    ChainsDict[int(s[4])].append([s[5],s[1],s[3]])  # store MKJ in L's list

    return ChainsDict





























def check_reducible_last_block(stack,block_length,central_face_length,ChainsDict):
    # For each reducible chain in ChainsDict[central_face_length], we check to see if it fits anywhere which overlaps the most recent block in the facial length list.  We do this by filling the rest of the given stack with ('x')s and then extending with another copy of the stack to create a longer list which can handle cyclic wrap-around.  Then, for each reducible chain, we compare every sublist of the chain's length which overlaps the new block to see if the sublist is a match.

    # We make a list of three copies of the current facial length list.  The second copy of the facial length list is where our inspections are centered.  The first and third copies are for possible cyclic wrap-around to the left and to the right.
    # The current facial length list is the stack followed by trailing ('x')s.
    wraparound_stack = ( stack + ['x']*(central_face_length-len(stack)) )*3
    
    # a and b are the indices of the first and last entries of the last block (in the second copy of the stack).
    a = central_face_length + len(stack) - 1 - block_length  # add comment
    b = central_face_length + len(stack) - 2  # add comment
    
    for reducible_chain in ChainsDict[central_face_length]:  # Inspect each reducible chain individually.
        chain_length = len(reducible_chain)
        for i in range(a-chain_length+1,b+1):  # i will be the starting index for the sublist.
            # Compare to the sublist starting at i.  Since the chain_length is at most 7 and the central_face_length is at least 7, the sublist won't overlap itself.
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
                
                # There is no other kind of match.  Different numbers do not match.  Also, a specified number in a reducible configuration does not match a '+' in a facial length list, because the '+' represents multiple numbers some of which are not equal to the specified number.  (For example, a 6 in the reducible configuration does not match the '+' in a facial length list because the '+' might be a 10-face.)
                else:
                    match_bool = False
                    break  # We found two entries which don't match, so this reducible chain does not match the sublist.
            
            if match_bool:
                return True  # If we find a match, then the new block has resulted in a reducible configuration.
            # If we did not match the sublist, then shift the sublist over and try again.
        # If that reducible chain does not show up anywhere overlapping our new block, then check the next reducible chain.
    return False  # If we don't find any match for any reducible chain, then the new block has not resulted in a reducible configuration.

















def facial_length_list_generator(central_face_length,ChainsDict):
    #Generates facial length lists around a specified face length which do not contain any configurations from the target set.
    #Facial length lists have entries from {3,4,5,'+'}, where '+' represents length 6+.
    #Facial length lists are generated via backtrack-and-search by building stacks of nonreducible blocks separated by a positive number of '+' entries.  This is more efficient than directly generating all the words from the alphabet {3,4,5,'+'} of a given length.  We also check at each forward step whether the new block forms any reducible configuration with the rest of the current stack.
    #We do not attempt to remove repeated words up to cyclic shifts or reversals.

    # First, we list the possible words from {3,4,5} which do not contain any reducible configurations.  According to Lemma 20-something?, that list consists of the first 17 items listed below.  This can be verified by computer by generating {3,4,5}-words while checking for the occurence of any configuration in the target set.
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
    
    # By the time the first block gets to the spacer, there is only one facial length list left.  This is the list of all '+' symbols.  Since the facial length list contains all (6+)-faces, the central face will certainly have nonnegative charge and there is no need to check it.
    return









def final_charge_lower_bound(central_face_length,facial_length_list):
    # Applies Lemma 25 for large faces; to obtain a lower bound for the final charge a central face has left after its adjacent small faces pull charge.
    # Assumes length at least 7 and full facial length list with entries from {3,4,5,'+'}.
    
    charge = central_face_length-6  # Set initial charge.
    
    # Apply Lemma 25.1:  A 3-face pulls 1 charge unless it is next to a 5-face -- then 3/2 charge.
    # threes is the list of indices with entries equal to 3.
    threes = [z for z in range(central_face_length) if facial_length_list[z]==3]
    for i in threes:
        # We check the entries immediately before or after the 3 at index i.  Python automatically interprets -1 as the last entry in the list, so we can simply use (i-1) for the index of the previous entry.  However, we mod out (i+1) by the list length just in case i is the last index for the list.
        if 5 in [facial_length_list[i-1],facial_length_list[(i+1)%central_face_length]]:  
            charge -= 3/2
        else:
            charge -= 1

    # Apply Lemma 25.2:  A 4-face pulls at most 2/3 charge if there is another 4-face two spots away or if the central face length is less than 9; otherwise,it pulls at most 1 charge.
    # fours is the list of indices with entries equal to 4.
    fours = [z for z in range(central_face_length) if facial_length_list[z]==4]
    for i in fours:
        # We check the entries two spots before or after the 4 at index i.  Again, we mod out by the list length for the cases that i is small or large.
        if (i-2)%central_face_length in fours or (i+2)%central_face_length in fours or central_face_length < 9:
            charge -= 2/3
        else:
            charge -= 1

    # Apply Lemma 25.3:  A 5-face pulls at most 1/4 charge if it has a neighboring 3-face, and it pulls at most 1/2 charge if it has two neighboring 5-faces and the central face length is at least 9; otherwise it pulls at most 1/3 charge.
    # fives is the list of indices with entries equal to 5.
    fives = [z for z in range(central_face_length) if facial_length_list[z]==5]
    for i in fives:
        # nbrs is the list containing the entry on the immediately left of i and the entry on the immediate right of i.  Again, we mod out for the (i+1) case.
        nbrs = [facial_length_list[i-1],facial_length_list[(i+1)%central_face_length]]
        if 3 in nbrs:  # If there is a 3 next to the 5
            charge -= 1/4
        elif nbrs == [5,5] and central_face_length>=9:  # If the 5 is next to two 5s and the central face is big enough
            charge -= 1/2
        else:
            charge -= 1/3

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
        print "Facial length lists which pull too much charge:"

        # How can these blocks form full stacks around this central face length?  (Without any reducible configurations.)
        gen_count = 0
        unhappy_count = 0
        for FLL in facial_length_list_generator(k,ChainsDict):
            gen_count += 1

            # Which ones pull too much charge?
            if final_charge_lower_bound(k,FLL) < 0:
                unhappy_count += 1
                s = ""
                for a in FLL:
                    s += str(a)
                print "        "+s+"  ( Lower Bound =",final_charge_lower_bound(k,FLL),")"

        unhappy_total += unhappy_count

        print "    (Of the %d facial length lists, %d have possibly negative charge.)"%(gen_count,unhappy_count)
        print "    Time to generate and check: "+timestring(time.clock()-lap)
        print

    print "In total, there are %d facial length lists which pull too much charge."%(unhappy_total)
    print "Total time: "+timestring(time.clock()-time_zero)





if __name__ == "__main__":
    run_lemma_28()




