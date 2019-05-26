# To-Do:
    # After paper is finished, check references to Lemmas/Observations/Definitions.


#------------------------------------------------------
# Contents
#------------------------------------------------------

#    run_lemma_20
#    (main command)



#------------------------------------------------------
# Libraries, Comments, Global Variables
#------------------------------------------------------

# These scripts were written with the intent of using Python 2.7 before the switch to Python 3.  There should not be many differences to implement our code in Python 3 once this change is made.

from checking_realizations import *
# We use checking_realizations.py as a module which contains the common code for Lemmas ??, ??, and ??.



#------------------------------------------------------
# The Routines
#------------------------------------------------------



def run_lemma_20():
    # Verifies that each configuration in the target set is reducible by checking that all realizations are core-choosable.
    
    begin = time.clock()
    
    for rc in TargetSet:
        check_all_realizations_of_configuration(rc)
        print "\n"*20
    
    print "Finished with all configurations in the target set!"
    print "Total time for entire target set: "+timestring(time.clock()-begin)












print __name__

if __name__ == "__main__":
    run_lemma_20()




