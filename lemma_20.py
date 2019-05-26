# To-Do:
    # After paper is finished, check references to Lemmas/Observations/Definitions.


#------------------------------------------------------
# Contents
#------------------------------------------------------

#    check_all_realizations_of_configuration
#    run_lemma_20
#    (main command)



#------------------------------------------------------
# Libraries, Comments, Global Variables
#------------------------------------------------------

# These scripts were written with the intent of using Python 2.7 before the switch to Python 3.  There should not be many differences to implement our code in Python 3 once this change is made.

from checking_realizations import *
# We use checking_realizations.py as a module which contains the common code for Lemmas ??, ??, and ??.



# For configuration names, we use notation slightly different from the notation used in the paper.  The configuration b:s_1...s_t is encoded as "cbas_1...s_t" where any "*" is replaced with "x".  (The "c" precedes the (c)entral face and the "a" precedes the facial length list (a)round the central face.)  For example, "9:3**4" is encoded as "c9a3xx4" and "4*555*4" is encoded as "cxa4x555x4".



#------------------------------------------------------
# The Routines
#------------------------------------------------------



def check_all_realizations_of_configuration(config_str):
    # Given a configuration as string in the c/a notation described in a comment above, generates and checks every realization for core-choosability.
    # This simply prepares the natural core subgraph and feeds it into check_all_realizations_from_initial_plane_graph.
    
    print "Configuration:",config_str
    print "Checking all realizations for core-choosability.\n"
    
    # The NaturalCoreSubgraph class prepares all the information we need about the natural core subgraph of the configuration.
    NCS = NaturalCoreSubgraph(config_str)
    
    check_all_realizations_from_initial_plane_graph(NCS,forbid_identifications={})










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




