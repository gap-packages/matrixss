###############################################################################
##
#W    code.gd     The Matrix Schreier Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
##    Dev start : 2004-01-10 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/lib/code_gd") := 
  "@(#)$Id$";


# Uses our version of Schreier-Sims to compute the order of a group
# Input:
#      G - a matrix group
DeclareGlobalFunction("MatrixGroupOrder");

# Compares results of the above function with the Order function for a bunch
# of classical matrix groups
# Input:
#      maxDegree - maximum matrix size to be tested
#      maxFieldSize - maximum finite field size to be tested
DeclareGlobalFunction("MatrixSchreierSimsTest");

# Check speed of routines against some classical groups and Ree groups
#
# Input:
#      maxDegree - maximum matrix size for classical groups to be tested
#      maxFieldSize - maximum finite field size to be tested
#      maxReeSize - maximum ReeGroup size
#      maxSuzukiSize - maximum SuzukiGroup size
DeclareGlobalFunction("MatrixSchreierSimsBenchmark");

# The Schreier-Sims algorithm for matrix groups
# Input:
#      G - a matrix group
# Output: list L with contents
#      L[1] - a base for G
#      L[2] - a strong generating set for G, corresponding to L[1]
#      L[3] - a list of orbits, ie SparseHashTables representing Schreier trees
#             the trees has the corresponding base points as roots
# Options:
#      SimpleSchreierTree - calculate coset representatives at the moment of
#                           creation of the Schreier trees, thus making them
#                           have height 1
#                           This should make the algorithm significantly faster
DeclareGlobalFunction("MatrixSchreierSims");


DeclareInfoClass("MatrixSchreierSimsInfo");

#E
