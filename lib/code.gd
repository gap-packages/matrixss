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
# of general linear groups
# Input:
#      maxDegree - maximum GL degree to be tested
#      maxFieldSize - maximum GL finite field size to be tested
DeclareGlobalFunction("MatrixSchreierSimsTest");

# The Schreier-Sims algorithm for matrix groups
# Input:
#      G - a matrix group
# Output: list L with contents
#      L[1] - a base for G
#      L[2] - a strong generating set for G, corresponding to L[1]
#      L[3] - a list of Schreier trees
DeclareGlobalFunction("MatrixSchreierSims");


DeclareInfoClass("MatrixSchreierSimsInfo");

#E
