###############################################################################
##
#W    test.gd     The Matrix Schreier Sims package - Test code                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
##    Dev start : 2004-07-01
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/lib/test_gd") := 
  "@(#)$Id$";

# Uses our version of Schreier-Sims to compute the order of a group
# Input:
#      G - a matrix group
DeclareGlobalFunction("MatrixGroupOrder");

# Uses our version of random Schreier-Sims to compute the order of a group
# Input:
#      G - a matrix group
DeclareGlobalFunction("RandomMatrixGroupOrder");

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

#E
