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
