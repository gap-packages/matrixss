###############################################################################
##
#W    stcs.gd  The Matrix Schreier Sims package 
#W             Schreier-Todd-Coxeter-Sims implementation
##
#H    File      : $RCSfile$
#H    Author    : Henrik B��rnhielm
##    Dev start : 2004-07-01 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/lib/stcs_gd") := 
  "@(#)$Id$";

# The Schreier-Todd-Coxeter-Sims algorithm for matrix groups
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
#      ExtendSchreierTree - Do not recompute Schreier trees at each run of a
#                           given level, but extend the Schreier trees from the
#                           last run at that level.
#      AlternatingActions - Always prepend a base point with the line that
#                           contains it, using the projective action on the 
#                           line.
#      CleverBasePoints   - Choose an initial list of base points using 
#                           BasisVectorsForMatrixAction, which is due to
#                           O'Brien and Murray.
DeclareGlobalFunction("MatrixSchreierToddCoxeterSims");

#E