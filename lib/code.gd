###############################################################################
#1
#W    code.gd     The Matrix Schreier Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
#H    Dev start : 2004-01-10 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
##
## These are the general declarations used by the package. Most notably the
## attribute `StabChainMatrixGroup' which is the core of the package 
## functionality.
##
###############################################################################

Revision.("matrixss/lib/code_gd") := 
  "@(#)$Id$";

###############################################################################
##
#A StabChainMatrixGroup(G)
##
## Declare new attribute for storing base, SGS and Schreier trees.
## The attribute is computed using the Schreier-Sims algorithm for finite
## matrix groups, which is the main content of the package.
##
## The attribute is a record with two components:
## \beginitems
## `SchreierStructure' & the main information structure, see "ssInfo".
##
## `SGS' & a list of the strong generators
## \enditems
##
## The corresponding attribute operations are aware of a few Options.
## \beginitems
## `SimpleSchreierTree' & calculate coset representatives at the moment of
##                      creation of the Schreier trees, thus making them
##                      have height 1.
##                      This should make the algorithm significantly 
##                      faster.
##
## `ExtendSchreierTree' & Do not recompute Schreier trees at each run of a
##                        given level, but extend the Schreier trees from 
##                        the last run at that level.
##
## `AlternatingActions' & Always prepend a base point with the line that
##                        contains it, using the projective action on the 
##                        line.
##
## `CleverBasePoints' & Choose an initial list of base points using 
##                      `BasisVectorsForMatrixAction', which is made by
##                      O'Brien and Murray.
##
## `ShallowSchreierTree' & Create the Schreier trees using the method described
##                         by Babai et al (1991) which guarantees logarithmic
##                         depth. This Option is only used when the option
##                         `SimpleSchreierTree' is <not> defined.
## 
## `Random' & Use probabilistic algorithm.
##
## `Linear' & Use nearly linear-time algorithm. This takes precedence over
##            `Random', if it is present.
##
## \enditems
##
###############################################################################
DeclareAttribute("StabChainMatrixGroup", IsMatrixGroup and IsFinite);

###############################################################################
##
#V MatrixSchreierSimsInfo
##
## The {\GAP} InfoClass used by the package, for debugging purposes.
##
###############################################################################
DeclareInfoClass("MatrixSchreierSimsInfo");

###############################################################################
##
#V MATRIXSS_DEBUGLEVEL
##
## The internal debugging level. This is really obsolete and the above info
## class should be used instead.
##
###############################################################################
MATRIXSS_DEBUGLEVEL := 0;

###############################################################################
##
#V MATRIXSS_BasePointStore
##
## A list of hopefully good base points, ie base points with small orbits.
## They are fetched with `BasisVectorsForMatrixAction' which is due to
## O{\rq}Brien and Murray, and as long as the list is non-empty, 
## new base points will be shifted from it.
##
###############################################################################
MATRIXSS_BasePointStore := [];

###############################################################################
##
#V MATRIXSS_SubProdGroups
##
## A Dictionary of SymmetricGroups, used when permuting random subproducts.
##
###############################################################################
MATRIXSS_SubProdGroups := NewDictionary(10, true, Integers);

###############################################################################
##
#F MatrixGroupOrderStabChain(ssInfo)
##
## Computes the order of the group defined by the given Schreier trees, see
## "ssInfo".
##
###############################################################################
DeclareGlobalFunction("MatrixGroupOrderStabChain");

###############################################################################
##
#F ElementToSLP(G, g)
##
## Computes an SLP of g in the generators of G, using StabChainMatrixGroup.
##
###############################################################################
DeclareGlobalFunction("ElementToSLP");

#E
