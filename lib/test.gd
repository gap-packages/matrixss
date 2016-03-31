###############################################################################
#1
#W    test.gd     The Matrix Schreier Sims package - Test & Benchmark routines
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-07-01
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
##
## These are the main routines for testing and benchmarking the package.
##
###############################################################################

Revision.("matrixss/lib/test_gd") := 
  "@(#)$Id$";

###############################################################################
##
#F MatrixSchreierSimsTest(maxDegree, maxFieldSize)
##
## Compares results of the above function with the built-in GAP Size method for
## a bunch of classical matrix groups. (`GL', `SL', etc)
## \beginitems
## `maxDegree' & maximum matrix size for classical matrix groups to be used
##             for testing
##
## `maxFieldSize' & maximum finite field size for classical matrix groups to be
##                used for testing
## \enditems
##
###############################################################################
DeclareGlobalFunction("MatrixSchreierSimsTest");

###############################################################################
##
#F MatrixSchreierSimsBenchmark(maxDegree, maxFieldSize, maxReeSize, maxSuzukiSize)
##
## Check speed of package routines against classical matrix groups and the
## matrix representations of Ree and Suzuki sporadic groups.
## \beginitems
## `maxDegree' & maximum matrix size for classical matrix groups to be used
##             for testing
##
## `maxFieldSize' & maximum finite field size for classical matrix groups to be
##                used for testing
##
## `maxReeSize' & maximum `ReeGroup' size, see "ref:Ree" in the reference 
##                manual.
##
## `maxSuzukiSize' & maximum `SuzukiGroup' size, see "ref:Sz" in the reference 
##                   manual.
## \enditems
##
###############################################################################
DeclareGlobalFunction("MatrixSchreierSimsBenchmark");

DeclareGlobalFunction("MatrixSchreierSimsSetBenchmark");

###############################################################################
#E
