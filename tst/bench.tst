###############################################################################
##
#W    bench.tst The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B��rnhielm
#H    Dev start : 2004-03-22 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

gap> START_TEST("$Id$");

gap> MatrixSchreierSimsBenchmark(6, 6, 27, 32 : AlternatingActions, SimpleSchreierTree, CleverBasePoints);
Benchmark completed

gap> MATRIXSS_RandomSchreierSimsBenchmark(6, 6, 27, 32);
Benchmark completed

gap> STOP_TEST("bench.tst", 10000);
