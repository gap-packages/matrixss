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

gap> MatrixSchreierSimsBenchmark(5, 5, 0, 0 : SimpleSchreierTree);
Benchmark completed

gap> MatrixSchreierSimsBenchmark(5, 5, 0, 0 : SimpleSchreierTree, UseRandomSS);
Benchmark completed

gap> MATRIXSS_RandomSchreierSimsBenchmark(5, 5, 0, 0);
Benchmark completed

gap> STOP_TEST("bench.tst", 10000);
