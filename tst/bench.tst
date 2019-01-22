###############################################################################
##
#W    bench.tst The Matrix Schreier-Sims package                
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-03-22 
##
###############################################################################

gap> START_TEST("bench.tst");

gap> MatrixSchreierSimsBenchmark(4, 4, 0, 32 : AlternatingActions, Random);
Benchmark completed

gap> MATRIXSS_RandomSchreierSimsBenchmark(4, 4, 0, 32);
Benchmark completed

gap> STOP_TEST("bench.tst", 10000);
