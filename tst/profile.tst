###############################################################################
##
#W    profile.tst The Matrix Schreier-Sims package                
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-03-22 
##
###############################################################################

gap> START_TEST("profile.tst");

gap> MatrixSchreierSimsBenchmark(5, 5, 0, 0 : Random);
Benchmark completed

gap> STOP_TEST("profile.tst", 10000);
