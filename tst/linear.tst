###############################################################################
##
#W    linear.tst  The Matrix Schreier-Sims package                
#W                Nearly linear time algorithm test
##
#H    Author    : Henrik B��rnhielm
#H    Dev start : 2004-01-24 
##
###############################################################################

gap> START_TEST("linear.tst");

gap> MatrixSchreierSimsTest(3, 3 : Linear);;
No order differences

gap> STOP_TEST("linear.tst", 0);

###############################################################################
#E
