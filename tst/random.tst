###############################################################################
##
#W    random.tst  The Matrix Schreier-Sims package                
#W                Probailistic algorithm test
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-01-24 
##
###############################################################################

gap> START_TEST("random.tst");;

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : Random, Verify, STCS);;
No order differences

#gap> MatrixSchreierSimsTest(2, 3, 2, 3 : Random, Verify, AlternatingActions);;
#No order differences

#gap> MatrixSchreierSimsTest(2, 3, 2, 3 : Random, Verify, STCS, AlternatingActions);;
#No order differences

gap> STOP_TEST("random.tst", 0);;

###############################################################################
#E
