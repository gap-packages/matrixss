###############################################################################
##
#W    random.tst  The Matrix Schreier-Sims package                
#W                Probailistic algorithm test
##
#H    File      : $RCSfile$
#H    Author    : Henrik B��rnhielm
#H    Dev start : 2004-01-24 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

gap> START_TEST("$Id$");;

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : Random, Verify, STCS);;
No order differences

#gap> MatrixSchreierSimsTest(2, 3, 2, 3 : Random, Verify, AlternatingActions);;
#No order differences

#gap> MatrixSchreierSimsTest(2, 3, 2, 3 : Random, Verify, STCS, AlternatingActions);;
#No order differences

gap> STOP_TEST("random.tst", 0);;
$Id$
GAP4stones: 0

###############################################################################
#E
