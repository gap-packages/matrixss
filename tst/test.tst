###############################################################################
##
#W    test.tst The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-01-24 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

gap> START_TEST("$Id$");

gap> MatrixSchreierSimsTest(5, 5 : Random, Verify, STCS, AlternatingActions);
No order differences

gap> STOP_TEST("test.tst", 10000);
