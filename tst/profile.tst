###############################################################################
##
#W    profile.tst The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-03-22 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

gap> START_TEST("$Id$");

gap> MatrixSchreierSimsBenchmark(4, 4, 0, 0 : SimpleSchreierTree);
Benchmark completed

gap> STOP_TEST("profile.tst", 10000);
