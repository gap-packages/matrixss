###############################################################################
##
#W    test.g    The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
#H    Dev start : 2004-01-24 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/tst/test_g") := 
  "@(#)$Id$";;

RequirePackage("matrixss");;

SetAssertionLevel(1);
MATRIXSS_DEBUGLEVEL := 1;;
SetInfoLevel(MatrixSchreierSimsInfo, MATRIXSS_DEBUGLEVEL);
ReadTest("tst/test.tst");;
QUIT;;
