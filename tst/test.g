###############################################################################
##
#W    test.g    The Matrix Schreier-Sims package                
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

Revision.("matrixss/tst/test_g") := 
  "@(#)$Id$";;

# Must not define method for Size, since we want to compare with GAP:s
# built-in results.
BindGlobal("MATRIXSS_TEST", true);;
LoadPackage("matrixss");;

SetAssertionLevel(2);
MATRIXSS_DEBUGLEVEL := 0;;
SetInfoLevel(MatrixSchreierSimsInfo, MATRIXSS_DEBUGLEVEL);
Print("Testing standard deterministic algorithm\n");
ReadTest("tst/standard.tst");;
Print("Testing standard probabilistic algorithm\n");
ReadTest("tst/random.tst");;
Print("Testing nearly linear time algorithm\n");
ReadTest("tst/linear.tst");;
QUIT;;

###############################################################################
#E
