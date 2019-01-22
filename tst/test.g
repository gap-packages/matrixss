###############################################################################
##
#W    test.g    The Matrix Schreier-Sims package                
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-01-24 
##
###############################################################################

# Must not define method for Size, since we want to compare with GAP:s
# built-in results.
BindGlobal("MATRIXSS_TEST", true);;
#Read("tst/input/suz");
#Read("tst/input/tits");
LoadPackage("matrixss");;

SetAssertionLevel(2);
MATRIXSS_DEBUGLEVEL := 9;;
SetInfoLevel(MatrixSchreierSimsInfo, MATRIXSS_DEBUGLEVEL);
#Print("Testing standard deterministic algorithm\n");
#ReadTest("tst/standard.tst");;
Print("Testing standard probabilistic algorithm\n");
ReadTest("tst/random.tst");;
#Print("Testing nearly linear time algorithm\n");
#ReadTest("tst/linear.tst");;
QUIT;;

###############################################################################
#E
