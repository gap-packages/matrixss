###############################################################################
##
#W    standard.tst The Matrix Schreier-Sims package                
#W                 Deterministic algorithm test
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-01-24 
##
###############################################################################

gap> START_TEST("standard.tst");;

gap> MatrixSchreierSimsTest(2, 3, 2, 3);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : SimpleSchreierTree);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : AlternatingActions);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : CleverBasePoints);;
No order differences 

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : ShallowSchreierTree);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : ExtendSchreierTree);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : ExtendSchreierTree, SimpleSchreierTree);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : ExtendSchreierTree, AlternatingActions);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : ExtendSchreierTree, AlternatingActions, SimpleSchreierTree);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : ExtendSchreierTree, AlternatingActions, SimpleSchreierTree, ShallowSchreierTree);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : ExtendSchreierTree, SimpleSchreierTree, ShallowSchreierTree);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : ExtendSchreierTree, AlternatingActions, ShallowSchreierTree);;
No order differences

gap> MatrixSchreierSimsTest(2, 3, 2, 3 : AlternatingActions, ShallowSchreierTree);;
No order differences

gap> STOP_TEST("standard.tst", 0);;

###############################################################################
#E
