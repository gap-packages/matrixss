###############################################################################
##
#W    profile.g    The Matrix Schreier-Sims package                
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-03-22
##
###############################################################################

BindGlobal("MATRIXSS_PROFILE", true);;
LoadPackage("matrixss");;

SetAssertionLevel(0);
MATRIXSS_DEBUGLEVEL := 2;;
SetInfoLevel(MatrixSchreierSimsInfo, MATRIXSS_DEBUGLEVEL);
functionNames := Filtered(NamesUserGVars(), function(element)
    return IsMatchingSublist(element, "MATRIXSS_");
end);;
ProfileFunctions(Filtered(List(functionNames, EvalString), IsFunction));
        
#ProfileOperations(true);
ProfileGlobalFunctions(true);

ReadTest("tst/standard.tst");;
DisplayProfile();

###############################################################################
#E
