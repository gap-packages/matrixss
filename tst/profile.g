###############################################################################
##
#W    profile.g    The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
#H    Dev start : 2004-03-22
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/tst/profile_g") := 
  "@(#)$Id$";;

RequirePackage("matrixss");;

SetAssertionLevel(0);
MATRIXSS_DEBUGLEVEL := 0;;
SetInfoLevel(MatrixSchreierSimsInfo, MATRIXSS_DEBUGLEVEL);
functionNames := Filtered(NamesUserGVars(), function(element)
    return IsMatchingSublist(element, "MATRIXSS_");
end);;
ProfileFunctions(Filtered(List(functionNames, EvalString), IsFunction));
#ProfileOperationsAndMethods(true);
ProfileOperations(true);
ProfileGlobalFunctions(true);

ReadTest("tst/profile.tst");;
DisplayProfile();
