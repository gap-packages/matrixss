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
SetAssertionLevel(0);
MATRIXSS_DEBUGLEVEL := 1;;
SetInfoLevel(MatrixSchreierSimsInfo, MATRIXSS_DEBUGLEVEL);
#ProfileFunctions(["MATRIXSS_ComputeSchreierTree",
#        "MATRIXSS_GetOrbit",
#        "MATRIXSS_GetOrbitSize",
#        "MATRIXSS_GetSchreierTreeEdge",
#        "MATRIXSS_GetSchreierTrees",
#        "MATRIXSS_IsPointInOrbit",
#        "MATRIXSS_MSSAction",
#        "MATRIXSS_Membership2",
#        "MATRIXSS_NewBasePoint",
#        "MATRIXSS_OrbitElement",
#        "MATRIXSS_SchreierSims",
#        "MATRIXSS_SchreierTree",
#        "MATRIXSS_Stabiliser"]);
ProfileOperationsAndMethods(true);
ProfileGlobalFunctions(true);
ReadTest("tst/matrixss.tst");;
DisplayProfile();
#QUIT;;
