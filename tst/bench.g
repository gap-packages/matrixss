###############################################################################
##
#W    bench.g    The Matrix Schreier-Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B��rnhielm
#H    Dev start : 2004-03-22
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/tst/bench_g") := 
  "@(#)$Id$";;

RequirePackage("matrixss");;

SetAssertionLevel(0);
MATRIXSS_DEBUGLEVEL := 0;;
SetInfoLevel(MatrixSchreierSimsInfo, MATRIXSS_DEBUGLEVEL);

MATRIXSS_RandomSchreierSimsBenchmark := 
  function(maxDegree, maxFieldSize, maxReeSize, maxSuzukiSize)
    local groups, group, test_time, size, group_time, g;
    
    test_time := Runtime();
    groups := MATRIXSS_GetBenchmarkGroups(maxDegree, maxFieldSize, 
                      maxReeSize, maxSuzukiSize);
    
    Print("Group\t\tTime [ms]\n\n");
    
    for group in groups do
        #Print("Checking group : ", group, "\n");
        
        group_time := Runtime();
        #DefaultStabChainOptions.random := 1;
        g := Group(GeneratorsOfGroup(group));
        Order(g);
        #Print("Time for this group : ", Runtime() - group_time, "\n");        
        Print(group, "\t", Runtime() - group_time, "\n");        
    od;
    
    Print("Total time for test : ", Runtime() - test_time, "\n");
    
    Print("Benchmark completed\n");
    return true;
end;;

ReadTest("tst/bench.tst");;
