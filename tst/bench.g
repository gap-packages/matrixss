###############################################################################
##
#W    bench.g    The Matrix Schreier-Sims package                
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-03-22
##
###############################################################################

# Must not define method for Size, since we want to compare with GAP:s
# performance.
BindGlobal("MATRIXSS_TEST", true);;
LoadPackage("matrixss");;

SetAssertionLevel(0);
MATRIXSS_DEBUGLEVEL := 0;;
SetInfoLevel(MatrixSchreierSimsInfo, MATRIXSS_DEBUGLEVEL);

MATRIXSS_RandomSchreierSimsBenchmark := 
  function(startDeg, stopDeg, startField, stopField, maxReeSize, maxSuzukiSize)
    local groups, group, test_time, size, group_time, g;
    
    test_time := Runtime();
    groups := MATRIXSS_GetBenchmarkGroups(startDeg, stopDeg, startField, 
                      stopField, maxReeSize, maxSuzukiSize);
    
    Print("Group\t\tTime [ms]\n\n");
    
    for group in [1 .. Length(groups)] do
        #DefaultStabChainOptions.random := 1;
        Print(groups[group], "\t", MATRIXSS_TimedCall(Size, 
                [Group(GeneratorsOfGroup(groups[group]))]), "\n");        
        Unbind(groups[group]);
    od;
    
    Print("Total time for test : ", Runtime() - test_time, "\n");
    
    Print("Benchmark completed\n");
    return true;
end;;

ReadTest("tst/bench.tst");;

###############################################################################
#E
