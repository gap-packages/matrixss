###############################################################################
#1
#W    test.gi     The Matrix Schreier Sims package - Test code                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
#H    Dev start : 2004-07-01
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
##
## These are auxiliary functions for test and benchmark.
##
###############################################################################

Revision.("matrixss/lib/test_gi") := 
  "@(#)$Id$";

###############################################################################
##
#F MATRIXSS_GetTestGroups(maxDegree, maxFieldSize)
##
## Creates a list of classical matrix groups to use when testing the package.
## The groups are `GL', `SL', `GO', `SO', `GU' and `SU'.
## \beginitems
## `maxDegree' & maximum matrix size for classical matrix groups to be used
##             for testing
##
## `maxFieldSize' & maximum finite field size for classical matrix groups to be
##                used for testing
## \enditems
##
###############################################################################
MATRIXSS_GetTestGroups := 
  function(maxDegree, maxFieldSize) 
    local degree, power, primeNr, prime, groups, groupTypes, type;
    
    # Use the following group creation functions to make some test groups
    groupTypes := 
      Immutable([GeneralLinearGroup, SpecialLinearGroup, 
              GeneralOrthogonalGroup, SpecialOrthogonalGroup]);
#              GeneralUnitaryGroup, SpecialUnitaryGroup]);
    
    # List of test groups
    groups := [];
    
    for degree in [2 .. maxDegree] do
        primeNr := 1;
        while primeNr <= 168 and Primes[primeNr] <= maxFieldSize do
            prime := Primes[primeNr];
            
            power := 1;
            while Primes[primeNr]^power <= maxFieldSize do
                for type in groupTypes do
                    if (type = GO or type = SO) and IsEvenInt(degree) then
                        Add(groups, type(1, degree, prime^power));
                        Add(groups, type(-1, degree, prime^power));
                    else
                        Add(groups, type(degree, prime^power));
                    fi;
                od;
                
                power := power + 1;
            od;
            
            primeNr := primeNr + 1;
        od;
    od;
    
    return groups;
end;

###############################################################################
##
#F MATRIXSS_GetBenchmarkGroups(maxDegree, maxFieldSize)
##
## Creates a list of classical matrix groups and sporadic groups to use when 
## benchmarking the package.
## The classical groups are `GL', `SL', `GO' and `SO'. The sporadic groups are
## `Ree' and `Sz'.
## \beginitems
## `maxDegree' & maximum matrix size for classical matrix groups to be used
##             for testing
##
## `maxFieldSize' & maximum finite field size for classical matrix groups to be
##                used for testing
##
## `maxReeSize' & maximum `ReeGroup' size, see "ref:Ree" in the reference 
##                manual.
##
## `maxSuzukiSize' & maximum `SuzukiGroup' size, see "ref:Sz" in the reference 
##                   manual.
## \enditems
##
###############################################################################
MATRIXSS_GetBenchmarkGroups := 
  function(maxDegree, maxFieldSize, maxReeGroupSize, maxSuzukiSize)
  local groups, size;
    
    # Use all test groups as benchmark groups
    groups := MATRIXSS_GetTestGroups(maxDegree, maxFieldSize);
    
    # Also use some sporadic matrix groups as benchmark groups
    size := 1;
    while 3^(1 + 2 * size) <= maxReeGroupSize do
        Add(groups, ReeGroup(3^(1 + 2 * size)));
        size := size + 1;
    od;
    
    size := 3;
    while 2^size <= maxSuzukiSize do
        Add(groups, SuzukiGroup(2^size));
        size := size + 2;
    od;
    
    return groups;
end;

###############################################################################
##
#F MATRIXSS_TimedCall(call, args)
##
## Runs the specified function with arguments and return the running time in 
## milliseconds, as given by Runtime, see "ref:Runtime" in the reference 
## manual.
## \beginitems
## `call' & function to call
##
## `args' & list of arguments to `call'
## \enditems
##
###############################################################################
MATRIXSS_TimedCall := function(call, args)
    local time;
    
    time := Runtime();
    CallFuncList(call, args);
    return Runtime() - time;
end;

InstallGlobalFunction(MatrixSchreierSimsTest, function(maxDegree, maxFieldSize)
    local groups, group, size1, size2;
    
    # Get list of test groups
    groups := MATRIXSS_GetTestGroups(maxDegree, maxFieldSize);
    
    # Compute order of all groups using GAP:s builtin Order and using our
    # Schreier-Sims algorithm
    for group in groups do
        MATRIXSS_DebugPrint(1, ["Checking group : ", group]);
        
        size1 := Size(group);
        size2 := MatrixGroupOrderStabChain(StabChainMatrixGroup(group).
                         SchreierStructure);
        
        if size1 <> size2 then
            Print("Group: ", group, " Correct order: ", size1, 
                  "Computed order: ", size2, "\n");
        fi;
    od;
    
    Print("No order differences\n");
    return true;
end);

InstallGlobalFunction(MatrixSchreierSimsBenchmark, function(maxDegree, 
        maxFieldSize, maxReeSize, maxSuzukiSize)
    local groups, group, test_time, size, schreierSims;
        
    # Get list of benchmark groups
    groups := MATRIXSS_GetBenchmarkGroups(maxDegree, maxFieldSize, 
                      maxReeSize, maxSuzukiSize);
    
    Print("Group\t\tTime [ms]\n\n");
    
    # Start test timer
    test_time := Runtime();
    
    for group in groups do
        # Time each run of Schreier-Sims  
        Print(group, "\t", MATRIXSS_TimedCall(StabChainMatrixGroup, [group]), 
              "\n");
    od;
        
    Print("Total time for test : ", Runtime() - test_time, "\n");
    Print("Benchmark completed\n");
    return true;
end);

###############################################################################
#E
