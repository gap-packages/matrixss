###############################################################################
##
#W    test.gi     The Matrix Schreier Sims package - Test code                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
##    Dev start : 2004-07-01
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/lib/test_gi") := 
  "@(#)$Id$";

InstallGlobalFunction(MatrixGroupOrder, function(G)
    local ret, orbit, order;
    
    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    # Compute SGS and base and orbits (ie Schreier trees)
    ret := MatrixSchreierSims(G);
    
    # Compute order of group using computed orbit sizes
    order := 1;
    for orbit in ret[3] do
        order := order * Size(orbit);
    od;
    
    return order;
end);

InstallGlobalFunction(RandomMatrixGroupOrder, function(G)
    local ret, orbit, order;
    
    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    # Compute SGS and base and orbits (ie Schreier trees)
    ret := MatrixRandomSchreierSims(G, 99/100);
    
    # Compute order of group using computed orbit sizes
    order := 1;
    for orbit in ret[3] do
        order := order * Size(orbit);
    od;
    
    return order;
end);

# Creates a list of matrix groups to test the Schreier-Sims algorithm
#      maxFieldSize - maximum finite field size to be tested
#      maxMatrixSize - maximum degree of matrix groups
MATRIXSS_GetTestGroups := 
  function(maxFieldSize, maxMatrixSize) 
    local degree, power, primeNr, prime, groups, groupTypes, type;
    
    # Use the following group creation functions to make some test groups
    groupTypes := 
      Immutable([GeneralLinearGroup, SpecialLinearGroup, 
              GeneralOrthogonalGroup, SpecialOrthogonalGroup]);
#              GeneralUnitaryGroup, SpecialUnitaryGroup]);
    
    # List of test groups
    groups := [];
    
    for degree in [2 .. maxMatrixSize] do
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

# Creates a list of matrix groups to benchmark the Schreier-Sims algorithm
#      maxClassicalGroupFieldSize - maximum finite field size to use for the
#         classical matrix groups
#      maxClassicalGroupDegree - maximum matrix size for classical 
#         matrix groups
#      maxReeSize - maximum ReeGroup size
#      maxSuzukiSize - maximum SuzukiGroup size
MATRIXSS_GetBenchmarkGroups := 
  function(maxClassicalGroupFieldSize, maxClassicalGroupDegree,
          maxReeGroupSize, maxSuzukiSize)
  local groups, size;
    
    # Use all test groups as benchmark groups
    groups := MATRIXSS_GetTestGroups(maxClassicalGroupFieldSize,
                      maxClassicalGroupDegree);
    
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

InstallGlobalFunction(MatrixSchreierSimsTest, function(maxDegree, maxFieldSize)
    local groups, group, size1, size2;
    
    # Get list of test groups
    groups := MATRIXSS_GetTestGroups(maxDegree, maxFieldSize);
    
    # Compute order of all groups using GAP:s builtin Order and using our
    # Schreier-Sims algorithm
    for group in groups do
        Print("Checking group : ", group, "\n");
        
        size1 := Order(group);
        if ValueOption("UseRandomSS") <> fail then
            size2 := RandomMatrixGroupOrder(group);
        else
            size2 := MatrixGroupOrder(group);
        fi;
        
        if size1 <> size2 then
            Print("Correct order: ", size1, "\n");
            Print("Computed order: ", size2, "\n");
            Print("Order difference!\n");
        fi;
    od;
    
    Print("No order differences\n");
    return true;
end);

InstallGlobalFunction(MatrixSchreierSimsBenchmark, function(maxDegree, 
        maxFieldSize, maxReeSize, maxSuzukiSize)
    local groups, group, test_time, size, group_time, schreierSims;
        
    # Get list of benchmark groups
    groups := MATRIXSS_GetBenchmarkGroups(maxDegree, maxFieldSize, 
                      maxReeSize, maxSuzukiSize);
    
    Print("Group\t\tTime [ms]\n\n");
    
    if ValueOption("UseRandomSS") <> fail then
    
        # Start test timer
        test_time := Runtime();
        
        for group in groups do
            # Time each run of Schreier-Sims        
            group_time := Runtime();
            MatrixRandomSchreierSims(group, 99/100);
            
            Print(group, "\t", Runtime() - group_time, "\n");
        od;
        
        Print("Total time for test : ", Runtime() - test_time, "\n");
    else
        # Start test timer
        test_time := Runtime();
        
        for group in groups do
            # Time each run of Schreier-Sims        
            group_time := Runtime();
            MatrixSchreierSims(group);
            
            Print(group, "\t", Runtime() - group_time, "\n");
        od;
        
        Print("Total time for test : ", Runtime() - test_time, "\n");
    fi;
    
    Print("Benchmark completed\n");
    return true;
end);

#E
