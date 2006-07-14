###############################################################################
#1
#W    test.gi     The Matrix Schreier Sims package - Test code                
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
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

MATRIXSS_GetRandomGenSet :=
  function(field, dimension, length)
    local generators;
    
    generators := [];
    while length > 0 do
        Add(generators, RandomInvertibleMat(dimension, field));
        length := length - 1;
    od;
    return Immutable(generators);
end;

MATRIXSS_BenchmarkRandomSets :=
  function(maxDegree, maxFieldSize, maxGenSetSize, nRunsPerType) 
    local degree, power, primeNr, prime, length, runs, generators, field, list;
    
    list := [];
    primeNr := 1;
    while primeNr <= 168 and Primes[primeNr] <= maxFieldSize do
        prime := Primes[primeNr];
        
        power := 1;
        while Primes[primeNr]^power <= maxFieldSize do
            
            field := GaloisField(prime, power);
            for degree in [2 .. maxDegree] do

                length := 1;
                while length <= maxGenSetSize do
                    
                    runs := 1;
                    Add(list, [field, degree]);
                    while runs <= nRunsPerType do
                        generators := 
                          MATRIXSS_GetRandomGenSet(field, degree, length);
                        Add(list[Length(list)], generators);
                        runs := runs + 1;
                    od;                    
                    length := length + 1;
                od;
            od;
            power := power + 1;
        od;
        primeNr := primeNr + 1;
    od;
    
    return list;
end;
                    
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
  function(startDeg, stopDeg, startField, stopField) 
    local degree, power, primeNr, prime, groups, groupTypes, type, gl, group,
          conj, tits;
    
    # Use the following group creation functions to make some test groups
    groupTypes := 
      Immutable([GeneralLinearGroup, SpecialLinearGroup, 
              GeneralOrthogonalGroup, SpecialOrthogonalGroup,
              GeneralUnitaryGroup, SpecialUnitaryGroup]);
    
    # List of test groups
    groups := [];
    
    for degree in [startDeg .. stopDeg] do
        primeNr := 1;
        while primeNr <= 168 and Primes[primeNr] <= stopField and 
          Primes[primeNr] >= startField do
            prime := Primes[primeNr];
            
            power := 1;
            while Primes[primeNr]^power <= stopField and 
              Primes[primeNr]^power >= startField do
                gl := GL(degree, prime^power);
                
                for type in groupTypes do
                    if (type = GO or type = SO) and IsEvenInt(degree) then
                        group := type(1, degree, prime^power);
                        conj := PseudoRandom(gl);
                        Add(groups, group^conj);
                        
                        group := type(-1, degree, prime^power);
                        conj := PseudoRandom(gl);
                        Add(groups, group^conj);
                    else
                        group := type(degree, prime^power);
                        conj := PseudoRandom(gl);
                        Add(groups, group^conj);
                    fi;
                od;
                
                power := power + 1;
            od;
            
            primeNr := primeNr + 1;
        od;
    od;
    
    #tits := TitsGroup();
    #Add(groups, tits[1]);
    
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
  local groups, size, group, gl;
    
    # Use all test groups as benchmark groups
    groups := MATRIXSS_GetTestGroups(maxDegree, maxFieldSize);
    
    # Also use some sporadic matrix groups as benchmark groups
    size := 1;
    while 3^(1 + 2 * size) <= maxReeGroupSize do
        group := ReeGroup(3^(1 + 2 * size));
        gl := GL(7, 3^(1 + 2 * size));
        
        Add(groups, group^PseudoRandom(gl));
        size := size + 1;
    od;
    
    size := 1;
    while 2^size <= maxSuzukiSize do
        group := SuzukiGroup(2^(1 + 2 * size));
        gl := GL(4, 2^(1 + 2 * size));
        
        Add(groups, group^PseudoRandom(gl));
        size := size + 1;
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

InstallGlobalFunction(MatrixSchreierSimsSetBenchmark, 
        function(maxDegree, maxFieldSize, maxLength, nRunsPerType)
        local sets, time, i, j;
    
    sets := MATRIXSS_BenchmarkRandomSets(maxDegree, maxFieldSize, maxLength,
                    nRunsPerType);
    
    for i in [1 .. Length(sets)] do
        time := 0;
        for j in [3 .. Length(sets[i])] do
            MATRIXSS_DebugPrint(1, ["Group with gen set ", sets[i][j]]);
            MATRIXSS_DebugPrint(1, ["Group identity : ", 
                    One(GL(sets[i][2], sets[i][1]))]);
            time := time + MATRIXSS_TimedCall(Size,
                            [Group(sets[i][j], One(GL(sets[i][2], 
                                    sets[i][1])))]);
        od;
        time := time / (j - 2);
        time := Int(time + 1/2);
        Print(Size(sets[i][1]), "\t", sets[i][2], "\t", 
              Length(sets[i][3]), "\t", time, "\n");
    od;
    Print("Benchmark completed\n");
end);             

InstallGlobalFunction(MatrixSchreierSimsTest, function(startDeg, stopDeg, 
        startField, stopField)
    local groups, group, size1, size2;
    
    # Get list of test groups
    groups := MATRIXSS_GetTestGroups(startDeg, stopDeg, startField, stopField);
    
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

InstallGlobalFunction(MatrixSchreierSimsBenchmark, function(startDeg, stopDeg,
        startField, stopField, maxReeSize, maxSuzukiSize)
    local groups, group, test_time, size, schreierSims;
        
    # Get list of benchmark groups
    groups := MATRIXSS_GetBenchmarkGroups(startDeg, stopDeg, startField, 
                      stopField, maxReeSize, maxSuzukiSize);
    
    Print("Group\t\tTime [ms]\n\n");
    
    # Start test timer
    test_time := Runtime();
    
    for group in [1 .. Length(groups)] do
        # Time each run of Schreier-Sims  
        Print(groups[group], "\t", 
              MATRIXSS_TimedCall(StabChainMatrixGroup, [groups[group]]), 
              "\n");
        Unbind(groups[group]);
    od;
        
    Print("Total time for test : ", Runtime() - test_time, "\n");
    Print("Benchmark completed\n");
    return true;
end);

###############################################################################
#E
