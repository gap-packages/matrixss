###############################################################################
##
#W    code.gi     The Matrix Schreier Sims package                
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
##    Dev start : 2004-01-10 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/lib/code_gi") := 
  "@(#)$Id$";

MATRIXSS_DEBUGLEVEL := 0;

# An implementation of the Schreier-Sims algorithm, for matrix groups
InstallGlobalFunction(MatrixSchreierSims, function(G)
    # Local functions
    local DebugPrint, PointAction, ProjectiveAction, GetOrbit, IsPointInOrbit,
          SchreierTree, ComputeSchreierTree, ExtendSchreierTree, SchreierSims, 
          NewBasePoint, ProjectiveNewBasePoint, Membership, OrbitElement, 
          GetSchreierTreeEdge, GetOrbitSize, GetSchreierGenerator, 
          GetPartialBaseSGS, ExtendBase, NewBasePoint2, IsIdentity,
          ProjectiveIsIdentity, IsConstantList, CreateInitialSchreierTree,
          CopySchreierTree,
    # Local variables
          ssInfo, generators, base, trees, points, element, gens, list, level,
          normedPoints, levelStruct, hash;

    DebugPrint := function(level, message)
        #Info(MatrixSchreierSimsInfo, level, message);
        #Info(MatrixSchreierSimsInfo, level, message, "\n");
        if level <= MATRIXSS_DEBUGLEVEL then
            CallFuncList(Print, Concatenation(message, ["\n"]));
        fi;
    end;

    # The action of a group element (a matrix) on a point (a row vector)
    # The action is from the right
    PointAction := function(element, point)
        if not IsMatrix(element) or not IsRowVector(point) then
            DebugPrint(1, ["Element : ", element]);
            DebugPrint(1, ["Point : ", point]);
            Error("<element> must be a matrix and <point> must be a row vector");
        fi;
        
        return point * element;
    end;
    
    # The projective action of a matrix on a row vector
    # The one-dimensional subspace corresponding to the point is represented
    # by the corresponding normed row vector
    ProjectiveAction := function(element, point)
        return NormedRowVector(PointAction(element, 
                       NormedRowVector(point)));
    end;
    
    # Identity check when using the PointAction
    IsIdentity := function(element, identity)
        return element = identity;
    end;
    
    # Checks if all list elements are equal
    IsConstantList := function(list)
        local i;
        
        for i in [2 .. Length(list)] do
            if not list[i] = list[1] then
                return false;
            fi;
        od;
        return true;
    end;
    
    # Identity check when using projective action (all scalar matrices are
    # considered equal to the identity)
    ProjectiveIsIdentity := function(element, identity)
        if IsDiagonalMat(element) and 
           IsConstantList(DiagonalOfMat(element)) then
            return true;
        else
            return false;
        fi;
    end;
    
    # return all points (as a list) in the orbit of the point 
    # which is root of the schreier tree
    # The list is not necessarily sorted, and it is mutable
    GetOrbit := function(schreierTree)
        return HashKeyEnumerator(schreierTree);
    end;

    # Check if a given point is in the orbit defined by the given Schreier tree
    IsPointInOrbit := function(schreierTree, point)
        if not IsBool(GetHashEntry(schreierTree, point)) then
            return true;
        else
            return false;
        fi;    
    end;

    # Get the label of the edge originating at the given point, and directed 
    # towards the root of the given Schreier tree
    GetSchreierTreeEdge := function(schreierTree, point)
        return GetHashEntry(schreierTree, point);
    end;

    # Get size of orbit defined by the given Schreier tree
    GetOrbitSize := function(schreierTree)
        return Size(schreierTree);
    end;

    # Helper function for SchreierTree.
    ComputeSchreierTree := function(generators, tree, root, action)
        local newElement, childElements, generator, child, elements, element;
        
        elements := [root];
        repeat
            childElements := [];
            for element in elements do
                for generator in generators do
                    newElement := action(generator[1], element);
                    
                    # check that the element is not already in the 
                    # Schreier tree (do not create cycles)
                    if not IsPointInOrbit(tree, newElement) then
                        AddHashEntry(tree, newElement, generator);
                        Add(childElements, newElement);
                    fi;
                od;    
            od;
            elements := childElements;
        until IsEmpty(elements);
        
        return tree;
    end;      
    
    CreateInitialSchreierTree := function(root, hash, identity)
        local tree;
        
        # Create Schreier vector
        tree := SparseHashTable(hash);
        
        # Make the root point to itself 
        AddHashEntry(tree, root, [identity, identity]);
        
        return tree;
    end;
    
    CopySchreierTree := function(tree, hash)
        local copyTree, keys, value, key;
        
        copyTree := SparseHashTable(hash);
        keys := HashKeyEnumerator(tree);
        for key in keys do
            value := GetHashEntry(tree, key);
            #SetHashEntryAtLastIndex(copyTree, value);
            AddHashEntry(copyTree, key, value);
        od;
        
        return [copyTree, keys];
    end;
        
    # Extends an existing Schreier tree by a given set of generators
    ExtendSchreierTree := 
      function(oldTree, generators, oldGenerators, action, hash)
      local tree, point, generator, newPoint, newPoints, orbit, list;
      
      # Make the root point to itself 
      #AddHashEntry(oldTree, root, [identity, identity]);
      
      list := CopySchreierTree(oldTree, hash);
      tree := list[1];
      orbit := ShallowCopy(list[2]);
      
      DebugPrint(4, ["Old orbit: ", orbit]);
      DebugPrint(4, ["Old gens: ", oldGenerators]);
      DebugPrint(4, ["Gens    : ", generators]);
      
      repeat
          newPoints := [];
          for point in orbit do
              for generator in generators do
                  
                  # Add edges for all new points and new generators
                  if not IsPointInOrbit(oldTree, point) or
                     not generator in oldGenerators then
                      newPoint := action(generator[1], point);
                      
                      if not IsPointInOrbit(tree, newPoint) then
                          AddHashEntry(tree, newPoint, generator);
                          Add(newPoints, newPoint);
                      fi;
                  fi;
              od;
          od;
          orbit := newPoints;
      until IsEmpty(orbit);
      
      return tree;
  end;    

    # Computes a Schreier tree 
    # root - root of the Schreier tree
    # generators - generators for group
    # points - point set where root comes from
    # action - the action used to create the tree
    # hash - hash function to be used
    #
    # This is just a computation of a spanning tree for a connected component
    SchreierTree := 
      function(generators, points, root, action, identity, hash)
        local tree;
        
        # Create Schreier vector
        tree := SparseHashTable(hash);
        
        # Make the root point to itself 
        AddHashEntry(tree, root, [identity, identity]);
        
        # Fill Schreier vector
        return Immutable(ComputeSchreierTree(generators, tree, root, action));
    end;

    # Compute the group element that connects the root of the Schreier tree to
    # a given point
    # this function assumes that the point actually is in the orbit described by
    # the given Schreier tree
    OrbitElement := 
      function(schreierTree, point, action, identity, IsIdentity)
        local element, edge;
        
        # the group element and its inverse
        element := [identity, identity];
        
        repeat
            edge := GetSchreierTreeEdge(schreierTree, point);
            
            Assert(1, not IsBool(edge), "Point not in orbit!\n");
            
            if IsIdentity(edge[1], identity) then
                return element;
            fi;
            
            point := action(edge[2], point);
            element[1] := edge[1] * element[1];
            element[2] := element[2] * edge[2];
        until false;
    end;

    # check if an element belongs to a group, using sifting
    # ssInfo - main information structure about our base
    # element - the element to check membership for
    # identity - group identity
    Membership := 
      function(ssInfo, element, identity)
        local level, residue, word, point;
        
        residue := element;
        
        # Find an expression of element in terms of the generators of the
        # groups in our stabiliser chain, using the Schreier trees
        for level in [1 .. Length(ssInfo)] do
            DebugPrint(5, ["residue: ", residue[1], "\nbase: ", 
                    ssInfo[level].partialBase, "\naction", 
                    ssInfo[level].action]);
            point := ssInfo[level].action(residue[1], 
                             ssInfo[level].partialBase);
            
            if not IsPointInOrbit(ssInfo[level].schreierTree, point) then
                return [residue, level];
            fi;
            
            word := OrbitElement(ssInfo[level].schreierTree, point, 
                            ssInfo[level].action, identity, 
                            ssInfo[level].IsIdentity);
            residue[1] := residue[1] * word[2];
            residue[2] := word[1] * residue[2];
        od;
        
        return [residue, Length(ssInfo) + 1];
    end; 
        
    # Find a point not in base that is moved by element
    # (element fixes the base)
    NewBasePoint := function(element, action, identity, field)
        local basis, point, basePoint, i, j, length;
        
        DebugPrint(3, ["Matrix that fixes whole base: ", element]);
        DebugPrint(5, ["Point field: ", field]);
        field := AsList(field);
        
        length := Length(element);
        for i in [1 .. length] do
            for j in [1 .. length] do
                DebugPrint(8, ["Checking matrix element: ", element[i][j]]);           
                # If the element is not a diagonal matrix
                if not i = j and not IsZero(element[i][j]) then
                    basePoint := ZeroMutable(field[1]);
                    DebugPrint(6, ["Basepoint: ", basePoint]);
                    basePoint[i] := 
                      One(FieldOfMatrixGroup(Group(identity)));
                    DebugPrint(6, ["Basepoint: ", basePoint]);
                    return Immutable(basePoint);
                fi;
            od;
        od;
        
        for i in [1 .. length] do
            for j in [1 .. length] do
                
                # If the element is not a scalar matrix
                if not i = j and not element[i][i] = element[j][j] then
                    basePoint := ZeroMutable(field[1]);
                    basePoint[i] := 
                      One(FieldOfMatrixGroup(Group(identity)));
                    basePoint[j] := 
                      One(FieldOfMatrixGroup(Group(identity)));
                    return Immutable(basePoint);
                fi;
            od;
        od;
        
        # If the element is not the identity matrix
        if not element[1][1] = 
           One(FieldOfMatrixGroup(Group(identity))) then
            basePoint := ZeroMutable(field[1]);
            basePoint[1] := 
              One(FieldOfMatrixGroup(Group(identity)));
            return Immutable(basePoint);
        fi;
        
        Error("No new base point found!\n");
        return fail;
    end;

    # Create a Schreier generator for the stabiliser in the group which has 
    # "generator" as one of its generators. The stabiliser fixes "point".
    GetSchreierGenerator := 
      function(schreierTree, generator, point, action, identity, IsIdentity)
        local element1, element2, edge, inv_edge;
        
        element1 := OrbitElement(schreierTree, point, action,
                            identity, IsIdentity);
        element2 := OrbitElement(schreierTree, 
                            action(generator[1], point), action, identity,
                            IsIdentity);
        
        edge := element1[1] * generator[1] * element2[2];
        inv_edge := element2[1] * generator[2] * element1[2];
        
        return [edge, inv_edge];
    end;
    
    
    # Add a new base point to the base, so that a given element is not in the
    # stabiliser of the point
    # ssInfo - main information structure for the current Schreier-Sims run
    # badElement - the element that fixes all current base points
    # identity - the group identity
    ExtendBase := function(ssInfo, badElement, identity)
        local newPoint, length, levelStruct;
        
        DebugPrint(3, ["Finding new base point"]);
        
        length := Length(ssInfo);
        
        # Find new base point
        newPoint := NewBasePoint(badElement[1],
                            ssInfo[length].action, 
                            identity, ssInfo[length].points);
        
        # Extend base
        levelStruct := rec(
                           partialBase := newPoint,
                           action := PointAction,
                           points := ssInfo[1].points,
                           hash := ssInfo[length].hash,
                           schreierTree := 
                           CreateInitialSchreierTree(newPoint, 
                                   ssInfo[length].hash, identity),
                           oldSGS := AsSortedList([]),
                           IsIdentity := IsIdentity);
        Add(ssInfo, levelStruct); 

        levelStruct := rec(
                           partialBase := NormedRowVector(newPoint),
                           action := ProjectiveAction,
                           points := ssInfo[length].points,
                           hash := ssInfo[length].hash,
                           schreierTree := 
                           CreateInitialSchreierTree(NormedRowVector(newPoint),
                                   ssInfo[length].hash, identity),
                           oldSGS := AsSortedList([]),
                           IsIdentity := ProjectiveIsIdentity);
        Add(ssInfo, levelStruct); 
    end;
    
    # The main Schreier-Sims function
    # ssInfo - main information structure for the current Schreier-Sims run
    # partialSGS - given partial strong generating set
    # level - the level of the call to Schreier-Sims
    # identity - the group identity
    SchreierSims := function(ssInfo, partialSGS, level, identity)
        local generator, point, orbit, strip, schreierGenerator, element, 
              action, recursiveLevel, schreierTree, SGS, oldSGS, points, 
              newPoint, oldSchreierTree, field, newBasePoint, oldOrbit;
        
        action := ssInfo[level].action;
        field  := ssInfo[level].points;
        
        DebugPrint(6, ["Schreier-Sims at level ", level]);
        
        # Find the generators from the partial SGS that fixes all points at
        # lower levels.
        SGS := ShallowCopy(partialSGS);
        for element in partialSGS do
            for point in ssInfo{[1 .. level - 1]} do
                if not action(element[1], point.partialBase) = 
                   point.partialBase then
                    RemoveSet(SGS, element);
                    break;
                fi;
            od;
        od;
        MakeImmutable(SGS);
                
        DebugPrint(4, ["Base point : ", ssInfo[level].partialBase]);
        DebugPrint(9, ["Hash func : ", ssInfo[level].hash]);
        
        # Compute schreier tree for current level
        oldSchreierTree := ssInfo[level].schreierTree;
        ssInfo[level].schreierTree := 
          ExtendSchreierTree(ssInfo[level].schreierTree, 
                  SGS, ssInfo[level].oldSGS, action, ssInfo[level].hash);
        
        orbit := Immutable(GetOrbit(ssInfo[level].schreierTree));
        
        DebugPrint(2, ["New Schreier Tree : ", ssInfo[level].schreierTree]);
        DebugPrint(2, ["Old Schreier Tree : ", oldSchreierTree]);
        Assert(3, not IsIdenticalObj(ssInfo[level].schreierTree, 
                oldSchreierTree));
        DebugPrint(4, ["Orbit size for level ", level, " is ", 
                Length(orbit)]); 
        
        # We now want to make sure that SGS also fixes the current level
        
        for point in orbit do
            for generator in SGS do
                
                # Avoid rechecking Schreier generators
                if not IsPointInOrbit(oldSchreierTree, point) or 
                   not generator in ssInfo[level].oldSGS then
                    
                    # Compute Schreier generator for current level
                    schreierGenerator := 
                      GetSchreierGenerator(ssInfo[level].schreierTree,
                              generator, point, action, identity,
                              ssInfo[level].IsIdentity);
                                        
                    DebugPrint(4, ["Schreier Generator : ", 
                            schreierGenerator]);
                    
                    Assert(3, not schreierGenerator[1] = identity, 
                           "Identity Schreier generator!");
                    if ssInfo[level].IsIdentity(schreierGenerator[1], 
                               identity) then
                        continue;
                    fi;
                    
                    DebugPrint(4, ["Schreier Generator : ", 
                            schreierGenerator]);
                    
                    # Check if Schreier generator is in stabiliser at 
                    # the current level
                    points := [level + 1 .. Length(ssInfo)];
                    strip := Membership(ssInfo{points},
                                     schreierGenerator, 
                                     identity);
                                        
                    # The drop-out level is in range 
                    # [1 .. Length(ssInfo.partialBase) - level] 
                    # but we want the range given by points
                    strip[2] := strip[2] + level;
                    
                    DebugPrint(3, ["Dropout level : ", strip[2]]);
                    
                    if not strip[1][1] = identity then
                        DebugPrint(3, ["Residue found"]);
                        
                        # We have found a Schreier generator which is not in
                        # the stabiliser of the current level, and so we must
                        # add the residue of this generator to our partial SGS
                        # in order to make it into a real SGS
                        
                        # Add residue to partial SGS
                        # This makes some levels incomplete and so we must
                        # recompute them recursively
                        AddSet(partialSGS, Immutable(strip[1]));
                        AddSet(partialSGS, Immutable(Reversed(strip[1])));
                        
                        # Possibly extend the base if the Schreier generator
                        # fixes all points in our base
                        if strip[2] = Length(ssInfo) + 1 then
                            ExtendBase(ssInfo, strip[1], identity);
                        fi;
                        
                        # We must not recompute all levels downward from the
                        # dropout level
                        for recursiveLevel in Reversed([level + 1 .. 
                                strip[2] + 1]) do
                            oldSGS := ssInfo[level].oldSGS;
                            ssInfo[level].oldSGS := SGS;
                            SchreierSims(ssInfo, partialSGS,
                                    recursiveLevel, identity);
                            ssInfo[level].oldSGS := oldSGS;
                        od;
                    fi;
                fi;
            od;
        od;
        
        ssInfo[level].oldSGS := SGS;
    end;
    
    # Construct a partial base and a partial SGS given a set of generators
    # for a group.
    # generators - given set of generators
    # action - action to use when finding new base points
    # IsIdentity - function to use when checking that a group element is the
    # identity (under the given action)
    # field - the vector space on which the group acts
    # identity - the group identity
    GetPartialBaseSGS := 
      function(generators, action, IsIdentity, field, identity)
      local bae, sgs, point, nFixedPoints, newPoint, element, list, 
            points, newSGS;
        
        points := [];
        newSGS := [];
        
        # we make a partial strong generating set which also contain
        # inverses of all elements
        for element in generators do
            if not IsIdentity(element, identity) then
                nFixedPoints := 0;
                for point in points do
                    if action(element, point) = point then
                        nFixedPoints := nFixedPoints + 1;
                    else
                        break;
                    fi;
                od;
                list := [element, Inverse(element)];
                
                DebugPrint(2, ["Matrix ", element, " fixes all points ", 
                        points]);

                if nFixedPoints = Length(points) then
                    newPoint := NewBasePoint(list[1], 
                                        action, identity, field);
                    Add(points, newPoint);
                fi;
                
                # Save reference to generator and its inverse
                # Then inverses need not be calculated later
                AddSet(newSGS, Immutable(list));
                AddSet(newSGS, Immutable(Reversed(list)));
            fi;
        od;
        
        return [points, newSGS];
    end;


    ### MAIN Schreier-Sims 

    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    # Get initial set of generators, to be extended to a partial SGS
    gens := GeneratorsOfGroup(G);
    
    # The vector space on which the group acts
    points := FullRowSpace(FieldOfMatrixGroup(G), DimensionOfMatrixGroup(G));
    
    # Compute initial partial SGS and base
    list := GetPartialBaseSGS(gens, PointAction, IsIdentity, points, 
                    Identity(G));
    
    # Set of lines, to be used by projective action
    normedPoints := NormedRowVectors(points);
    
    # Initial base and partial SGS
    generators := list[2];
    base := list[1];
    
    DebugPrint(3, ["Partial base : ", base]);
    DebugPrint(3, ["Partial sgs : ", generators]);

    # Main structure holding information needed by the algorithm
    ssInfo := [];
    
    hash := SparseIntKey(points, base[1]);
    Assert(1, not IsBool(hash));
    
    # Fill ssInfo with initial data
    for level in [1 .. Length(base)] do
        levelStruct := 
          rec(
              partialBase := base[level],
              action := PointAction,
              points := points,
              hash := hash,
              schreierTree := 
              CreateInitialSchreierTree(base[level], hash, Identity(G)),
              oldSGS := AsSortedList([]),
              IsIdentity := IsIdentity);
        Add(ssInfo, levelStruct); 
        
        levelStruct := 
          rec(
              partialBase := NormedRowVector(base[level]),
              action := ProjectiveAction,
              points := normedPoints,
              hash := hash,
              schreierTree := 
              CreateInitialSchreierTree(NormedRowVector(base[level]), hash,
                      Identity(G)),
              oldSGS := AsSortedList([]),
              IsIdentity := ProjectiveIsIdentity);
        Add(ssInfo, levelStruct); 
    od;
    
    
    DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    
    # Call Schreier-Sims algorithm for each level (starting from top)
    for level in Reversed([1 .. Length(ssInfo)]) do
        SchreierSims(ssInfo, generators, level, Identity(G));
    od;
    
    DebugPrint(2, ["Matrix Schreier-Sims done"]);
    
    # Create output structure
    list := [[], generators, []];
    for levelStruct in ssInfo do
        Add(list[1], levelStruct.partialBase);
        Add(list[3], levelStruct.schreierTree);
    od;
    
    return Immutable(list);
end);


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

# Creates a list of matrix groups to test the Schreier-Sims algorithm
#      maxFieldSize - maximum finite field size to be tested
#      maxMatrixSize - maximum degree of matrix groups
MATRIXSS_GetTestGroups := 
  function(maxFieldSize, maxMatrixSize) 
    local degree, power, primeNr, prime, groups, groupTypes, type;
    
    # Use the following group creation functions to make some test groups
    groupTypes := 
      Immutable([GeneralLinearGroup, SpecialLinearGroup, 
              GeneralOrthogonalGroup, SpecialOrthogonalGroup,
              GeneralUnitaryGroup, SpecialUnitaryGroup]);
    
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
        size2 := MatrixGroupOrder(group);
            
        if not size1 = size2 then
            Print("Correct order: ", size1, "\n");
            Print("Computed order: ", size2, "\n");
            Print("Order difference!\n");
        fi;
    od;
    
    Print("No order differences\n");
    return true;
end);

InstallGlobalFunction(MatrixSchreierSimsBenchmark, function(maxDegree, maxFieldSize, maxReeSize, maxSuzukiSize)
    local groups, group, test_time, size, group_time;
        
    # Get list of benchmark groups
    groups := MATRIXSS_GetBenchmarkGroups(maxDegree, maxFieldSize, 
                      maxReeSize, maxSuzukiSize);
    
    Print("Group\t\tTime [ms]\n\n");
    
    # Start test timer
    test_time := Runtime();
    
    for group in groups do
        # Time each run of Schreier-Sims        
        group_time := Runtime();
        MatrixSchreierSims(group);
        
        Print(group, "\t", Runtime() - group_time, "\n");
    od;

    Print("Total time for test : ", Runtime() - test_time, "\n");
    
    Print("Benchmark completed\n");
    return true;
end);

#E