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

MATRIXSS_DEBUGLEVEL := 3;

# An implementation of the Schreier-Sims algorithm, for matrix groups
InstallGlobalFunction(MatrixSchreierSims, function(G)
    # Local functions
    local DebugPrint, PointAction, ProjectiveAction, GetOrbit, IsPointInOrbit,
          SchreierTree, ComputeSchreierTree, ExtendSchreierTree, SchreierSims, 
          NewBasePoint, ProjectiveNewBasePoint, Membership, OrbitElement, 
          GetSchreierTreeEdge, GetOrbitSize, GetSchreierGenerator, 
          GetPartialBaseSGS, ExtendBase, NewBasePoint2, IsIdentity,
          ProjectiveIsIdentity,
          # Local variables
          ssInfo, generators, base, trees, points, element, gens, list, level,
          normedPoints;

    DebugPrint := function(level, message)
        #Info(MatrixSchreierSimsInfo, level, message);
        #Info(MatrixSchreierSimsInfo, level, message, "\n");
        if level <= MATRIXSS_DEBUGLEVEL then
            CallFuncList(Print, Concatenation(message, ["\n"]));
        fi;
    end;

    # The action of a group element g (a matrix) on a point p (a row vector)
    # The action is from the right
    # Used in Matrix Schreier-Sims
    PointAction := function(element, point)
        
        if not IsMatrix(element) or not IsRowVector(point) then
            DebugPrint(1, ["Element : ", element]);
            DebugPrint(1, ["Point : ", point]);
            Error("<element> must be a matrix and <point> must be a row vector");
        fi;
        
        return point * element;
    end;

    ProjectiveAction := function(element, point)
        return NormedRowVector(PointAction(element, 
                       NormedRowVector(point)));
    end;
    
    IsIdentity := function(element, identity)
        return element = identity;
    end;
    
    ProjectiveIsIdentity := function(element, identity)
        local order;
        
        order := ProjectiveOrder(element);
        return order[1] = 1;
    end;
    
    # return all points (as a list) in the orbit of the point 
    # which is root of the schreier tree
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
                    
                    # check that the element is not already in the Schreier tree
                    # (do not create cycles)
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

    ExtendSchreierTree := 
      function(oldTree, root, generators, oldGenerators, action, points, 
              identity)
        local tree, point, generator, newPoint, newPoints, oldOrbit, orbit;
        
        # Make the root point to itself 
        AddHashEntry(oldTree, root, [identity, identity]);

        oldOrbit := GetOrbit(oldTree);
        tree := oldTree;
        orbit := GetOrbit(tree);
        
        DebugPrint(4, ["Old orbit: ", oldOrbit]);
        DebugPrint(4, ["Old gens: ", oldGenerators]);
        DebugPrint(4, ["Gens    : ", generators]);
        
        repeat
            newPoints := [];
            for point in orbit do
                for generator in generators do
                    if not point in oldOrbit or
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
    # schreierTrees - trees with the base points as roots, and with generators from
    # the stabiliser of the previos point
    # base - base for the group under consideration
    # sgs - strong generating set for the stabiliser chain corresponding to base
    # element - the element to check membership for
    # action - the action to use
    Membership := 
      function(schreierTrees, base, element, action, identity, IsIdentity)
        local level, residue, word, point;
        
        residue := element;
        for level in [1 .. Length(base)] do
            DebugPrint(5, ["residue: ", residue[1], "\nbase: ", 
                    base[level], "\naction", action[level]]);
            point := action[level](residue[1], base[level]);
            
            if not IsPointInOrbit(schreierTrees[level], point) then
                return [residue, level];
            fi;
            
            word := OrbitElement(schreierTrees[level], point, 
                            action[level], identity, IsIdentity[level]);
            residue[1] := residue[1] * word[2];
            residue[2] := word[1] * residue[2];
        od;
        
        return [residue, Length(base) + 1];
    end; 
    
    # find a point not in base that is moved by element
    # (element fixes the base)
    NewBasePoint2 := function(base, element, action, identity, field)
        local basis, point, basePoint, i, j, dimension;
        
        basePoint := fail;
        if IsVectorSpace(field) then
            dimension := Dimension(field);
            basis := BasisVectors(CanonicalBasis(field));
        else
            dimension := Length(field[1]);
            for i in [1 .. dimension] do
                point := ZeroMutable(field[1]);
                point[i] := One(FieldOfMatrixGroup(Group(identity)));
                Add(basis, point);
            od;
        fi;
        
        DebugPrint(2, ["Basis vectors : ", basis]);
        
        for point in basis do
            if not action(element, point) = point and 
               not point in base then
                basePoint := point;
                break;
            fi;
        od;
        
        if IsBool(basePoint) then
            Error("No new base point found!\n");
        fi;
        
        return basePoint;
    end;
    
    # find a point not in base that is moved by element
    # (element fixes the base)
    NewBasePoint := function(base, element, action, identity, field)
        local basis, point, basePoint, i, j, length;
        
        DebugPrint(2, ["Matrix that fixes whole base: ", element]);
        DebugPrint(5, ["Point field: ", field]);
        field := AsList(field);
        
        length := Length(element);
        for i in [1 .. length] do
            for j in [1 .. length] do
                DebugPrint(5, ["Checking matrix element: ", element[i][j]]);           
                if not i = j and not IsZero(element[i][j]) then
                    basePoint := ZeroMutable(field[1]);
                    DebugPrint(2, ["Basepoint: ", basePoint]);
                    basePoint[i] := 
                      One(FieldOfMatrixGroup(Group(identity)));
                    DebugPrint(2, ["Basepoint: ", basePoint]);
                    #basePoint := NormedRowVector(basePoint);
                    #Assert(1, not action(element, basePoint) = basePoint);
                    return Immutable(basePoint);
                fi;
            od;
        od;
        
        for i in [1 .. length] do
            for j in [1 .. length] do
                if not i = j and not element[i][i] = element[j][j] then
                    basePoint := ZeroMutable(field[1]);
                    basePoint[i] := 
                      One(FieldOfMatrixGroup(Group(identity)));
                    basePoint[j] := 
                      One(FieldOfMatrixGroup(Group(identity)));
                    #basePoint := NormedRowVector(basePoint);
                    #Assert(1, not action(element, basePoint) = basePoint);
                    return Immutable(basePoint);
                fi;
            od;
        od;
        
        if not element[1][1] = 
           One(FieldOfMatrixGroup(Group(identity))) then
            basePoint := ZeroMutable(field[1]);
            basePoint[1] := 
              One(FieldOfMatrixGroup(Group(identity)));
            #basePoint := NormedRowVector(basePoint);
            #Assert(1, not action(element, basePoint) = basePoint);
            return Immutable(basePoint);
        fi;
        
        Error("No new base point found!\n");
        return fail;
    end;

    ProjectiveNewBasePoint := 
      function(base, element, action, identity, field)
        local point;
        
        point := NormedRowVector(NewBasePoint(base, element, action, 
                         identity, field));
        Assert(1, not point in base, "Point is already in base!");
        Assert(1, not action(element, point) = point, "Point is fixed!");
        return point;
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
    
    ExtendBase := function(ssInfo, badElement, identity)
        local newPoint, length;
        
        DebugPrint(3, ["Finding new base point"]);
        
        length := Length(ssInfo.partialBase);
        
        newPoint := ProjectiveNewBasePoint(ssInfo.partialBase, 
                            badElement[1],
                            ssInfo.action[length], 
                            identity, ssInfo.points[length]);
        
        Add(ssInfo.partialBase, newPoint);
        Add(ssInfo.action, ProjectiveAction);
        Add(ssInfo.oldSGS, []);
        Add(ssInfo.points, NormedRowVectors(ssInfo.points[1]));
        Add(ssInfo.newBasePoint, ProjectiveNewBasePoint);
        Add(ssInfo.hash, ssInfo.hash[length]);
        Add(ssInfo.schreierTrees, SparseHashTable(ssInfo.hash[length]));    
        Add(ssInfo.IsIdentity, ProjectiveIsIdentity);
    end;
                          
    SchreierSims := function(ssInfo, level, identity)
        local generator, point, orbit, strip, schreierGenerator, element, 
              action, recursiveLevel, schreierTree, SGS, oldSGS, points, 
              newPoint, oldSchreierTree, field, newBasePoint, oldOrbit;
        
        action := ssInfo.action[level];
        field  := ssInfo.points[level];
        newBasePoint := ssInfo.newBasePoint[level];
        
        DebugPrint(1, ["Schreier-Sims at level ", level]);
        
        SGS := ShallowCopy(ssInfo.partialSGS);
        for element in ssInfo.partialSGS do
            for point in ssInfo.partialBase{[1 .. level - 1]} do
                if not action(element[1], point) = point then
                    RemoveSet(SGS, element);
                    break;
                fi;
            od;
        od;
        MakeImmutable(SGS);
        
        DebugPrint(4, ["Base point : ", ssInfo.partialBase[level]]);
        DebugPrint(4, ["Hash func : ", ssInfo.hash[level]]);
        
        oldOrbit := Immutable(GetOrbit(ssInfo.schreierTrees[level]));
        #ssInfo.schreierTrees[level] := 
        #  SchreierTree(SGS, field, ssInfo.partialBase[level], 
        #          action, identity, ssInfo.hash[level]);
        ssInfo.schreierTrees[level] := 
          ExtendSchreierTree(ssInfo.schreierTrees[level], 
                  ssInfo.partialBase[level], SGS, ssInfo.oldSGS[level], 
                  action, field, identity);
        
        orbit := Immutable(GetOrbit(ssInfo.schreierTrees[level]));
        
        DebugPrint(4, ["Orbit size for level ", level, " is ", 
                GetOrbitSize(ssInfo.schreierTrees[level])]);
        
        for point in orbit do
            for generator in SGS do
                if not point in oldOrbit or 
                   not generator in ssInfo.oldSGS[level] then
                    schreierGenerator := 
                      GetSchreierGenerator(ssInfo.schreierTrees[level],
                              generator, point, action, identity,
                              ssInfo.IsIdentity[level]);
                    
                    Assert(3, not schreierGenerator[1] = identity, 
                           "Identity Schreier generator!");
                    if ssInfo.IsIdentity[level](schreierGenerator[1], 
                               identity) then
                        continue;
                    fi;
                    
                    DebugPrint(4, ["Schreier Generator : ", 
                            schreierGenerator]);
                    
                    points := [level + 1 .. Length(ssInfo.partialBase)];
                    strip := Membership(ssInfo.schreierTrees{points},
                                     ssInfo.partialBase{points},
                                     schreierGenerator, ssInfo.action{points}, 
                                     identity, ssInfo.IsIdentity{points});
                    
                    # The drop-out level is in range 
                    # [1 .. Length(ssInfo.partialBase) - level] 
                    # but we want the range given by points
                    strip[2] := strip[2] + level;
                    
                    DebugPrint(3, ["Dropout level : ", strip[2]]);
                    
                    if not strip[1][1] = identity then
                        DebugPrint(3, ["Residue found"]);
                        
                        # Add residue to partial SGS
                        # This makes some levels incomplete and so we must
                        # recompute them recursively
                        AddSet(ssInfo.partialSGS, Immutable(strip[1]));
                        AddSet(ssInfo.partialSGS, Immutable(Reversed(strip[1])));
                        
                        if strip[2] = Length(ssInfo.partialBase) + 1 then
                            ExtendBase(ssInfo, strip[1], identity);
                        fi;
                        
                        for recursiveLevel in Reversed([level + 1 .. strip[2]]) do
                            oldSGS := ssInfo.oldSGS[level];
                            ssInfo.oldSGS[level] := SGS;
                            SchreierSims(ssInfo, 
                                    recursiveLevel, identity);
                            ssInfo.oldSGS[level] := oldSGS;
                        od;
                    fi;
                fi;
            od;
        od;
        
        ssInfo.oldSGS[level] := SGS;
    end;

    GetPartialBaseSGS := 
      function(base, sgs, generators, action, IsIdentity, field, 
              identity)
      local gens, point, nFixedPoints, newPoint, element, list, 
            points, newSGS;
        
        points := ShallowCopy(base);
        gens := UnionSet(generators, sgs);
        #RemoveSet(gens, identity);
        newSGS := [];
        
        # we make a partial strong generating set which also contain
        # inverses of all elements
        for element in gens do
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
                    newPoint := ProjectiveNewBasePoint(points, list[1], 
                                        action, identity, field);
                    Add(points, newPoint);
                fi;
                
                # Save reference to generator and its inverse
                # Then inverses need not be calculated later
                AddSet(newSGS, Immutable(list));
                #if not list[1] = list[2] then
                AddSet(newSGS, Immutable(Reversed(list)));
                #fi;
            fi;
        od;
        
        return [points, newSGS];
    end;


    ### MAIN Schreier-Sims 

    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    gens := GeneratorsOfGroup(G);
    points := FullRowSpace(FieldOfMatrixGroup(G), DimensionOfMatrixGroup(G));
    normedPoints := NormedRowVectors(points);
    
    list := GetPartialBaseSGS([], [], gens, 
                    ProjectiveAction, ProjectiveIsIdentity, points, 
                    Identity(G));
    
    generators := list[2];
    base := list[1];
    
    ssInfo := rec(
                  partialSGS := generators,
                  partialBase := base,
                  action := [],
                  schreierTrees := [],
                  oldSGS := [],
                  points := [],
                  newBasePoint := [],
                  hash := [],
                  IsIdentity := []
                  );
    
    for level in [1 .. Length(base)] do
        Add(ssInfo.action, ProjectiveAction);
        Add(ssInfo.points, normedPoints);
        Add(ssInfo.hash, SparseIntKey(points, base[1]));
        Assert(1, not IsBool(ssInfo.hash[level]));
        Add(ssInfo.schreierTrees, SparseHashTable(ssInfo.hash[level]));
        Add(ssInfo.oldSGS, []);
        Add(ssInfo.newBasePoint, ProjectiveNewBasePoint);
        Add(ssInfo.IsIdentity, ProjectiveIsIdentity);
    od;
    
    ssInfo.action[1] := PointAction;
    ssInfo.points[1] := points;
    ssInfo.newBasePoint[1] := NewBasePoint;
    ssInfo.hash[1] :=  SparseIntKey(points, BasisVectors(CanonicalBasis(points))[1]);
    Assert(1, not IsBool(ssInfo.hash[1]));
    ssInfo.schreierTrees[1] := SparseHashTable(ssInfo.hash[1]);
    ssInfo.IsIdentity[1] := IsIdentity;
      
    DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    DebugPrint(3, ["Partial base : ", base]);
    DebugPrint(3, ["Partial sgs : ", generators]);
    
    for level in Reversed([1 .. Length(base)]) do
        SchreierSims(ssInfo, level, Identity(G));
    od;
    
    return [ssInfo.partialBase, ssInfo.partialSGS, ssInfo.schreierTrees];
end);



InstallGlobalFunction(MatrixGroupOrder, function(G)
    local ret, schreierTrees, tree, order;
    
    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
        
    ret := MatrixSchreierSims(G);
    
    order := 1;
    for tree in ret[3] do
        order := order * Size(tree);
    od;
    
    return order;
end);

MATRIXSS_GetTestGroups := 
  function(maxFieldSize, maxMatrixSize) 
    local degree, power, primeNr, prime, groups, groupTypes, type;
    
    groupTypes := 
      Immutable([GeneralLinearGroup, SpecialLinearGroup, 
              GeneralOrthogonalGroup, SpecialOrthogonalGroup,
              GeneralUnitaryGroup, SpecialUnitaryGroup]);
    
    groups := [];
    
    degree := 2;
    while degree <= maxMatrixSize do
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
        
        degree := degree + 1;
    od;
    
    return groups;
end;

MATRIXSS_GetBenchmarkGroups := 
  function(maxClassicalGroupFieldSize, maxClassicalGroupDegree,
          maxReeGroupSize, maxSuzukiSize)
  local groups, size;
    
    groups := MATRIXSS_GetTestGroups(maxClassicalGroupFieldSize,
                      maxClassicalGroupDegree);
    
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
    
    groups := MATRIXSS_GetTestGroups(maxDegree, maxFieldSize);
    
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
    
    test_time := Runtime();
    groups := MATRIXSS_GetBenchmarkGroups(maxDegree, maxFieldSize, 
                      maxReeSize, maxSuzukiSize);
    
    Print("Group\t\tTime [ms]\n\n");
    #Print("-----------------\n");
    
    for group in groups do
        #Print("Checking group : ", group, "\n");
        
        group_time := Runtime();
        MatrixSchreierSims(group);
        #Print("Time for this group : ", Runtime() - group_time, "\n");
        Print(group, "\t", Runtime() - group_time, "\n");
    od;

    Print("Total time for test : ", Runtime() - test_time, "\n");
    
    Print("Benchmark completed\n");
    return true;
end);

#E