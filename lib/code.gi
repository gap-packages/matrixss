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

MATRIXSS_DEBUGLEVEL := 2;

MATRIXSS_DebugPrint := function(level, message)
    #Info(MatrixSchreierSimsInfo, level, message);
    #Info(MatrixSchreierSimsInfo, level, message, "\n");
    if level <= MATRIXSS_DEBUGLEVEL then
        CallFuncList(Print, Concatenation(message, ["\n"]));
    fi;
end;

# The action of a group element g (a matrix) on a point p (a row vector)
# The action is from the right
# Used in Matrix Schreier-Sims
MATRIXSS_MSSAction := function(element, point)
    
    if not IsMatrix(element) or not IsRowVector(point) then
        MATRIXSS_DebugPrint(1, ["Element : ", element]);
        MATRIXSS_DebugPrint(1, ["Point : ", point]);
        Error("<element> must be a matrix and <point> must be a row vector");
    fi;
    
    return point * element;
end;

# return all points (as a list) in the orbit of the point 
# which is root of the schreier tree
MATRIXSS_GetOrbit := function(schreierTree)
    return HashKeyEnumerator(schreierTree);
end;

# Check if a given point is in the orbit defined by the given Schreier tree
MATRIXSS_IsPointInOrbit := function(schreierTree, point)
    if not IsBool(GetHashEntry(schreierTree, point)) then
        return true;
    else
        return false;
    fi;    
end;

# Get the label of the edge originating at the given point, and directed 
# towards the root of the given Schreier tree
MATRIXSS_GetSchreierTreeEdge := function(schreierTree, point)
    return GetHashEntry(schreierTree, point);
end;

# Get size of orbit defined by the given Schreier tree
MATRIXSS_GetOrbitSize := function(schreierTree)
    return Size(schreierTree);
end;

# Helper function for SchreierTree.
MATRIXSS_ComputeSchreierTree := function(generators, tree, root, action)
    local newElement, childElements, generator, child, elements, element;
    
    elements := [root];
    repeat
        childElements := [];
        for element in elements do
            for generator in generators do
                newElement := action(generator[1], element);
                
                # check that the element is not already in the Schreier tree
                # (do not create cycles)
                if not MATRIXSS_IsPointInOrbit(tree, newElement) then
                    AddHashEntry(tree, newElement, generator);
                    Add(childElements, newElement);
                fi;
            od;    
        od;
        elements := childElements;
    until IsEmpty(elements);
    
    return tree;
end;      

MATRIXSS_ExtendSchreierTree := 
  function(oldTree, root, generators, oldGenerators, action)
    local tree, point, generator, newPoint, points, newPoints;
    
    tree := ShallowCopy(oldTree);
    points := MATRIXSS_GetOrbit(tree);
    repeat
        newPoints := [];
        for point in points do
            for generator in generators do
                if not MATRIXSS_IsPointInOrbit(oldTree, point) or
                   not generator in oldGenerators then
                    newPoint := action(generator[1], point);
                    if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                        AddHashEntry(tree, newPoint, generator);
                        Add(newPoints, newPoint);
                    fi;
                fi;
            od;
        od;
        points := newPoints;
    until IsEmpty(points);
    
    return tree;
end;    

# Computes a Schreier tree 
# root - root of the Schreier tree
# generators - generators for group
# points - point set where root comes from
# action - the action used to create the tree
#
# This is just a computation of a spanning tree for a connected component
MATRIXSS_SchreierTree := function(generators, points, root, action, identity)
    local tree;
    
    # Create Schreier vector
    tree := SparseHashTable(SparseIntKey(points, root));
    
    # Make the root point to itself 
    AddHashEntry(tree, root, [identity, identity]);
    
    # Fill Schreier vector
    return MATRIXSS_ComputeSchreierTree(generators, tree, root, action);
end;

# Compute the group element that connects the root of the Schreier tree to
# a given point
# this function assumes that the point actually is in the orbit described by
# the given Schreier tree
MATRIXSS_OrbitElement := 
  function(schreierTree, point, action, identity, onlyInverse)
    local element, edge;
    
    # the group element and its inverse
    element := [identity, identity];
    
    repeat
        edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
        
        Assert(1, not IsBool(edge), "Point not in orbit!\n");
        
        if edge[1] = identity then
            return element;
        fi;
        
        point := action(edge[2], point);
        if not onlyInverse then
            element[1] := edge[1] * element[1];
        fi;
        
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
MATRIXSS_Membership := 
  function(schreierTrees, base, element, action, identity)
    local level, residue, word, point;
    
    residue := element;
    for level in [1 .. Length(base)] do
        MATRIXSS_DebugPrint(5, ["residue: ", residue[1], "\nbase: ", 
                base[level], "\naction", action]);
        point := action(residue[1], base[level]);
        
        if not MATRIXSS_IsPointInOrbit(schreierTrees[level], point) then
            return [residue, level];
        fi;
        
        word := MATRIXSS_OrbitElement(schreierTrees[level], point, action, 
                        identity, false);
        residue[1] := residue[1] * word[2];
        residue[2] := word[1] * residue[2];
    od;
    
    return [residue, Length(base) + 1];
end; 

# find a point not in base that is moved by element
# (element fixes the base)
MATRIXSS_NewBasePoint := function(base, element, action, identity, field)
    local basis, point, basePoint;
    
    basePoint := fail;
    basis := BasisVectors(CanonicalBasis(field));
    
    for point in basis do
        if not action(element[1], point) = point and not point in base then
            basePoint := point;
            break;
        fi;
    od;
    
    if IsBool(basePoint) then
        Error("No new base point found!\n");
    fi;
    
    return basePoint;
end;

# Create a Schreier generator for the stabiliser in the group which has 
# "generator" as one of its generators. The stabiliser fixes "point".
MATRIXSS_GetSchreierGenerator := 
  function(schreierTree, generator, point, action, identity)
    local element1, element2, edge, inv_edge;
    
    element1 := MATRIXSS_OrbitElement(schreierTree, point, action,
                        identity, false);
    element2 := MATRIXSS_OrbitElement(schreierTree, 
                        action(generator[1], point), action, identity, false);
            
    edge := element1[1] * generator[1] * element2[2];
    inv_edge := element2[1] * generator[2] * element1[2];
    
    return [edge, inv_edge];
end;

MATRIXSS_SchreierSims := function(ssInfo, level, identity, field)
    local generator, point, orbit, strip, schreierGenerator, element, action, 
          recursiveLevel, schreierTree, SGS, oldSGS, newSGS, points, newPoint, i;
    
    action := ssInfo.action[level];
    
    newSGS := [];
    for element in ssInfo.partialSGS do
        for point in ssInfo.partialBase{[1 .. level - 1]} do
            if not action(element[1], point) = point then
                AddSet(newSGS, element);
                break;
            fi;
        od;
    od;
#    MATRIXSS_DebugPrint(4, ["Bad SGS for level ", level, " is ", newSGS]);

    SGS := Immutable(Difference(ssInfo.partialSGS, newSGS));
#    oldSGS := ssInfo.fixingGens[level];
    
#    SGS := Immutable(ssInfo.fixingGens[level]);
#    oldSGS := Intersection(SGS, oldGens);
    
#    if IsEmpty(ssInfo.schreierTrees[level]) then
        ssInfo.schreierTrees[level] := 
          MATRIXSS_SchreierTree(SGS, field, ssInfo.partialBase[level], 
                  ssInfo.action[level], identity);
#        ssInfo.oldSchreierTrees[level] := ssInfo.schreierTrees[level];
#    else
        #MATRIXSS_DebugPrint(4, ["SGS for level ", level, " is ", SGS]);
#        ssInfo.oldSchreierTrees[level] := ssInfo.schreierTrees[level];
#        ssInfo.schreierTrees[level] := 
#          MATRIXSS_ExtendSchreierTree(ssInfo.oldSchreierTrees[level],
#                  ssInfo.partialBase[level], SGS, oldSGS, action);
#    fi;
    
    #ssInfo.fixingGens[level] := SGS;
    
    #ssInfo.schreierTrees[level] := 
    #  MATRIXSS_SchreierTree(SGS, field, ssInfo.partialBase[level], 
    #          ssInfo.action[level], identity);

    orbit := MATRIXSS_GetOrbit(ssInfo.schreierTrees[level]);
    MakeImmutable(orbit);
    
    MATRIXSS_DebugPrint(4, ["Orbit size for level ", level, " is ", 
            MATRIXSS_GetOrbitSize(ssInfo.schreierTrees[level])]);
    
    for point in orbit do
        #if MATRIXSS_IsPointInOrbit(ssInfo.oldSchreierTrees[level], point) then
        #    newSGS := newGens;
        #else
        #    newSGS := SGS;
        #fi;
        
        for generator in SGS do
            #if not MATRIXSS_IsPointInOrbit(ssInfo.oldSchreierTrees[level],
            #           point) or not generator in oldGens then
                schreierGenerator := 
                  MATRIXSS_GetSchreierGenerator(ssInfo.schreierTrees[level],
                          generator, point, ssInfo.action[level], identity);
                
                Assert(3, not schreierGenerator[1] = identity, 
                       "Identity Schreier generator!");
                if schreierGenerator[1] = identity then
                    continue;
                fi;
                
                MATRIXSS_DebugPrint(4, ["Schreier Generator : ", 
                        schreierGenerator]);
                
                points := [level + 1 .. Length(ssInfo.partialBase)];
                strip := MATRIXSS_Membership(ssInfo.schreierTrees{points},
                                 ssInfo.partialBase{points},
                                 schreierGenerator, 
                                 ssInfo.action[level], identity);
                
                # The drop-out level is in range 
                # [1 .. Length(ssInfo.partialBase) - level] 
                # but we want the range given by points
                strip[2] := strip[2] + level;
                
                MATRIXSS_DebugPrint(3, ["Dropout level : ", strip[2]]);
                
                if not strip[1][1] = identity then
                    MATRIXSS_DebugPrint(3, ["Residue found"]);

                    AddSet(ssInfo.partialSGS, Immutable(strip[1]));
                    AddSet(ssInfo.partialSGS, Immutable(Reversed(strip[1])));
                    
                    #for i in [1 .. Minimum(strip[2], 
                    #        Length(ssInfo.partialBase))] do
                    #    AddSet(ssInfo.fixingGens[i], Immutable(strip[1]));
                    #    AddSet(ssInfo.fixingGens[i], 
                    #           Immutable(Reversed(strip[1])));
                    #od;

                    if strip[2] = Length(ssInfo.partialBase) + 1 then
                        MATRIXSS_DebugPrint(3, ["Finding new base point"]);
                        
                        newPoint := 
                          MATRIXSS_NewBasePoint(ssInfo.partialBase, strip[1],
                                  ssInfo.action[level], identity, field);
                        
                        Add(ssInfo.partialBase, newPoint);
                        Add(ssInfo.action, MATRIXSS_MSSAction);
                        #Add(ssInfo.fixingGens, []);
                        #Add(ssInfo.schreierTrees, []);
                            #MATRIXSS_SchreierTree([], 
                            #    field, newPoint, 
                            #    MATRIXSS_MSSAction, identity));
                    fi;
                    
                    
                    for recursiveLevel in Reversed([level + 1 .. strip[2]]) do
                        #newSGS := ShallowCopy(ssInfo.partialSGS);
                        #RemoveSet(newSGS, strip[1]);
                        #RemoveSet(newSGS, Reversed(strip[1]));
                        MATRIXSS_SchreierSims(ssInfo, 
                                recursiveLevel, identity, field);
                    od;
                fi;
            #fi;
        od;
    od;
end;

MATRIXSS_GetPartialBaseSGS := 
  function(base, sgs, generators, action, field, identity)
    local gens, point, nFixedPoints, newPoint, element, list, points, newSGS;
    
    points := ShallowCopy(base);
    gens := UnionSet(generators, sgs);
    RemoveSet(gens, identity);
    newSGS := [];
    
    # we make a partial strong generating set which also contain
    # inverses of all elements
    for element in gens do
        nFixedPoints := 0;
        for point in points do
            if action(element, point) = point then
                nFixedPoints := nFixedPoints + 1;
            else
                break;
            fi;
        od;
        list := [element, Inverse(element)];
        if nFixedPoints = Length(points) then
            newPoint := MATRIXSS_NewBasePoint(points, list, action, 
                                identity, field);
            Add(points, newPoint);
        fi;

        # Save reference to generator and its inverse
        # Then inverses need not be calculated later
        AddSet(newSGS, Immutable(list));
        if not element * element = identity then
            AddSet(newSGS, Immutable(Reversed(list)));
        fi;
    od;
    
    return [points, newSGS];
end;

# An implementation of the Schreier-Sims algorithm, for matrix groups
InstallGlobalFunction(MatrixSchreierSims, function(G)
    local ssInfo, generators, base, trees, points, element, gens, list, level;

    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    gens := GeneratorsOfGroup(G);
    points := FullRowSpace(FieldOfMatrixGroup(G), DimensionOfMatrixGroup(G));
    base := BasisVectors(CanonicalBasis(points));
    
    list := MATRIXSS_GetPartialBaseSGS([], [], gens, MATRIXSS_MSSAction, 
                    points, Identity(G));
    
    generators := list[2];
    base := list[1];

    ssInfo := rec(
                  partialSGS := generators,
                  partialBase := base,
                  action := [],
                  schreierTrees := [],
                  oldSchreierTrees := [],
                  fixingGens := []
                  );
    
    #ssInfo.fixingGens[1] := ssInfo.partialSGS;
    
    for level in [1 .. Length(base)] do
        Add(ssInfo.action, MATRIXSS_MSSAction);
        
        #Add(ssInfo.schreierTrees, []);
        #Add(ssInfo.schreierTrees, 
        #  MATRIXSS_SchreierTree(ssInfo.fixingGens[level], points, 
        #          ssInfo.partialBase[level], ssInfo.action[level], 
        #          Identity(G)));
        
        #ssInfo.fixingGens[level + 1] := [];
        #for element in ssInfo.fixingGens[level] do
        #    if ssInfo.action[level](element[1], ssInfo.partialBase[level]) =
        #       ssInfo.partialBase[level] then
        #        AddSet(ssInfo.fixingGens[level + 1], element);
        #    fi;
        #od;
        #MakeImmutable(ssInfo.fixingGens[level + 1]);                
    od;
    
    MATRIXSS_DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    MATRIXSS_DebugPrint(3, ["Partial base : ", base]);
    MATRIXSS_DebugPrint(3, ["Partial sgs : ", generators]);
    
    for level in Reversed([1 .. Length(base)]) do
        MATRIXSS_SchreierSims(ssInfo, level, Identity(G), points);
    od;
    
    return [ssInfo.partialBase, ssInfo.partialSGS, ssInfo.schreierTrees];
end);
    
InstallGlobalFunction(MatrixGroupOrder, function(G)
    local ret, schreierTrees, tree, order;
    
    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
        
    MATRIXSS_DebugPrint(2, ["Entering MatrixGroupOrder"]);
    MATRIXSS_DebugPrint(2, ["Input group : ", G]);
    
    ret := MatrixSchreierSims(G);
    
    MATRIXSS_DebugPrint(3, ["Base : ", ret[1]]);
    MATRIXSS_DebugPrint(3, ["SGS : ", ret[2]]);
    
    order := 1;
    for tree in ret[3] do
        MATRIXSS_DebugPrint(3, ["Orbit size : ", MATRIXSS_GetOrbitSize(tree)]);
        order := order * MATRIXSS_GetOrbitSize(tree);
    od;
    
    MATRIXSS_DebugPrint(2, ["Order : ", order]);

    return order;
end);

MATRIXSS_GetTestGroups := 
  function(maxFieldSize, maxMatrixSize) 
    local degree, power, primeNr, prime, groups, groupTypes, type;
    
    groupTypes := 
      Immutable([GeneralLinearGroup, SpecialLinearGroup, GeneralUnitaryGroup,
              SpecialUnitaryGroup, GeneralOrthogonalGroup, 
              SpecialOrthogonalGroup]);
                          
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
            Print("Correct order: ", size1);
            Print("Computed order: ", size2);
            Error("Order difference!");
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
    
    for group in groups do
        Print("Checking group : ", group, "\n");
        
        group_time := Runtime();
        size := MatrixGroupOrder(group);
        Print("Time for this group : ", Runtime() - group_time, "\n");        
    od;
    
    Print("Total time for test : ", Runtime() - test_time, "\n");
    
    Print("Benchmark completed\n");
    return true;
end);

#E