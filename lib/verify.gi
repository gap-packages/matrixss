###############################################################################
##
#W    verify.gi  The Matrix Schreier Sims package 
#W               The Verify routine by Sims.
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-08-14 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/lib/verify_gi") := 
  "@(#)$Id$";

MATRIXSS_BaseChange := function(ssInfo, partialSGS, level, identity)
    local SGS, point, orbitSize, generator, orbit, schreierGen, inverseGen,
          recursiveLevel;
    
    MATRIXSS_DebugPrint(2, ["Switching places of base points at level ",
            level, " and ", level + 1]);
    
    if level > 1 then
        SGS := ShallowCopy(ssInfo[level - 1].partialSGS);
    else
        SGS := ShallowCopy(partialSGS);
    fi;
    MakeImmutable(SGS);
    
    MATRIXSS_DebugPrint(4, ["Orbit size 1 : ", 
            MATRIXSS_GetOrbitSize(ssInfo[level].schreierTree.Tree)]);
    MATRIXSS_DebugPrint(4, ["Orbit size 2 : ", 
            MATRIXSS_GetOrbitSize(ssInfo[level + 1].schreierTree.Tree)]);
    
    orbitSize := MATRIXSS_GetOrbitSize(ssInfo[level].schreierTree.Tree) * 
                 MATRIXSS_GetOrbitSize(ssInfo[level + 1].schreierTree.Tree);
    
    point := ssInfo[level];
    ssInfo[level] := ssInfo[level + 1];
    ssInfo[level + 1] := point;
    
    MATRIXSS_DebugPrint(4, ["Schreier structure : ", ssInfo]);
    MATRIXSS_DebugPrint(4, ["Level : ", level]);
    
    ssInfo[level].schreierTree := 
      MATRIXSS_GetSchreierTree(fail, ssInfo[level].partialBase, SGS, [], 
              ssInfo[level].action, ssInfo[level].hash, identity);
    orbitSize := orbitSize / 
                 MATRIXSS_GetOrbitSize(ssInfo[level].schreierTree.Tree);
    MATRIXSS_DebugPrint(3, ["Wanted orbit size : ", orbitSize]);
    Assert(1, IsInt(orbitSize));
    
    orbit := MATRIXSS_GetOrbit(ssInfo[level].schreierTree.Tree);
    ssInfo[level + 1].schreierTree := 
      MATRIXSS_CreateInitialSchreierTree(ssInfo[level + 1].partialBase,
              ssInfo[level + 1].hash, identity);
    for point in orbit do
        for generator in SGS do
            schreierGen := MATRIXSS_GetSchreierGenerator(
                                   ssInfo[level].schreierTree.Tree, generator, 
                                   point, ssInfo[level].action, identity);
            if not MATRIXSS_IsPointInOrbit(
                       ssInfo[level + 1].schreierTree.Tree,
                       ssInfo[level + 1].action(ssInfo[level + 1].partialBase,
                               schreierGen[1])) then
                inverseGen := Immutable(Reversed(schreierGen));
                        
                AddSet(partialSGS, schreierGen);
                AddSet(partialSGS, inverseGen);
                        
                for recursiveLevel in [1 .. level + 1] do
                    if ssInfo[recursiveLevel].action(
                               ssInfo[recursiveLevel].partialBase,
                               schreierGen[1]) = 
                       ssInfo[recursiveLevel].partialBase then
                        AddSet(ssInfo[recursiveLevel].partialSGS, schreierGen);
                        AddSet(ssInfo[recursiveLevel].partialSGS, inverseGen);
                    else
                        break;
                    fi;
                od;
                
                # recompute SGS ?
                
                ssInfo[level + 1].schreierTree := 
                  MATRIXSS_GetSchreierTree(ssInfo[level + 1].schreierTree, 
                          ssInfo[level + 1].partialBase, 
                          ssInfo[level].partialSGS, 
                          ssInfo[level + 1].oldSGS, 
                          ssInfo[level + 1].action, 
                          ssInfo[level + 1].hash, identity);
                if MATRIXSS_GetOrbitSize(ssInfo[level + 1].schreierTree.Tree) =
                   orbitSize then
                    break;
                fi;
            fi;
        od;
        if MATRIXSS_GetOrbitSize(ssInfo[level + 1].schreierTree.Tree) =
           orbitSize then
            break;
        fi;
    od;
    MATRIXSS_DebugPrint(4, ["Actual orbit size : ", 
            MATRIXSS_GetOrbitSize(ssInfo[level + 1].schreierTree.Tree)]);
    Assert(1, MATRIXSS_GetOrbitSize(ssInfo[level + 1].schreierTree.Tree) =
           orbitSize);
end;

MATRIXSS_StabiliserGens := function(ssInfo, partialSGS, point, action, 
                                   dictinfo, identity)
    local ssInfoNew, generators, preBase, postBase, level, generator, SGS,
          ssStruct, i;
    
    MATRIXSS_DebugPrint(2, ["Fetching generators for stabiliser of ", point]);
    
    ssInfoNew := ShallowCopy(ssInfo);
    generators := ShallowCopy(partialSGS);
    level := 0;
    for level in [1 .. Length(ssInfoNew)] do
        ssInfoNew[level] := ShallowCopy(ssInfo[level]);
        #ssInfoNew[level].partialSGS := ShallowCopy(ssInfo[level].partialSGS);
        #generators := Difference(generators, ssInfoNew[level].partialSGS);
        if ForAll(generators, generator -> 
                  action(point, generator[1]) = point) then 
            level := level - 1;
            break;
        fi;
        generators := ssInfoNew[level].partialSGS;
    od;
    
    if level = 0 then
        MATRIXSS_DebugPrint(2, ["Whole group stabilises point"]);
        return Immutable(generators);
    fi;
    #Assert(1, level >= 0);
    #ssInfoNew[level + 1].partialSGS := 
    #  ShallowCopy(ssInfo[level + 1].partialSGS);
    preBase := ssInfoNew{[1 .. level]};
    postBase := ssInfoNew{[level + 1 .. Length(ssInfoNew)]};
    MATRIXSS_AugmentBase(preBase, point, action, dictinfo, identity : 
            AlternatingActions := fail);
    ssInfoNew := Concatenation(preBase, postBase);
    
    for i in [1 .. Length(ssInfoNew)] do
        ssInfoNew[i] := ShallowCopy(ssInfoNew[i]);
    od;
    
    generators := ShallowCopy(partialSGS);
    for ssStruct in ssInfoNew do
        ssStruct.schreierTree := 
          MATRIXSS_GetSchreierTree(fail, ssStruct.partialBase, 
                  generators, [], ssStruct.action, ssStruct.hash, identity);
        ssStruct.partialSGS := [];
        Perform(ShallowCopy(generators), function(generator) 
            if ssStruct.action(ssStruct.partialBase, 
                       generator[1]) = ssStruct.partialBase then 
                AddSet(ssStruct.partialSGS, generator);
            else
                RemoveSet(generators, generator);
            fi;
        end);
    od;
    
    generators := ShallowCopy(partialSGS);    
    while level > 0 do
        MATRIXSS_BaseChange(ssInfoNew, generators, level, identity);
        level := level - 1;
    od;
    
    SGS := [];
    for generator in generators do
        if ssInfoNew[1].action(ssInfoNew[1].partialBase, generator[1]) =
           ssInfoNew[1].partialBase then
            AddSet(SGS, generator);
        fi;
    od;
    #Assert(1, SGS = ssInfoNew[1].partialSGS);
    
    return Immutable(SGS);
end;    
    
MATRIXSS_DecomposeOrbit := function(schreierTree, root, generators, action,
                                   hash, identity)
    local orbits, orbit, points, point, tree, index;
    
    MATRIXSS_DebugPrint(2, ["Decomposing orbit based at ", root]);
    
    orbits := [];
    orbit := MATRIXSS_GetOrbit(schreierTree.Tree);
    point := root;
    index := 0;
    MATRIXSS_DebugPrint(3, ["Original orbit size : ", Size(orbit)]);
    repeat
        tree := MATRIXSS_GetSchreierTree(fail, point, generators, generators, 
                        action, hash, identity);
        Add(orbits, Immutable(rec(SchreierTree := tree, Point := point)));
        Perform(orbits, function(i) MATRIXSS_DebugPrint(4, ["Orbit size : ", 
                Size(i.SchreierTree.Tree)]); end);
        MATRIXSS_DebugPrint(3, ["Total size of decomposed orbits : ",
                Sum(List(orbits, i -> Size(i.SchreierTree.Tree)))]);
        if Sum(List(orbits, i -> Size(i.SchreierTree.Tree))) = 
           Size(schreierTree.Tree) then
            break;
        fi;
        repeat
            index := index + 1;
        until ForAll(orbits, i -> not MATRIXSS_IsPointInOrbit(
                      i.SchreierTree.Tree, orbit[index]));
        point := orbit[index];
    until false;
    
    return orbits;
end;

MATRIXSS_VerifySingleGenerator := 
  function(generators, schreierTree, point, action, hash, subGenerator, 
          ssInfo, SGS, identity)
  local subGens, orbits, representatives, stabPoint, orbit, stabGens,
        subRepresentatives, subOrbits, subOrbit, stabRepresentatives, level,
        orbitStabGens, generator, representative, residue, subLevel;
    
    subGens := Difference(generators, AsSet([subGenerator]));
    
    MATRIXSS_DebugPrint(5,  ["Decomposing ", schreierTree, " with root at ", 
            point, " with respect to ", subGens]);
    orbits := MATRIXSS_DecomposeOrbit(schreierTree, point, subGens,
                      action, hash, identity);
    representatives := [];
    for orbit in orbits do
        Add(representatives, MATRIXSS_OrbitElement(schreierTree.Tree, 
                orbit.Point, action, identity));
    od;
    
    stabPoint := action(point, subGenerator[2]);
    stabGens := MATRIXSS_StabiliserGens(ssInfo, SGS, stabPoint, action, 
                        hash, identity);
    subOrbits := [];
    subRepresentatives := [];
    stabRepresentatives := [];
    for level in [1 .. Length(orbits)] do
        
        MATRIXSS_DebugPrint(3,  ["Decomposing ", orbits[level].SchreierTree, 
                " with root at ", orbits[level].Point, " with respect to ", 
                stabGens]);
        
        Add(subOrbits, MATRIXSS_DecomposeOrbit(orbits[level].SchreierTree, 
                orbits[level].Point, stabGens, action, hash, identity));
        Add(subRepresentatives, []);
        Add(stabRepresentatives, []);
        for subOrbit in subOrbits[Length(subOrbits)] do
            residue := MATRIXSS_OrbitElement(orbits[level].SchreierTree.Tree, 
                               subOrbit.Point, action, identity);
            Add(subRepresentatives[Length(subRepresentatives)],
                Immutable([representatives[level][1] * residue[1],
                        residue[2] * representatives[level][2]]));
            Add(stabRepresentatives[Length(stabRepresentatives)],
                MATRIXSS_OrbitElement(schreierTree.Tree, 
                        action(subOrbit.Point, subGenerator[1]), action, 
                        identity));
        od;
    od;
    
    # check 1, 2 and 3
    # check membership in H = <generators, subGenerator> using ssInfo
    
    for level in [1 .. Length(orbits)] do
        orbitStabGens := MATRIXSS_StabiliserGens(ssInfo, SGS, 
                                 orbits[level].Point, action, hash, identity);
        for generator in orbitStabGens do
            residue := MATRIXSS_Membership(ssInfo, 
                               Immutable([representatives[level][1] * 
                                       generator[1] * 
                                       representatives[level][2],
                                       representatives[level][1] * 
                                       generator[2] * 
                                       representatives[level][2]]),
                               identity);
            if residue[1][1] <> identity then
                MATRIXSS_DebugPrint(2, ["Cond1, Residue found : ", 
                        residue[1]]);
                return Immutable(residue[1]);
            fi;
        od;
        
        for subLevel in [1 .. Length(subOrbits[level])] do
            residue := 
              MATRIXSS_Membership(ssInfo, Immutable(
                      [subRepresentatives[level][subLevel][1] *
                       subGenerator[1] *
                       stabRepresentatives[level][subLevel][2],
                       stabRepresentatives[level][subLevel][2] *
                       subGenerator[2] *
                       stabRepresentatives[level][subLevel][1]]),
                      identity);
            if residue[1][1] <> identity then
                MATRIXSS_DebugPrint(2, ["Cond2, Residue found : ", 
                        residue[1]]);
                return Immutable(residue[1]);
            fi;
        od;
    od;
    
    for generator in stabGens do
        residue := MATRIXSS_Membership(ssInfo, 
                           Immutable([subGenerator[2] * generator[1] *
                                   subGenerator[1],
                                   subGenerator[2] * generator[2] *
                                   subGenerator[1]]), identity);
        if residue[1][1] <> identity then
            MATRIXSS_DebugPrint(2, ["Cond3, Residue found : ", 
                    residue[1]]);
            return Immutable(residue[1]);
        fi;
    od;
    
    return Immutable([identity, identity]);
end;

MATRIXSS_IsBlockOfImprimitivity := 
  function(schreierTree, generators, block, action, identity)
    local orbit, partition, image, point, testPoint, representative,
          generator, set, index;
    
    MATRIXSS_DebugPrint(2, ["Checking block of imprimitivity"]);

    orbit := AsSet(ShallowCopy(MATRIXSS_GetOrbit(schreierTree.Tree)));
    orbit := Difference(orbit, block);
    partition := [rec(Block := block, 
                      Element := Immutable([identity, identity]))];
    repeat
        MATRIXSS_DebugPrint(5, ["Orbit size : ", Length(orbit)]);
        testPoint := orbit[1];
        representative := MATRIXSS_OrbitElement(schreierTree.Tree, testPoint,
                                  action, identity);
        image := [];
        for point in block do
            testPoint := action(point, representative[1]);
            if not testPoint in orbit then
                for image in partition do
                    if testPoint in image.Block then
                        MATRIXSS_DebugPrint(3, ["Not a block"]);
                        return Immutable([representative[1] * image.Element[2],
                                       image.Element[1] * representative[2]]);
                    fi;
                od;
            fi;
            Add(image, testPoint);
        od;
        image := AsSet(image);
        MATRIXSS_DebugPrint(5, ["Image size : ", Length(image)]);
        orbit := Difference(orbit, image);
        Add(partition, rec(Block := image, 
                                    Element := Immutable(representative)));
    until IsEmpty(orbit);
    
    for block in partition do
        for generator in generators do
            image := [];
            for point in block.Block do
                Add(image, action(point, generator[1]));
            od;
            image := rec(Block := AsSet(image), 
                         Element := Immutable([block.Element[1] * generator[1],
                                 generator[2] * block.Element[2]]));
            for set in partition do
                if image.Block <> set.Block and 
                   not IsEmpty(Intersection(image.Block, set.Block)) then
                    MATRIXSS_DebugPrint(3, ["Gens not ok"]);
                    return Immutable([image.Element[1] * set.Element[2],
                                   set.Element[1] * image.Element[2]]);
                fi;
            od;
        od;
    od;
    
    return Immutable([identity, identity]);
end;
    
MATRIXSS_VerifyMultipleGenerators :=   
  function(generators, schreierTree, point, action, hash, subGenerators, 
          ssInfo, SGS, identity, points, IsIdentity, field)
  local lastSchreierTree, level, residue, gens, ssInfoNew, orbit,
        dictinfo, element, newSchreierTree, testPoint, i;
    
    gens := Difference(generators, subGenerators);
    gens := Union(gens, AsSet([subGenerators[1]]));
    MATRIXSS_DebugPrint(3, ["Verifying single generator"]);
    
    schreierTree := MATRIXSS_GetSchreierTree(fail, point, gens, [], action, 
                            hash, identity);
    residue := MATRIXSS_VerifySingleGenerator(gens, schreierTree, 
                       point, action, hash, subGenerators[1], ssInfo, SGS, 
                       identity);
    if residue[1] <> identity then
        return residue;
    fi;
    
    ssInfoNew := [rec(
                      partialSGS := SGS, 
                      partialBase := point,
                      action := action,
                      points := points,
                      hash := hash,
                      schreierTree := schreierTree,
                      oldSGS := AsSSortedList([]),
                      relations := [],
                      IsIdentity := IsIdentity)];
    Append(ssInfoNew, ssInfo);
    
    for level in [2 .. Length(subGenerators)] do
        
        MATRIXSS_DebugPrint(3, ["Verifying ", level, " generators"]);
        residue := MATRIXSS_Membership(ssInfoNew, subGenerators[level], 
                           identity);
        if residue[1][1] <> identity then
            MATRIXSS_DebugPrint(3, ["Residue found"]);
            
            if action(point, residue[1][1]) = point then
                return residue[1];
            fi;
            
            MATRIXSS_DebugPrint(3, ["Residue does not fix base point"]);
            
            lastSchreierTree := schreierTree;
            gens := Union(gens, AsSet([subGenerators[level]]));
            schreierTree := MATRIXSS_GetSchreierTree(fail, point, gens, [], 
                                    action, hash, identity);
            
            # check if orbit of lastSchreierTree is a block of imprimitivity
            orbit := 
              AsSet(ShallowCopy(MATRIXSS_GetOrbit(lastSchreierTree.Tree)));
            residue := MATRIXSS_IsBlockOfImprimitivity(
                               schreierTree, gens, 
                               orbit, action, identity);
            if residue[1] <> identity then
                MATRIXSS_DebugPrint(2, ["Residue from primitivity check : ",
                        residue]);
                i := 1;
                repeat
                    testPoint := orbit[i];
                    if MATRIXSS_IsPointInOrbit(lastSchreierTree.Tree,
                               action(testPoint, residue[1])) then
                        break;
                    fi;
                    i := i + 1;
                until false;
                return MATRIXSS_GetSchreierGenerator(lastSchreierTree.Tree,
                               residue, testPoint, action, identity);
            fi;
            
            # Ugly hack to make Dictionaries work correctly
            # Make sure that the hash function uses the correct finite field,
            # and not a smaller one. It happened earlier that the hash function
            # chose the prime field instead of the whole field.
            element := ShallowCopy(identity[1]);
            element[1] := PrimitiveRoot(field);
            dictinfo := [[element], true, [[element]]];
            
            MATRIXSS_DebugPrint(2, ["Schreier tree for point ", orbit]);
            MATRIXSS_DebugPrint(2, ["Using dictinfo ", dictinfo]);
            newSchreierTree := MATRIXSS_GetSchreierTree(fail, orbit, 
                                       gens, [], OnSets, dictinfo, identity);
            
            MATRIXSS_DebugPrint(2, ["Verifying single gen"]);
            residue := MATRIXSS_VerifySingleGenerator(gens, newSchreierTree,
                               orbit, OnSets, dictinfo, subGenerators[level],
                               ssInfoNew, 
                               ShallowCopy(Difference(gens, 
                                       AsSet([subGenerators[level]]))), 
                               identity);
            if residue[1] <> identity then
                element := MATRIXSS_OrbitElement(lastSchreierTree.Tree, 
                                   action(point, residue[1]), action, 
                                   identity);
                return Immutable([residue[1] * element[2],
                               element[1] * residue[2]]);
            fi;
            
            ssInfoNew[1].partialSGS := ShallowCopy(gens);
            ssInfoNew[1].schreierTree := schreierTree;
        fi;
    od;
    
    return Immutable([identity, identity]);
end;
                
MATRIXSS_VerifyLevel := function(ssInfo, partialSGS, level, identity)
    local generators, SGS, residue;
    
    if level > 1 then
        SGS := ShallowCopy(ssInfo[level - 1].partialSGS);
    else
        SGS := ShallowCopy(partialSGS);
    fi;
    MakeImmutable(SGS);
    
    MATRIXSS_DebugPrint(2, ["Verify at level ", level]);
    return MATRIXSS_VerifyMultipleGenerators(SGS, 
                   ssInfo[level].schreierTree.Tree, ssInfo[level].partialBase,
                   ssInfo[level].action, ssInfo[level].hash,
                   Difference(SGS, ssInfo[level].partialSGS), 
                   ssInfo{[level + 1 .. Length(ssInfo)]}, 
                   ssInfo[level].partialSGS, identity, 
                   ssInfo[level].points, ssInfo[level].IsIdentity,
                   FieldOfMatrixGroup(Group(List(SGS, i -> i[1]), identity)));
end;

InstallGlobalFunction(MatrixSchreierSimsVerify, function(ssInfo, SGS, identity)
    local level, residue;
    
    for level in Reversed([1 .. Length(ssInfo)]) do
        residue := MATRIXSS_VerifyLevel(ssInfo, SGS, level, identity);
        if residue[1] <> identity then
            MATRIXSS_DebugPrint(2, ["Residue found at level ", level]);
            return rec(Residue := residue, Level := level);
        fi;
    od;
    
    
    return rec(Residue := identity, Level := 0);
end);

###############################################################################
#E
