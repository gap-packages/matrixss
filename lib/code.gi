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

# return generating set of a point stabiliser
# generators - generating set for group to compute in
# root - point to compute stabiliser of
# schreierTree - should have root as root
# action - the action to use
# identity - the identity element in the group
MATRIXSS_Stabiliser := 
  function(generators, root, schreierTree, action, identity)
    local stabiliser, orbit, element1, element2, 
          generator, point, edge, inv_edge, list;
    
    stabiliser := [];
    orbit := MATRIXSS_GetOrbit(schreierTree);
    
    # Compute Schreier generators of the stabiliser
    for generator in generators do
        for point in orbit do
            if not MATRIXSS_GetSchreierTreeEdge(schreierTree, point) = 
               generator[1] then
                element1 := MATRIXSS_OrbitElement(schreierTree, point, 
                                    action, identity, false);
                element2 := MATRIXSS_OrbitElement(schreierTree, 
                                    action(generator[1], point), 
                                    action, identity, false);
                
                edge := element1[1] * generator[1] * element2[2];
                inv_edge := element2[1] * generator[2] * element1[2];
                
                # Important to make sure no generator is the identity
                if not edge = identity then
                    list := [edge, inv_edge];
                    AddSet(stabiliser, Immutable(list));
                    #AddSet(stabiliser, Immutable(Reversed(list)));
                fi;
            fi;
        od;
    od;
    
    return stabiliser;
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
        MATRIXSS_DebugPrint(3, ["residue: ", residue[1], "\nbase: ", base[level],
                "\naction", action]);
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

# ssInfo := rec(
#           schreierTrees[1..k]
#           partialSGS[1..n]
#           partialBase[1..k]
#           action[1..k]


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

MATRIXSS_SchreierSims1 := function(ssInfo, level, identity, field)
    local generator, point, orbit, strip, schreierGenerator, element, action,
          recursiveLevel, schreierTree, SGS, newSGS, points, newPoint;
    
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
    MATRIXSS_DebugPrint(4, ["Bad SGS for level ", level, " is ", newSGS]);

    SGS := Difference(ssInfo.partialSGS, newSGS);
    MakeImmutable(SGS);
    
    MATRIXSS_DebugPrint(4, ["SGS for level ", level, " is ", SGS]);
    #ssInfo.oldSchreierTrees[level] := ssInfo.schreierTrees[level];
    #ssInfo.schreierTrees[level] := 
    #  MATRIXSS_ExtendSchreierTree(ssInfo.oldSchreierTrees[level],
    #          ssInfo.partialBase[level], ssInfo.partialSGS,
    #          ssInfo.oldSGS, ssInfo.action[level]);
    
    ssInfo.schreierTrees[level] := 
      MATRIXSS_SchreierTree(SGS, field, ssInfo.partialBase[level], 
              ssInfo.action[level], identity);

    orbit := MATRIXSS_GetOrbit(ssInfo.schreierTrees[level]);
    MakeImmutable(orbit);
    
    MATRIXSS_DebugPrint(4, ["Orbit size for level ", level, " is ", 
            MATRIXSS_GetOrbitSize(ssInfo.schreierTrees[level])]);
    
    for point in orbit do
        for generator in SGS do
            #if not MATRIXSS_IsPointInOrbit(ssInfo.oldSchreierTrees[level],
            #           point) or
            #   not generator in ssInfo.oldSGS then
                schreierGenerator := 
                  MATRIXSS_GetSchreierGenerator(ssInfo.schreierTrees[level],
                          generator, point, ssInfo.action[level], identity);
                
                MATRIXSS_DebugPrint(4, ["Schreier Generator : ", 
                        schreierGenerator]);
                
                points := [level .. Length(ssInfo.partialBase)];

                strip := 
                  MATRIXSS_Membership(ssInfo.schreierTrees{points},
                          ssInfo.partialBase{points},
                          schreierGenerator, 
                          ssInfo.action[level], identity);
                
                MATRIXSS_DebugPrint(3, ["Dropout level : ", strip[2]]);
                
                if not strip[1][1] = identity then
                    MATRIXSS_DebugPrint(3, ["Residue found"]);
                    #ssInfo.oldSGS := ssInfo.partialSGS;
                    AddSet(ssInfo.partialSGS, Immutable(strip[1]));
                    AddSet(ssInfo.partialSGS, Immutable(Reversed(strip[1])));
                    if strip[2] = Length(ssInfo.partialBase) + 1 then
                        MATRIXSS_DebugPrint(2, ["Finding new base point"]);
                        newPoint := 
                          MATRIXSS_NewBasePoint(ssInfo.partialBase, strip[1],
                                  ssInfo.action[level], identity, field);
                        Add(ssInfo.partialBase, newPoint);
                        Add(ssInfo.action, MATRIXSS_MSSAction);
                        Add(ssInfo.schreierTrees,
                            MATRIXSS_SchreierTree([], field, newPoint, 
                                    MATRIXSS_MSSAction, identity));
                    fi;
                    
                    for recursiveLevel in Reversed([level + 1 .. strip[2]]) do
                        #ssInfo.oldSGS := ShallowCopy(ssInfo.partialSGS);
                        #RemoveSet(ssInfo.partialSGS, strip[1]);
                        #RemoveSet(ssInfo.partialSGS, Reversed(strip[1]));
                        MATRIXSS_SchreierSims1(ssInfo, recursiveLevel,
                                identity, field);
                    od;
                fi;
            #fi;
        od;
    od;
end;

# Compute Schreier trees for the stabiliser chain generated by sgs and
# associated to base.
# Tree j is created with base point j as root, using generators in sgs that fix
# base points 1 .. j-1 (or all generators, in case of j = 1)
MATRIXSS_GetSchreierTrees := function(sgs, base, field, action, identity)
    local generators, point, element, schreierTrees, gens;
    
    # Create mutable generator list
    generators := [];
    UniteSet(generators, sgs);
    
    schreierTrees := [];
    for point in base do
        Add(schreierTrees, 
            MATRIXSS_SchreierTree(generators, field, point, action, identity));
        
        gens := ShallowCopy(generators);
        for element in generators do
            if not action(element[1], point) = point then
                RemoveSet(gens, element);
            fi;
        od;
        generators := gens;
    od;
    
    # If sgs is a strong generating set, then no generator fixes all points
    Assert(3, IsEmpty(generators), "Not a real base!\n");
    
    return schreierTrees;
end;

# An implementation of the Schreier-Sims algorithm, for matrix groups
MatrixSchreierSims1 := function(G)
    local ssInfo, generators, base, trees, points, element, gens, list, level;

    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    gens := GeneratorsOfGroup(G);
    points := FullRowSpace(FieldOfMatrixGroup(G), DimensionOfMatrixGroup(G));
    base := BasisVectors(CanonicalBasis(points));
    
    # we make generators a partial strong generating set which also contain
    # inverses of all elements
    generators := [];
    for element in gens do
        if not element = Identity(G) then
            # Save reference to generator and its inverse
            # Then inverses need not be calculated later
            list := [element, Inverse(element)];
            AddSet(generators, Immutable(list));
            AddSet(generators, Immutable(Reversed(list)));
        fi;
    od;
    
    ssInfo := rec(
                  partialSGS := generators,
                  oldSGS := [],
                  partialBase := base,
                  action := [],
                  schreierTrees := [],
                  oldSchreierTrees := []
                  );
    
    
    for level in [1 .. Length(base)] do
        Add(ssInfo.action, MATRIXSS_MSSAction);
        Add(ssInfo.schreierTrees, 
          MATRIXSS_SchreierTree(ssInfo.oldSGS, points, 
                  ssInfo.partialBase[level], ssInfo.action[level], 
                  Identity(G)));
    od;
    
    MATRIXSS_DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    MATRIXSS_DebugPrint(3, ["Partial base : ", base]);
    MATRIXSS_DebugPrint(3, ["Partial sgs : ", generators]);
    
    for level in Reversed([1 .. Length(base)]) do
        MATRIXSS_SchreierSims1(ssInfo, level, Identity(G), points);
    od;
    
    #ssInfo.schreierTrees := 
    #  MATRIXSS_GetSchreierTrees(ssInfo.partialSGS, ssInfo.partialBase, points, 
    #          MATRIXSS_MSSAction, Identity(G));
    
    return [ssInfo.partialBase, ssInfo.partialSGS, ssInfo.schreierTrees];
end;
    

# Recursive Schreier-Sims over matrix group generated by S
# partialBase - a partial base
# partialSGS - a partial strong generating set
MATRIXSS_SchreierSims := function(partialBase, partialSGS, schreierTrees, 
                                 action, identity, field)
                local newBase, newSGS, newSchreierTrees, ret, element, point, 
                      orbit, generator, 
                      numFixedBasePoints, stabiliser, generators;
    
    if IsEmpty(partialSGS) then
        MATRIXSS_DebugPrint(2, ["Basis case, no generators"]);
        return [partialBase, partialSGS, schreierTrees];
    fi;
    
    repeat        
        MATRIXSS_DebugPrint(3, ["At current level"]);
        MATRIXSS_DebugPrint(3, ["Partial base : ", partialBase]);
        MATRIXSS_DebugPrint(3, ["Partial sgs : ", partialSGS]);
        
        MATRIXSS_DebugPrint(2, ["Getting partial base and sgs for next lower level"]);
        newBase := partialBase{[2 .. Length(partialBase)]};
        newSchreierTrees := schreierTrees{[2 .. Length(schreierTrees)]};
        newSGS  := [];
        for element in partialSGS do
            if action(element[1], partialBase[1]) = partialBase[1] then
                AddSet(newSGS, element);
            fi;
        od;
        
        MATRIXSS_DebugPrint(3, ["Partial base : ", newBase]);
        MATRIXSS_DebugPrint(3, ["Partial sgs : ", newSGS]);
        
        MATRIXSS_DebugPrint(2, ["Updating to complete base and sgs"]);
        ret              := MATRIXSS_SchreierSims(newBase, newSGS, 
                                    newSchreierTrees, 
                                    action, identity, field);
        newBase          := ret[1];
        newSGS           := ret[2];
        newSchreierTrees := ret[3];
        
        MATRIXSS_DebugPrint(3, ["Base : ", newBase]);
        MATRIXSS_DebugPrint(3, ["SGS : ", newSGS]);

        partialBase := partialBase{[1]};
        Append(partialBase, newBase);
        partialSGS := UnionSet(partialSGS, newSGS);
        #MakeImmutable(partialSGS);
        
        # recompute Schreier trees, when base and SGS has changed
        # only recompute first tree? we only need first one in this recursion
        schreierTrees := MATRIXSS_GetSchreierTrees(partialSGS, partialBase, 
                                 field, action, identity);
        
        # get Schreier generators of stabiliser
        stabiliser := MATRIXSS_Stabiliser(partialSGS, partialBase[1],
                              schreierTrees[1], action, identity);

        # check each Schreier generators if they are in the group generated by
        # newSGS
        generators := [];
        for element in stabiliser do                
            if MATRIXSS_Membership(newSchreierTrees, newBase, newSGS, 
                       element, action, identity) then
                AddSet(generators, element);
            else
                break;
            fi;
        od;
        
        if Length(generators) = Length(stabiliser) then
            MATRIXSS_DebugPrint(3, ["No extension needed"]);
            return [partialBase, partialSGS, schreierTrees];
        else
            MATRIXSS_DebugPrint(3, ["Extending SGS with : ", element]);
            
            # element is an element in the stabiliser that is not in the
            # group generated by newSGS
            AddSet(partialSGS, element);
            
            # Check if new generator fixes all points
            
            numFixedBasePoints := 0;
            for point in partialBase do
                if action(element[1], point) = point then
                    numFixedBasePoints := numFixedBasePoints + 1;
                else
                    break;
                fi;
            od;
            
            if numFixedBasePoints = Length(partialBase) then
                # find point that is moved by the Schreier generator
                # add point and compute its Schreier tree
                
                MATRIXSS_DebugPrint(6, ["Base : ", partialBase]);
                point := MATRIXSS_NewBasePoint(partialBase, element,
                                 action, identity, field);
                Add(partialBase, point);
                
                MATRIXSS_DebugPrint(3, ["Adding new base point : ", point]);           
                MATRIXSS_DebugPrint(8, ["Base : ", partialBase]);
                
                # recompute Schreier trees when base has changed
                # do we really need to recompute all?
                schreierTrees := MATRIXSS_GetSchreierTrees(partialSGS, 
                                         partialBase, 
                                         field, action, identity);             
            fi;
        fi;
        
    until false;
    
end;

# An implementation of the Schreier-Sims algorithm, for matrix groups
MatrixSchreierSims := function(G)
    local generators, base, trees, points, element, gens, list;

    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    gens := GeneratorsOfGroup(G);
    points := FullRowSpace(FieldOfMatrixGroup(G), DimensionOfMatrixGroup(G));
    base := BasisVectors(CanonicalBasis(points));
    
    # we make generators a partial strong generating set which also contain
    # inverses of all elements
    generators := [];
    for element in gens do
        if not element = Identity(G) then
            # Save reference to generator and its inverse
            # Then inverses need not be calculated later
            list := [element, Inverse(element)];
            AddSet(generators, Immutable(list));
            AddSet(generators, Immutable(Reversed(list)));
        fi;
    od;
                
    # Compute Schreier trees
    trees := MATRIXSS_GetSchreierTrees(generators, base, points, 
                     MATRIXSS_MSSAction, Identity(G));
    
    MATRIXSS_DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    MATRIXSS_DebugPrint(3, ["Partial base : ", base]);
    MATRIXSS_DebugPrint(3, ["Partial sgs : ", generators]);
    
    return MATRIXSS_SchreierSims(base, generators, trees, MATRIXSS_MSSAction, 
                   Identity(G), points);
end;

InstallGlobalFunction(MatrixGroupOrder, function(G)
    local ret, schreierTrees, tree, order;
    
    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    SetAssertionLevel(0);
    SetInfoLevel(MatrixSchreierSimsInfo, MATRIXSS_DEBUGLEVEL);
    
    MATRIXSS_DebugPrint(2, ["Entering MatrixGroupOrder"]);
    MATRIXSS_DebugPrint(2, ["Input group : ", G]);
    
    ret := MatrixSchreierSims1(G);
    
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

InstallGlobalFunction(MatrixSchreierSimsTest, function(maxDegree, maxFieldSize)
    local degree, power, size1, size2, primeNr, prime, 
          deg_time, prime_time, field_time, test_time;
    
    test_time := Runtime();
    degree := 2;
    while degree <= maxDegree do
        deg_time := Runtime();
        MATRIXSS_DebugPrint(1, ["Checking orders for degree : ", degree]);
        
        primeNr := 1;
        while primeNr <= 168 and Primes[primeNr] <= maxFieldSize do
            prime_time := Runtime();
            prime := Primes[primeNr];
            MATRIXSS_DebugPrint(1, ["\tChecking orders for prime : ", prime]);
            
            power := 1;
            while Primes[primeNr]^power <= maxFieldSize do
                field_time := Runtime();
                MATRIXSS_DebugPrint(1, ["\t\tChecking order for field size : ", 
                        prime^power]);
                
                size1 := Order(GL(degree, prime^power));
                size2 := MatrixGroupOrder(GL(degree, prime^power));
            
                if not size1 = size2 then
                    Print("Prime = ", prime, " power = ", power, "\n");
                    Print("Degree = ", degree, "\n");
                    Error("Order difference!");
                fi;
                
                power := power + 1;
                MATRIXSS_DebugPrint(1, ["\t\tTime for this field size : ", 
                        Runtime() - field_time]);
            od;
            
            primeNr := primeNr + 1;
            MATRIXSS_DebugPrint(1, ["\tTime for this prime : ", 
                    Runtime() - prime_time]);
        od;
        
        degree := degree + 1;
        MATRIXSS_DebugPrint(1, ["Time for this degree : ", Runtime() - deg_time]);
    od;
    
    Print("No order differences\n");
    Print("Total time for test : ", Runtime() - test_time, "\n");
    return true;
end);

#E