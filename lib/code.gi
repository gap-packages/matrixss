###############################################################################
#1
#W    code.gi     The Matrix Schreier Sims package - General functions
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
#H    Dev start : 2004-01-10 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
##
## These are the common functions used by all the variants of the Schreier-Sims
## algorithm implemented in the package.
##
###############################################################################

Revision.("matrixss/lib/code_gi") := 
  "@(#)$Id$";

###############################################################################
##
#F MATRIXSS_DebugPrint(level, message)
##
## Internal function for printing debug messages. Uses the internal variable
## `MATRIXSS_DEBUGLEVEL', see "MATRIXSS_DEBUGLEVEL", to determine if the 
## message should be printed.
##
###############################################################################
MATRIXSS_DebugPrint := function(level, message)
    #Info(MatrixSchreierSimsInfo, level, message);
    #Info(MatrixSchreierSimsInfo, level, message, "\n");
    if level <= MATRIXSS_DEBUGLEVEL then
        CallFuncList(Print, Concatenation(message, ["\n"]));
    fi;
end;

###############################################################################
##
#F MATRIXSS_PointAction(point, element)
##
## The action of a group element (a matrix) on a point (a row vector).
## The action is from the right 
## \beginitems
## point & The point (row vector) to act on.
##
## element & The group element (matrix) that acts.
## \enditems
##
###############################################################################
MATRIXSS_PointAction := OnRight;

###############################################################################
##
#F MATRIXSS_ProjectiveAction(point, element)
##
## The projective action of a matrix on a row vector.
## The one-dimensional subspace corresponding to the point is represented
## by the corresponding normed row vector
## \beginitems
## point & The point to act on. Must be a *normed* row vector.
##
## element & The group element (matrix) that acts.
## \enditems
##
###############################################################################
MATRIXSS_ProjectiveAction := OnLines;

###############################################################################
##
#F MATRIXSS_IsIdentity(element, identity)
##
## Identity check when using normal point action 
## \beginitems
## `element' & the group element to check if it is equal to identity
##
## `identity' & the group identity (the identity matrix)
## \enditems
##
###############################################################################
MATRIXSS_IsIdentity := function(element, identity)
    return element = identity;
end;

###############################################################################
##
#F MATRIXSS_ProjectiveIsIdentity(element, identity)
##
## Identity check when using projective action (all scalar matrices are
## considered equal to the identity)
## \beginitems
## `element' & the group element to check if it is equal to identity
##
## `identity' & the group identity (the identity matrix)
## \enditems
##
###############################################################################
MATRIXSS_ProjectiveIsIdentity := function(element, identity)
    return ForAll(identity, i -> i = OnLines(i, element));
end;

###############################################################################
##
#F MATRIXSS_GetOrbit(schreierTree)
##
## Return all points (as a list) in the orbit of the point which is root of 
## the Schreier tree, ie return all keys in the Dictionary.
## The list is not necessarily sorted, and it is mutable.
##
###############################################################################
MATRIXSS_GetOrbit := function(schreierTree)
    return HashKeyEnumerator(schreierTree);
end;

###############################################################################
##
#F MATRIXSS_IsPointInOrbit(schreierTree, point)
##
## Check if the given point is in the orbit defined by the given Schreier tree.
##
###############################################################################
MATRIXSS_IsPointInOrbit := function(schreierTree, point)
    MATRIXSS_DebugPrint(9, ["Lookup ", point, " in Schreier tree ",
            schreierTree]);
    if not IsBool(LookupDictionary(schreierTree, point)) then
        return true;
    else
        return false;
    fi;    
end;

###############################################################################
##
#F MATRIXSS_GetSchreierTreeEdge(schreierTree, point)
##
## Get the label of the edge originating at the given point, and directed 
## towards the root of the given Schreier tree.
##
###############################################################################
MATRIXSS_GetSchreierTreeEdge := function(schreierTree, point)
    return LookupDictionary(schreierTree, point);
end;

###############################################################################
##
#F MATRIXSS_GetOrbitSize(schreierTree)
##
## Get size of orbit defined by the given Schreier tree.
##
###############################################################################
MATRIXSS_GetOrbitSize := function(schreierTree)
    return Size(schreierTree);
end;

###############################################################################
##
#F MATRIXSS_CreateInitialSchreierTree(root, dictinfo, identity)
##
## Create a Schreier tree containing only the root.
## \beginitems
## `root' & The base point that is to be the root of the Schreier tree.
##
## `dictinfo' & Used when creating the Dictionary that is the Schreier tree
##
## `identity' & the group identity element
## \enditems
##
###############################################################################
MATRIXSS_CreateInitialSchreierTree := function(root, dictinfo, identity)
    local tree;
    
    # Create Schreier vector
    tree := NewDictionary(dictinfo[1], dictinfo[2], dictinfo[3]);
    
    # Make the root point to itself 
    AddDictionary(tree, root, Immutable([identity, identity]));
    
    return tree;
end;

###############################################################################
##
#F MATRIXSS_CopySchreierTree(tree, dictinfo)
##
## Creates a copy of a whole Schreier tree, ie of makes a copy of the 
## Dictionary.
## \beginitems
## `tree' & the Dictionary to copy
## 
## `dictinfo' & the dictinfo that was used when creating `tree'
## \enditems
##
###############################################################################
MATRIXSS_CopySchreierTree := function(tree, dictinfo)
    local copyTree, keys, value, key;
    
    # Make a copy of the hashtable
    # A simple ShallowCopy does not work for hash tables, so we must copy
    # all keys and values explicitly
    copyTree := NewDictionary(dictinfo[1], dictinfo[2], dictinfo[3]);
    keys := HashKeyEnumerator(tree);
    for key in keys do
        value := LookupDictionary(tree, key);
        AddDictionary(copyTree, key, value);
    od;
    
    return [copyTree, keys];
end;

###############################################################################
##
#F MATRIXSS_ExtendSchreierTree(oldTree, generators, oldGenerators, action, dictinfo)
##
# Extends an existing Schreier tree by a given set of generators
## \beginitems
## `oldTree' & The Schreier tree to extend, ie a Dictionary.
##
## `generators' & The generators for the group that gives rise to the orbit
##                represented by the Schreier tree.
##
## `oldGenerators' & The current generators (edge-labels) of `oldTree'.
##
## `action' & The action of the group on the point set.
##
## `dictinfo' & The Dictionary info used when `oldTree' was created.
## \enditems
##
###############################################################################
MATRIXSS_ExtendSchreierTree := 
  function(oldTree, generators, oldGenerators, action, dictinfo)
    local tree, point, generator, newPoint, newPoints, orbit, list, element;
    
    list := MATRIXSS_CopySchreierTree(oldTree, dictinfo);
    tree := list[1];
    orbit := ShallowCopy(list[2]);
    
    MATRIXSS_DebugPrint(4, ["Old orbit: ", orbit]);
    MATRIXSS_DebugPrint(4, ["Old gens: ", oldGenerators]);
    MATRIXSS_DebugPrint(4, ["Gens    : ", generators]);
    
    if ValueOption("SimpleSchreierTree") = fail then
        repeat
            newPoints := [];
            for point in orbit do
                for generator in generators do
                    
                    # Add edges for all new points and new generators
                    if not MATRIXSS_IsPointInOrbit(oldTree, point) or
                       not generator in oldGenerators then
                        newPoint := action(point, generator[1]);
                        
                        if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                            AddDictionary(tree, newPoint, generator);
                            Add(newPoints, newPoint);
                        fi;
                    fi;
                od;
            od;
            orbit := newPoints;
        until IsEmpty(orbit);
    else
        repeat
            newPoints := [];
            for point in orbit do
                for generator in generators do
                    
                    # Add edges for all new points and new generators
                    if not MATRIXSS_IsPointInOrbit(oldTree, point) or
                       not generator in oldGenerators then
                        newPoint := action(point, generator[1]);
                        
                        # Make Schreier tree have height 1
                        if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                            element := 
                              ShallowCopy(MATRIXSS_GetSchreierTreeEdge(tree, 
                                      point));
                            element[1] := element[1] * generator[1];
                            element[2] := generator[2] * element[2];
                            AddDictionary(tree, newPoint, 
                                    Immutable(element));
                            Add(newPoints, newPoint);
                        fi;
                    fi;
                od;
            od;
            orbit := newPoints;
        until IsEmpty(orbit);
    fi;          
    
    return tree;
end;    

###############################################################################
##
#F MATRIXSS_ComputeSchreierTree(tree, generators, action)
##
## Fill a Schreier tree that contains only the root.
## \beginitems
## `tree' & The Schreier tree to fill, ie a Dictionary.
##
## `generators' & The generators for the group that gives rise to the orbit
##                represented by the Schreier tree.
##
## `action' & The action of the group on the point set.
## \enditems
##
###############################################################################
MATRIXSS_ComputeSchreierTree := 
  function(tree, generators, action)
    local point, generator, newPoint, newPoints, orbit, element;
    
    orbit := MATRIXSS_GetOrbit(tree);
    
    if ValueOption("SimpleSchreierTree") = fail then
        repeat
            newPoints := [];
            for point in orbit do
                for generator in generators do
                    
                    newPoint := action(point, generator[1]);
                    
                    if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                        AddDictionary(tree, newPoint, generator);
                        Add(newPoints, newPoint);
                    fi;
                od;
            od;
            orbit := newPoints;
        until IsEmpty(orbit);
    else
        repeat
            newPoints := [];
            for point in orbit do
                for generator in generators do
                    
                    newPoint := action(point, generator[1]);
                    
                    # Make Schreier tree have height 1
                    if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                        element := 
                          ShallowCopy(MATRIXSS_GetSchreierTreeEdge(tree, 
                                  point));
                        element[1] := element[1] * generator[1];
                        element[2] := generator[2] * element[2];
                        AddDictionary(tree, newPoint, Immutable(element));
                        Add(newPoints, newPoint);
                    fi;
                od;
            od;
            orbit := newPoints;
        until IsEmpty(orbit);
    fi;          
    
    return tree;
end;    

###############################################################################
##
#F MATRIXSS_OrbitElement(schreierTree, point, action, identity, IsIdentity)
##
## Compute the group element that connects the root of the Schreier tree to
## a given point. This function assumes that the point actually is in the 
## orbit described by the given Schreier tree.
## \beginitems
## `schreierTree' & Schreier tree for the orbit to use
## 
## `point' & the point to check if it is in the orbit
##
## `action' & the action that was used to create the Schreier tree
##
## `identity' & the group identity (the identity matrix)
##
## `IsIdentity' & function to use when checking if a group element is equal
##                to the identity
## \enditems
##
###############################################################################
MATRIXSS_OrbitElement := 
  function(schreierTree, point, action, identity, IsIdentity)
    local element, edge;
    
    if ValueOption("SimpleSchreierTree") = fail then
        # the group element and its inverse
        element := [identity, identity];
        
        repeat
            edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
            
            Assert(1, not IsBool(edge), "Point not in orbit!\n");
            
            if IsIdentity(edge[1], identity) then
                return element;
            fi;
            
            point := action(point, edge[2]);
            element[1] := edge[1] * element[1];
            element[2] := element[2] * edge[2];
        until false;
    else
        # In this case the tree has height 1, so we are done with one
        # single lookup
        
        edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
        
        Assert(1, not IsBool(edge), "Point not in orbit!\n");
        return edge;
    fi;
end;

###############################################################################
##
#F MATRIXSS_OrbitElement_ToddCoxeter(schreierTree, point, action, identity, IsIdentity, freeGroup, genMap)
##
## Special version of `MATRIXSS_OrbitElement', see "MATRIXSS_OrbitElement", 
## that also calculates the Word in the generators of the group element it 
## returns.
## 
## More specifically, it computes the Word of the generators of the 
## corresponding free group.
## \beginitems
## `schreierTree' & Schreier tree for the orbit to use
## 
## `point' & the point to check if it is in the orbit
##
## `action' & the action that was used to create the Schreier tree
##
## `identity' & the group identity (the identity matrix)
##
## `IsIdentity' & function to use when checking if a group element is equal
##                to the identity
## 
## `freeGroup' & corresponding free group to the group whose generators form
##               the set of edge labels of `schreierTree'
##
## `genMap' & list of 2 lists of the same length, the first being the edge 
##            labels of `schreierTree' (the generators of the corresponding 
##            group), and the second being the corresponding generators of
##            `freeGroup'
## \enditems
##
###############################################################################
MATRIXSS_OrbitElement_ToddCoxeter := 
  function(schreierTree, point, action, identity, IsIdentity, freeGroup,
          genMap)
    local element, edge, word;
    
    if ValueOption("SimpleSchreierTree") = fail then
        # the group element and its inverse
        element := [identity, identity];
        word := [Identity(freeGroup), Identity(freeGroup)];
        repeat
            edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
            
            Assert(1, not IsBool(edge), "Point not in orbit!\n");
            
            if IsIdentity(edge[1], identity) then
                return [element, word];
            fi;
            
            point := action(point, edge[2]);
            element[1] := edge[1] * element[1];
            element[2] := element[2] * edge[2];
            
            MATRIXSS_DebugPrint(8, ["Looking up ", edge[1], " in ",
                    genMap[1]]);
            word[1] := genMap[2][Position(genMap[1], edge[1])] * word[1];
            word[2] := word[2] * genMap[2][Position(genMap[1], edge[2])];
        until false;
    else
        # In this case the tree has height 1, so we are done with one
        # single lookup
        
        edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
        
        Assert(1, not IsBool(edge), "Point not in orbit!\n");
        return [edge, [genMap[2][Position(genMap[1], edge[1])],
                       genMap[2][Position(genMap[1], edge[2])]]];
    fi;
end;

###############################################################################
##
#F MATRIXSS_Membership(ssInfo, element, identity)
##
## Check if an element belongs to a group, using sifting
## \beginitems
## ssInfo & Main information structure about our stabiliser chain. The Schreier
##          trees is used during the sifting.
##
## element & the element to check for membership 
##
## identity & group identity
## \enditems
##
###############################################################################
MATRIXSS_Membership := 
  function(ssInfo, element, identity)
    local level, residue, word, point;
    
    residue := element;
    
    # Find an expression of element in terms of the generators of the
    # groups in our stabiliser chain, using the Schreier trees
    for level in [1 .. Length(ssInfo)] do
        MATRIXSS_DebugPrint(9, ["residue: ", residue[1], "\nbase: ", 
                ssInfo[level].partialBase, "\naction", 
                ssInfo[level].action]);
        point := ssInfo[level].action(ssInfo[level].partialBase, 
                         residue[1]);
        
        if not MATRIXSS_IsPointInOrbit(ssInfo[level].schreierTree, 
                   point) then
            return [Immutable(residue), level];
        fi;
        
        word := MATRIXSS_OrbitElement(ssInfo[level].schreierTree, point, 
                        ssInfo[level].action, identity, 
                        ssInfo[level].IsIdentity);
        residue[1] := residue[1] * word[2];
        residue[2] := word[1] * residue[2];
    od;
    
    level := Length(ssInfo) + 1;
    if ValueOption("AlternatingActions") <> fail then
        level := level + 1;
    fi;

    return [Immutable(residue), level];
end; 

###############################################################################
##
#F MATRIXSS_Membership_ToddCoxeter(ssInfo, element, identity, freeGroup)
##
## Special version of `MATRIXSS_Membership', see "MATRIXSS_Membership", that 
## also expresses the sifted group element as a word in the generators of a 
## given free group.
## \beginitems
## ssInfo & Main information structure about our stabiliser chain. The Schreier
##          trees is used during the sifting.
##
## element & the element to check for membership 
##
## identity & group identity
##
## freeGroup & the free group in which the sifted element will be expressed
## \enditems
##
###############################################################################
MATRIXSS_Membership_ToddCoxeter := 
  function(ssInfo, element, identity, freeGroup)
    local level, residue, representative, point, word, gens1, gens2;
    
    # Apart from calculating the group element, calculate the word of the
    # element in the generators
    word := [Identity(freeGroup), Identity(freeGroup), freeGroup];
    residue := [element, word];
    
    # Find an expression of element in terms of the generators of the
    # groups in our stabiliser chain, using the Schreier trees
    for level in [1 .. Length(ssInfo)] do
        MATRIXSS_DebugPrint(9, ["residue: ", residue[1], "\nbase: ", 
                ssInfo[level].partialBase, "\naction", 
                ssInfo[level].action]);
        point := ssInfo[level].action(ssInfo[level].partialBase, 
                         residue[1][1]);
        
        MATRIXSS_DebugPrint(9, ["Check if point ", point, " is in orbit"]);
        MATRIXSS_DebugPrint(9, ["Orbit is ", ssInfo[level].schreierTree]);
        
        MATRIXSS_DebugPrint(9, ["residue : ", residue]);
        if not MATRIXSS_IsPointInOrbit(ssInfo[level].schreierTree, 
                   point) then
            return [Immutable(residue), level];
        fi;
        
        representative := 
          MATRIXSS_OrbitElement_ToddCoxeter(ssInfo[level].schreierTree, 
                  point, ssInfo[level].action, identity, 
                  ssInfo[level].IsIdentity, ssInfo[level].freeGroup,
                  ssInfo[level].genMap);
        
        MATRIXSS_DebugPrint(9, ["residue : ", residue]);
        MATRIXSS_DebugPrint(9, ["representative : ", representative]);
        
        residue[1][1] := residue[1][1] * representative[1][2];
        residue[1][2] := representative[1][1] * residue[1][2];
        
        gens1 := GeneratorsOfGroup(ssInfo[level].freeGroup);
        gens2 := GeneratorsOfGroup(freeGroup);
        
        # Map the words to the same free group
        if Length(gens2) >= Length(gens1) and Length(gens1) > 0 then
            word := MappedWord(representative[2][2], gens1, 
                            gens2{[1 .. Length(gens1)]});
            residue[2][1] := residue[2][1] * word;
            
            word := MappedWord(representative[2][1], gens1, 
                            gens2{[1 .. Length(gens1)]});
            residue[2][2] := word * residue[2][2];
        else
            MATRIXSS_DebugPrint(2, ["1 : Gens1 : ", gens1]);
            MATRIXSS_DebugPrint(2, ["1 : Gens2 : ", gens2]);
        fi;
    od;
    
    if residue[1][1] <> identity then
        level := Length(ssInfo) + 1;
        if ValueOption("AlternatingActions") <> fail then
            level := level + 1;
        fi;
    fi;
    
    return [Immutable(residue), level];
end; 

###############################################################################
##
#F MATRIXSS_NewBasePoint(element, identity, field)
##
## Find a point not in base that is moved by the given element 
## (which fixes the base)
## \beginitems
## `element' & the bad element that fixes the whole base
## 
## `identity' & the group identity (the identity matrix)
## 
## `field'    & the vector space on which the group acts
## \enditems
##
###############################################################################
MATRIXSS_NewBasePoint := function(element, identity, field)
    local basis, point, basePoint, i, j, length;
    
    if ValueOption("CleverBasePoints") <> fail then
        if not IsEmpty(MATRIXSS_BasePointStore) then
            point := MATRIXSS_BasePointStore[1];
            MATRIXSS_BasePointStore := 
              MATRIXSS_BasePointStore{[2 .. Length(MATRIXSS_BasePointStore)]};
            return point;
        fi;
    fi;
    
    MATRIXSS_DebugPrint(3, ["Matrix that fixes whole base: ", element]);
    MATRIXSS_DebugPrint(5, ["Point field: ", field]);
    if not IsList(field) then
        field := BasisVectors(CanonicalBasis((field)));
    fi;
    
    length := Length(element);
    for i in [1 .. length] do
        for j in [1 .. length] do
            MATRIXSS_DebugPrint(8, ["Checking matrix element: ", element[i][j]]);           
            # If the element is not a diagonal matrix
            if i <> j and not IsZero(element[i][j]) then
                basePoint := ZeroMutable(field[1]);
                MATRIXSS_DebugPrint(6, ["Basepoint: ", basePoint]);
                basePoint[i] := 
                  One(FieldOfMatrixGroup(Group(identity)));
                MATRIXSS_DebugPrint(6, ["Basepoint: ", basePoint]);
                return Immutable(basePoint);
            fi;
        od;
    od;
    
    for i in [1 .. length] do
        for j in [1 .. length] do
            
            # If the element is not a scalar matrix
            if i <> j and element[i][i] <> element[j][j] then
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
    if element[1][1] <> One(FieldOfMatrixGroup(Group(identity))) then
        basePoint := ZeroMutable(field[1]);
        basePoint[1] := 
          One(FieldOfMatrixGroup(Group(identity)));
        return Immutable(basePoint);
    fi;
    
    Error("No new base point found!\n");
    return fail;
end;

###############################################################################
##
#F MATRIXSS_GetSchreierGenerator(schreierTree, generator, point, action, identity, IsIdentity)
##
## Creates a Schreier generator for the stabiliser in the group which has 
## `generator' as one of its generators. The stabiliser fixes `point' under
## `action'.
##
###############################################################################
MATRIXSS_GetSchreierGenerator := 
  function(schreierTree, generator, point, action, identity, IsIdentity)
    local element1, element2, edge, inv_edge;
    
    element1 := MATRIXSS_OrbitElement(schreierTree, point, action,
                        identity, IsIdentity);
    element2 := MATRIXSS_OrbitElement(schreierTree, 
                        action(point, generator[1]), action, identity,
                        IsIdentity);
    
    edge := element1[1] * generator[1] * element2[2];
    inv_edge := element2[1] * generator[2] * element1[2];
    
    return [edge, inv_edge];
end;

###############################################################################
##
#F MATRIXSS_GetSchreierGenerator_ToddCoxeter(schreierTree, generator, point, action, identity, IsIdentity, freeGroup, genMap)
##
## Special version of `MATRIXSS_GetSchreierGenerator', see 
## "MATRIXSS_GetSchreierGenerator", that also returns the Schreier generator
## as a Word in the generators of `freeGroup', using `genMap' to map the
## generators to the free group. See "MATRIXSS_OrbitElement_ToddCoxeter".
##
###############################################################################
MATRIXSS_GetSchreierGenerator_ToddCoxeter := 
  function(schreierTree, generator, point, action, identity, IsIdentity,
          freeGroup, genMap)
    local element1, element2, edge, inv_edge;
    
    element1 := MATRIXSS_OrbitElement_ToddCoxeter(schreierTree, point, action,
                        identity, IsIdentity, freeGroup, genMap);
    element2 := MATRIXSS_OrbitElement_ToddCoxeter(schreierTree, 
                        action(point, generator[1]), action, identity,
                        IsIdentity, freeGroup, genMap);
    
    edge := [element1[1][1] * generator[1] * element2[1][2],
             element1[2][1] * genMap[2][Position(genMap[1], generator[1])] *
             element2[2][2]];
    inv_edge := [element2[1][1] * generator[2] * element1[1][2],
                 element2[2][1] * 
                 genMap[2][Position(genMap[1], generator[2])] *
                 element1[2][2]];
    
    return [[edge[1], inv_edge[1]], [edge[2], inv_edge[2], freeGroup]];
end;

###############################################################################
##
#F MATRIXSS_ExtendBase(ssInfo, badElement, identity)
##
## Add a new base point to the base, so that the given element is not in the
## stabiliser of the point
## \beginitems
## `ssInfo' & main information structure for the current Schreier-Sims run
##
## `badElement' & the element that fixes all current base points
##
## `identity' & the group identity
## \enditems
##
###############################################################################
MATRIXSS_ExtendBase := function(ssInfo, badElement, identity)
    local newPoint, length, levelStruct;
    
    MATRIXSS_DebugPrint(3, ["Finding new base point"]);
    
    length := Length(ssInfo);
    
    # Find new base point
    newPoint := MATRIXSS_NewBasePoint(badElement[1], identity, 
                        ssInfo[length].points);
    
    MATRIXSS_DebugPrint(3, ["Extending base"]);
    
    # Extend base
    levelStruct := 
      rec(
          partialSGS := [],
          partialBase := newPoint,
          action := MATRIXSS_PointAction,
          points := ssInfo[1].points,
          hash := ssInfo[length].hash,
          schreierTree := MATRIXSS_CreateInitialSchreierTree(newPoint, 
                  ssInfo[length].hash, identity),
          oldSGS := AsSSortedList([]),
          relations := [],
          IsIdentity := MATRIXSS_IsIdentity);
    Add(ssInfo, levelStruct); 

    if ValueOption("AlternatingActions") <> fail then
        levelStruct := 
          rec(
              partialSGS := [],
              partialBase := NormedRowVector(newPoint),
              action := MATRIXSS_ProjectiveAction,
              points := ssInfo[length].points,
              hash := ssInfo[length].hash,
              schreierTree := MATRIXSS_CreateInitialSchreierTree(
                      NormedRowVector(newPoint), ssInfo[length].hash, 
                      identity),
              oldSGS := AsSSortedList([]),
              relations := [],
              IsIdentity := MATRIXSS_ProjectiveIsIdentity);
        Add(ssInfo, levelStruct); 
    fi;
end;

###############################################################################
##
#F MATRIXSS_GetPartialBaseSGS(generators, ssInfo, identity, field)
##
## Constructs a partial base and a partial SGS given a set of generators
## for a group.
## \beginitems
## `generators' & given set of generators
##
## `ssInfo'     & main information structure for algorithm
##
## `identity'   & group identity element (the identity matrix)
##
## `field'      & the vector space on which the group acts
## \enditems
##
###############################################################################
MATRIXSS_GetPartialBaseSGS := 
  function(generators, ssInfo, identity, field)
    local newPoint, element, gen, invGen, newSGS, level, dictinfo, 
          levelStruct, point;
    
    newSGS := [];
    
    # we make a partial strong generating set which also contain
    # inverses of all elements
    for element in generators do
        if element = identity then
            continue;
        fi;
        
        MATRIXSS_DebugPrint(3, ["Considering generator ", element]);
        
        gen := Immutable([element, Inverse(element)]);
        invGen := Immutable(Reversed(gen));
        level := 1;
        while level <= Length(ssInfo) do
            MATRIXSS_DebugPrint(9, ["ssInfo at level ", level, " is ", 
                    ssInfo[level]]);
            point := ssInfo[level].partialBase;
            if ssInfo[level].action(point, element) = point then
                AddSet(ssInfo[level].partialSGS, gen);
                AddSet(ssInfo[level].partialSGS, invGen);
            else
                break;
            fi;
            level := level + 1;
        od;
        
        if level >= Length(ssInfo) then
            MATRIXSS_DebugPrint(8, ["Matrix ", element, " fixes all points "]);
            
            if Length(ssInfo) > 0 then
                MATRIXSS_ExtendBase(ssInfo, gen, identity);
            else              
                # Get initial point
                newPoint := MATRIXSS_NewBasePoint(gen[1], identity, field);
                
                dictinfo := [newPoint, true, field];
                
                levelStruct := 
                  rec(
                      partialSGS := [],
                      partialBase := newPoint,
                      action := MATRIXSS_PointAction,
                      points := field,
                      hash := dictinfo,
                      schreierTree := 
                      MATRIXSS_CreateInitialSchreierTree(newPoint, 
                              dictinfo, identity),
                      oldSGS := AsSSortedList([]),
                      relations := [],
                      IsIdentity := MATRIXSS_IsIdentity);
                Add(ssInfo, levelStruct); 

                if ValueOption("AlternatingActions") <> fail then
                    levelStruct := 
                      rec(
                          partialSGS := [],
                          partialBase := NormedRowVector(newPoint),
                          action := MATRIXSS_ProjectiveAction,
                          points := NormedRowVectors(field),
                          hash := dictinfo,
                          schreierTree := 
                          MATRIXSS_CreateInitialSchreierTree(
                                  NormedRowVector(newPoint), dictinfo, 
                                  identity),
                          oldSGS := AsSSortedList([]),
                          relations := [],
                          IsIdentity := MATRIXSS_ProjectiveIsIdentity);
                    
                    Add(ssInfo, levelStruct); 
                fi;
            fi;
            
            MATRIXSS_DebugPrint(7, ["Adding ", gen, " to SGS"]);
            
        fi;
        
        # Save reference to generator and its inverse
        # Then inverses need not be calculated later
        AddSet(newSGS, gen);
        AddSet(newSGS, invGen);
    od;
    
    return newSGS;
end;

###############################################################################
##
#F MATRIXSS_ComputeOrder(ssInfo)
##
## Computes the order of the group defined by the given Schreier trees, see
## "ssInfo".
##
###############################################################################
MATRIXSS_ComputeOrder := function(ssInfo)
    local order, levelStruct;
    
    order := 1;
    for levelStruct in ssInfo do
        order := order * MATRIXSS_GetOrbitSize(levelStruct.schreierTree);
    od;
    
    return order;
end;

###############################################################################
##
#M Size(G)
##
## A method for Size for finite matrix groups, that uses the implementation of
## Schreier-Sims algorithm in this package, ie it uses the StabChainMatrixGroup
## attribute to compute the order of `G'.
##
## This method is only installed if MATRIXSS_TEST is not defined when the 
## package is loaded.
##
###############################################################################
if MATRIXSS_TEST = fail then
    InstallMethod(Size, "for finite matrix groups", 
            [IsMatrixGroup and IsFinite], function(G)
        local ret, orbit, order;
        
        # Compute SGS and base and orbits (ie Schreier trees)
        ret := StabChainMatrixGroup(G);
        
        # Compute order of group using computed orbit sizes
        return MATRIXSS_ComputeOrder(ret);
    end);
fi;

###############################################################################
#E