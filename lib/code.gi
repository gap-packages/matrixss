###############################################################################
##
#W    code.gi     The Matrix Schreier Sims package - General functions
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

# Ugly hack to take advantage of the clever base point selection strategy
# by O'Brien and Murray. This global variable holds the list of future base
# points.
MATRIXSS_BasePointStore := [];

MATRIXSS_DebugPrint := function(level, message)
    #Info(MatrixSchreierSimsInfo, level, message);
    #Info(MatrixSchreierSimsInfo, level, message, "\n");
    if level <= MATRIXSS_DEBUGLEVEL then
        CallFuncList(Print, Concatenation(message, ["\n"]));
    fi;
end;

# The action of a group element (a matrix) on a point (a row vector)
# The action is from the right 
MATRIXSS_PointAction := OnRight;

# The projective action of a matrix on a row vector
# The one-dimensional subspace corresponding to the point is represented
# by the corresponding normed row vector
MATRIXSS_ProjectiveAction := OnLines;

# Identity check when using the PointAction
MATRIXSS_IsIdentity := function(element, identity)
    return element = identity;
end;

# Identity check when using projective action (all scalar matrices are
# considered equal to the identity)
MATRIXSS_ProjectiveIsIdentity := function(element, identity)
    return ForAll(identity, i -> i = OnLines(i, element));
end;

# return all points (as a list) in the orbit of the point 
# which is root of the schreier tree
# The list is not necessarily sorted, and it is mutable
MATRIXSS_GetOrbit := function(schreierTree)
    return HashKeyEnumerator(schreierTree);
end;

# Check if a given point is in the orbit defined by the given Schreier tree
MATRIXSS_IsPointInOrbit := function(schreierTree, point)
    MATRIXSS_DebugPrint(9, ["Lookup ", point, " in Schreier tree ",
            schreierTree]);
    if not IsBool(LookupDictionary(schreierTree, point)) then
        return true;
    else
        return false;
    fi;    
end;

# Get the label of the edge originating at the given point, and directed 
# towards the root of the given Schreier tree
MATRIXSS_GetSchreierTreeEdge := function(schreierTree, point)
    return LookupDictionary(schreierTree, point);
end;

# Get size of orbit defined by the given Schreier tree
MATRIXSS_GetOrbitSize := function(schreierTree)
    return Size(schreierTree);
end;

# Create a Schreier tree containing only the root
MATRIXSS_CreateInitialSchreierTree := function(root, dictinfo, identity)
    local tree;
    
    # Create Schreier vector
    tree := NewDictionary(dictinfo[1], dictinfo[2], dictinfo[3]);
    
    # Make the root point to itself 
    AddDictionary(tree, root, Immutable([identity, identity]));
    
    return tree;
end;

# Create a Schreier tree containing only the root
MATRIXSS_CreateInitialSchreierTree_NoInverse := function(root, dictinfo, identity)
    local tree;
    
    # Create Schreier vector
    tree := NewDictionary(dictinfo[1], dictinfo[2], dictinfo[3]);
    
    # Make the root point to itself 
    AddDictionary(tree, root, identity);
    
    return tree;
end;

# Creates a copy of a whole Schreier tree
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

# Extends an existing Schreier tree by a given set of generators
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

# Fill a Schreier tree that contains only the root
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

# Compute the group element that connects the root of the Schreier tree to
# a given point
# this function assumes that the point actually is in the orbit described by
# the given Schreier tree
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

# Compute the group element that connects the root of the Schreier tree to
# a given point
# this function assumes that the point actually is in the orbit described by
# the given Schreier tree
MATRIXSS_OrbitElement_ToddCoxeter := 
  function(schreierTree, point, action, identity, IsIdentity, freeGroupHomo)
    local element, edge, word;
    
    if ValueOption("SimpleSchreierTree") = fail then
        # the group element and its inverse
        element := [identity, identity];
        word := [Identity(PreImage(freeGroupHomo)),
                 Identity(PreImage(freeGroupHomo))];
        repeat
            edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
            
            Assert(1, not IsBool(edge), "Point not in orbit!\n");
            
            if IsIdentity(edge[1], identity) then
                return [element, word];
            fi;
            
            point := action(point, edge[2]);
            element[1] := edge[1] * element[1];
            element[2] := element[2] * edge[2];
            
            word[1] := PreImagesRepresentative(freeGroupHomo, 
                               edge[1]) * word[1];
            word[2] := word[2] * PreImagesRepresentative(freeGroupHomo, 
                               edge[2]);
        until false;
    else
        # In this case the tree has height 1, so we are done with one
        # single lookup
        
        edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
        
        Assert(1, not IsBool(edge), "Point not in orbit!\n");
        return [edge, [PreImagesRepresentative(freeGroupHomo, edge[1]),
                       PreImagesRepresentative(freeGroupHomo, edge[2])]];
    fi;
end;

# check if an element belongs to a group, using sifting
# ssInfo - main information structure about our base
# element - the element to check membership for
# identity - group identity
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

# check if an element belongs to a group, using sifting
# ssInfo - main information structure about our base
# element - the element to check membership for
# identity - group identity
MATRIXSS_Membership_ToddCoxeter := 
  function(ssInfo, element, identity)
    local level, residue, representative, point, mappedWord, gens1, gens2;
    
    residue := [element, [Identity(PreImage(ssInfo[1].freeGroupHomo)),
                       Identity(PreImage(ssInfo[1].freeGroupHomo)),
                       ssInfo[1].freeGroupHomo]];
    
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
        
        MATRIXSS_DebugPrint(9, ["homo : ", ssInfo[level].freeGroupHomo]);
        
        representative := 
          MATRIXSS_OrbitElement_ToddCoxeter(ssInfo[level].schreierTree, 
                  point, ssInfo[level].action, identity, 
                  ssInfo[level].IsIdentity, ssInfo[level].freeGroupHomo);
        
        MATRIXSS_DebugPrint(9, ["residue : ", residue]);
        MATRIXSS_DebugPrint(9, ["representative : ", representative]);
        
        residue[1][1] := residue[1][1] * representative[1][2];
        residue[1][2] := representative[1][1] * residue[1][2];
        
        gens1 := GeneratorsOfGroup(PreImages(ssInfo[level].freeGroupHomo));
        gens2 := GeneratorsOfGroup(PreImages(ssInfo[1].freeGroupHomo));
        Assert(1, Length(gens1) >= Length(gens2));
          
        mappedWord := 
          MappedWord(representative[2][2], gens1, gens2{[1 .. Length(gens1)]});
        residue[2][1] := residue[2][1] * mappedWord;
        
        mappedWord := 
          MappedWord(representative[2][1], gens1, gens2{[1 .. Length(gens1)]});
        residue[2][2] := mappedWord * residue[2][2];
    od;
    
    if residue[1][1] <> identity then
        level := Length(ssInfo) + 1;
        if ValueOption("AlternatingActions") <> fail then
            level := level + 1;
        fi;
    fi;
    
    return [Immutable(residue), level];
end; 

# Find a point not in base that is moved by element
# (element fixes the base)
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

# Create a Schreier generator for the stabiliser in the group which has 
# "generator" as one of its generators. The stabiliser fixes "point".
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

# Create a Schreier generator for the stabiliser in the group which has 
# "generator" as one of its generators. The stabiliser fixes "point".
MATRIXSS_GetSchreierGenerator_ToddCoxeter := 
  function(schreierTree, generator, point, action, identity, IsIdentity,
          freeGroupHomo)
    local element1, element2, edge, inv_edge;
    
    element1 := MATRIXSS_OrbitElement_ToddCoxeter(schreierTree, point, action,
                        identity, IsIdentity, freeGroupHomo);
    element2 := MATRIXSS_OrbitElement_ToddCoxeter(schreierTree, 
                        action(point, generator[1]), action, identity,
                        IsIdentity, freeGroupHomo);
    
    edge := [element1[1][1] * generator[1] * element2[1][2],
             element1[2][1] * 
             PreImagesRepresentative(freeGroupHomo, generator[1]) * 
             element2[2][2]];
    inv_edge := [element2[1][1] * generator[2] * element1[1][2],
                 element2[2][1] * 
                 PreImagesRepresentative(freeGroupHomo, generator[2]) * 
                 element1[2][2]];
    
    return [[edge[1], inv_edge[1]], [edge[2], inv_edge[2], freeGroupHomo]];
end;


# Add a new base point to the base, so that a given element is not in the
# stabiliser of the point
# ssInfo - main information structure for the current Schreier-Sims run
# badElement - the element that fixes all current base points
# identity - the group identity
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
          #genMap := NewDictionary(identity, true, Group(generators)),
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
              #genMap := NewDictionary(identity, true, Group(generators)),
              IsIdentity := MATRIXSS_ProjectiveIsIdentity);
        Add(ssInfo, levelStruct); 
    fi;
end;

# Construct a partial base and a partial SGS given a set of generators
# for a group.
# generators - given set of generators
# action - action to use when finding new base points
# IsIdentity - function to use when checking that a group element is the
# identity (under the given action)
# field - the vector space on which the group acts
# identity - the group identity
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
                      #genMap := NewDictionary(identity, true, Group(generators)),
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
                       #   genMap := NewDictionary(identity, true, Group(generators)),
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




#E