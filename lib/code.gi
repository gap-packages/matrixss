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
    return true;
end;

###############################################################################
##
#F MATRIXSS_PointAction(point, element)
##
## The action of a group element (a matrix) on a point (a row vector).
## The action is from the right 
## \beginitems
## `point' & The point (row vector) to act on.
##
## `element' & The group element (matrix) that acts.
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
## `point' & The point to act on. Must be a *normed* row vector.
##
## `element' & The group element (matrix) that acts.
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
    MATRIXSS_DebugPrint(9, ["Lookup ", point]);
    MATRIXSS_DebugPrint(9, ["Schreier tree ", schreierTree]);
    if KnowsDictionary(schreierTree, point) then
        return true;
    else
        return false;
    fi;    
end;

###############################################################################
##
#F MATRIXSS_RandomOrbitPoint(schreierTree)
##
## Returns a random point in the orbit given by `schreierTree'.
##
###############################################################################
MATRIXSS_RandomOrbitPoint := function(schreierTree)
    return RandomHashKey(schreierTree);
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
    local tree, edge;
    
    # Create Schreier vector
    tree := CallFuncList(NewDictionary, dictinfo);
    
    # Make the root point to itself 
    edge := [identity, identity];
    AddDictionary(tree, root, 
            Immutable(rec(Edge  := edge,
                          Depth := 0)));
    
    MATRIXSS_DebugPrint(9, ["Tree: ", tree]);
    return rec(Tree := tree, Labels := [], Height := 0);
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
  function(schreierTree, point, action, identity)
    local element, edge;
    
    if ValueOption("SimpleSchreierTree") = fail then
        # the group element and its inverse
        element := [identity, identity];
        
        repeat
            edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
            Assert(1, not IsBool(edge), "Point not in orbit!\n");
            
            MATRIXSS_DebugPrint(9, ["Found edge : ", edge]);
            if edge.Depth = 0 then
                return element;
            fi;
            
            point := action(point, edge.Edge[2]);
            element[1] := edge.Edge[1] * element[1];
            element[2] := element[2] * edge.Edge[2];
        until false;
    else
        # In this case the tree has height 1, so we are done with one
        # single lookup
        
        edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
        
        Assert(1, not IsBool(edge), "Point not in orbit!\n");
        return edge.Edge;
    fi;
end;

###############################################################################
##
#F MATRIXSS_RandomCosetRepresentative(schreierTree, action, identity)
##
## Return a random coset representative from the transversal defined by 
## `schreierTree'.
##
###############################################################################
MATRIXSS_RandomCosetRepresentative :=   
  function(schreierTree, action, identity)
    return MATRIXSS_OrbitElement(schreierTree, 
                   MATRIXSS_RandomOrbitPoint(schreierTree),
                   action, identity);
end;

###############################################################################
##
#F MATRIXSS_MonotoneTree(root, elements, action, identity, dictinfo)
##
## Create a monotone Schreier tree with given root and edge labels.
## \beginitems
## `root' & The root point of the tree.
##
## `elements' & The elements that, together with its inverses, will form the
##              edge labels of the tree.
##
## `action' & The action of the group on the point set.
##
## `identity' & the group identity (the identity matrix)
##
## `dictinfo' & The dictionary info for the tree, used to create hash function.
##
## \enditems
##
###############################################################################
MATRIXSS_MonotoneTree := function(root, elements, action, identity, dictinfo)
    local tree, orbit, level, newPoint, newPoints, gens, element,
          point, depth;
    
    # Add inverses to edge labels
    gens := [];
    for element in Reversed(elements) do
        Add(gens, Reversed(element));
    od;
    Append(gens, elements);
    MakeImmutable(gens);
    depth := 0;
    
    MATRIXSS_DebugPrint(4, ["Building monotone tree using gens : ", gens]);
    tree := MATRIXSS_CreateInitialSchreierTree(root, dictinfo, identity);
    orbit := MATRIXSS_GetOrbit(tree.Tree);
    for element in gens do
        newPoints := [];
        depth := depth + 1;
        
        for point in orbit do
            newPoint := action(point, element[1]);
            
            if not MATRIXSS_IsPointInOrbit(tree.Tree, newPoint) then
                AddDictionary(tree.Tree, newPoint, rec(
                        Edge  := Immutable(ShallowCopy(element)),
                        Depth := depth));
                Add(newPoints, newPoint);
            fi;
        od;
        Append(orbit, newPoints);
    od;         
    
    MATRIXSS_DebugPrint(4, ["Monotone tree built"]);
    return [tree.Tree, depth];
end;

###############################################################################
##
#F MATRIXSS_CreateShallowSchreierTree(orbitTree, root, generators, labels, action, identity, hash)
##
## Create a shallow Schreier tree, ie with at most logarithmic height.
## \beginitems
## `orbitTree' & Given tree representing the same orbit as the shallow Schreier
##               to be computed.
##
## `root' & The root point of the tree.
##
## `generators' & From this set will any new edge labels be taken.
##
## `labels' & The elements that, together with its inverses, will form the
##              edge labels of the tree.
##
## `action' & The action of the group on the point set.
##
## `identity' & the group identity (the identity matrix)
##
## `hash' & The dictionary info for the tree, used to create hash function.
##
## \enditems
##
###############################################################################
MATRIXSS_CreateShallowSchreierTree := 
  function(orbitTree, root, generators, labels, action, identity, hash)
    local tree, element, point, orbit, newPoint, generator, height, ret;
    
    tree := MATRIXSS_CreateInitialSchreierTree(root, hash, identity).Tree;
    height := 0;
    orbit := MATRIXSS_GetOrbit(tree);
    
    repeat
        element := fail;
        MATRIXSS_DebugPrint(3, ["Search for element to include"]);
        for point in orbit do
            for generator in generators do
                newPoint := action(point, generator[1]);
                if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                    element := MATRIXSS_OrbitElement(tree, point, 
                                       action, identity);
                    element[1] := element[1] * generator[1];
                    element[2] := generator[2] * element[2];
                    break;
                fi;
            od;
            if element <> fail then
                MATRIXSS_DebugPrint(3, ["Element found : ", element]);
                break;
            fi;
        od;
        if element <> fail then
            AddSet(labels, Immutable(element));
            ret := MATRIXSS_MonotoneTree(root, labels, action,
                           identity, hash);
            tree := ret[1];
            height := ret[2];
            orbit := MATRIXSS_GetOrbit(tree);
        fi;
    until Size(tree) = Size(orbitTree) or element = fail;
    MATRIXSS_DebugPrint(4, ["Size1 : ", Size(tree), " Size2 : ", 
            Size(orbitTree), " element : ", element]);
    Assert(1, Size(tree) = Size(orbitTree));
    return rec(Tree := tree, Labels := labels, Height := height);
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
  function(oldTree, root, generators, oldGenerators, action, identity, 
          dictinfo)
    local tree, point, generator, newPoint, newPoints, orbit, list, element,
          edge, height, labels;
            
    list := MATRIXSS_CopySchreierTree(oldTree.Tree, dictinfo);
    tree := list[1];
    orbit := ShallowCopy(list[2]);
    labels := ShallowCopy(oldTree.Labels);
    
    MATRIXSS_DebugPrint(4, ["Old orbit: ", orbit]);
    MATRIXSS_DebugPrint(4, ["Old gens: ", oldGenerators]);
    MATRIXSS_DebugPrint(4, ["Gens    : ", generators]);
    
    if ValueOption("SimpleSchreierTree") = fail then
        height := 0;
        #labels := oldTree.Labels;
        repeat
            newPoints := [];
            for point in orbit do
                edge := MATRIXSS_GetSchreierTreeEdge(tree, point);
                for generator in generators do
                    
                    # Add edges for all new points and new generators
                    if not MATRIXSS_IsPointInOrbit(oldTree.Tree, point) or
                       not generator in oldGenerators then
                        newPoint := action(point, generator[1]);
                        
                        if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                            AddDictionary(tree, newPoint, 
                                    rec(Edge := Immutable(ShallowCopy(generator)),
                                                Depth := edge.Depth + 1));
                            height := Maximum([height, edge.Depth + 1]);
                            Add(newPoints, newPoint);
                            AddSet(labels, generator);
                        fi;
                    fi;
                od;
            od;
            orbit := newPoints;
        until IsEmpty(orbit);
        
        if ValueOption("ShallowSchreierTree") <> fail then 
            return MATRIXSS_CreateShallowSchreierTree(tree, root, 
                           generators, oldTree.Labels, action, identity, 
                           dictinfo);
        fi;
    else
        repeat
            newPoints := [];
            for point in orbit do
                edge := MATRIXSS_GetSchreierTreeEdge(tree, point);
                for generator in generators do
                    
                    # Add edges for all new points and new generators
                    if not MATRIXSS_IsPointInOrbit(oldTree.Tree, point) or
                       not generator in oldGenerators then
                        newPoint := action(point, generator[1]);
                        
                        # Make Schreier tree have height 1
                        if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                            element := ShallowCopy(edge);
                            element[1] := element[1] * generator[1];
                            element[2] := generator[2] * element[2];
                            AddDictionary(tree, newPoint, 
                                    Immutable(rec(Edge := Immutable(element),
                                                  Depth := 1)));
                            Add(newPoints, newPoint);
                            AddSet(labels, generator);
                        fi;
                    fi;
                od;
            od;
            orbit := newPoints;
        until IsEmpty(orbit);
        height := 1;
    fi;          
    
    return rec(Tree := tree, Labels := Immutable(labels), Height := height);
end;    


###############################################################################
##
#F MATRIXSS_ComputeSchreierTree(tree, generators, action, root, hash, identity)
##
## Fill a Schreier tree that contains only the root.
## \beginitems
## `tree' & The Schreier tree to fill.
##
## `generators' & The generators for the group that gives rise to the orbit
##                represented by the Schreier tree.
##
## `action' & The action of the group on the point set.
##
## `root' & The root point of the tree.
##
## `hash' & The dictionary info for the tree, used to create hash function.
##
## `identity' & the group identity (the identity matrix)
##
## \enditems
##
###############################################################################
MATRIXSS_ComputeSchreierTree := 
  function(tree, generators, action, root, hash, identity)
    local point, generator, newPoint, newPoints, orbit, element, 
          elements, level, newElements, orbitTree, depth, edge, height, ret,
          labels;
        
    orbit := MATRIXSS_GetOrbit(tree);
    depth := 0;
    Assert(1, Length(orbit) = 1, "Tree not empty!\n");
        
    if ValueOption("SimpleSchreierTree") <> fail then
        labels := [];
        repeat
            newPoints := [];
            for point in orbit do
                edge := MATRIXSS_GetSchreierTreeEdge(tree, point);
                for generator in generators do
                    
                    newPoint := action(point, generator[1]);
                    
                    # Make Schreier tree edge labels be coset representatives
                    # rather than generators
                    if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                        element := ShallowCopy(edge.Edge);
                        element[1] := element[1] * generator[1];
                        element[2] := generator[2] * element[2];
                        AddDictionary(tree, newPoint, 
                                Immutable(rec(Edge := 
                                        element,
                                              Depth := 1)));
                        AddSet(labels, element);
                        Add(newPoints, newPoint);
                    fi;
                od;
            od;
            orbit := newPoints;
        until IsEmpty(orbit);
        return rec(Tree := tree, Labels := Immutable(labels), Height := 1);
    else
        labels := ShallowCopy(generators);
        repeat
            newPoints := [];
            depth := depth + 1;
            for point in orbit do
                for generator in generators do
                    
                    newPoint := action(point, generator[1]);
                    
                    if not MATRIXSS_IsPointInOrbit(tree, newPoint) then
                        AddDictionary(tree, newPoint, 
                                rec(Edge := generator,
                                    Depth := depth));
                        Add(newPoints, newPoint);
                    fi;
                od;
            od;
            orbit := newPoints;
        until IsEmpty(orbit);
        height := depth;
        
        if ValueOption("ShallowSchreierTree") <> fail then 
            return MATRIXSS_CreateShallowSchreierTree(tree, root, generators, 
                           [], action, identity, hash);
        fi;
        
        return rec(Tree := tree, Labels := Immutable(labels), 
                   Height := height);
    fi;
end;    

###############################################################################
##
#F MATRIXSS_GetSchreierTree(oldTree, root, generators, oldGenerators, action, hash, identity)
##
## Returns a Schreier tree. This routine encapsulates the other Schreier tree
## functions.
## \beginitems
## `oldTree' & The Schreier tree to extend, in case there should be an 
##             extension.
##
## `root' & The root point of the tree.
##
## `generators' & The generators for the group that gives rise to the orbit
##                represented by the Schreier tree.
##
## `oldGenerators' & The current generators (edge-labels) of `oldTree'.
##
## `action' & The action of the group on the point set.
##
## `hash' & The dictionary info for the tree, used to create hash function.
##
## `identity' & the group identity (the identity matrix)
##
## \enditems
##
###############################################################################
MATRIXSS_GetSchreierTree := 
  function(oldTree, root, generators, oldGenerators, action, hash, identity)
    
    if ValueOption("ExtendSchreierTree") <> fail then
        return MATRIXSS_ExtendSchreierTree(oldTree, root, generators, 
                       oldGenerators, action, identity, hash);
    else
        oldTree := MATRIXSS_CreateInitialSchreierTree(root, hash, identity);
        return MATRIXSS_ComputeSchreierTree(oldTree.Tree, generators, action, 
                       root, hash, identity);
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
  function(schreierTree, point, action, identity, freeGroup,
          genMap)
    local element, edge, word;
    
    if ValueOption("SimpleSchreierTree") = fail then
        # the group element and its inverse
        element := [identity, identity];
        word := [Identity(freeGroup), Identity(freeGroup)];
        repeat
            edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
            Assert(1, not IsBool(edge), "Point not in orbit!\n");
            
            if edge.Depth = 0 then
                return [element, word];
            fi;
            
            point := action(point, edge.Edge[2]);
            element[1] := edge.Edge[1] * element[1];
            element[2] := element[2] * edge.Edge[2];
            
            MATRIXSS_DebugPrint(4, ["Looking up ", edge.Edge[1], " in ",
                    genMap.Generators]);
            word[1] := genMap.FreeGenerators[Position(genMap.Generators, 
                               edge.Edge[1])] * word[1];
            word[2] := word[2] * 
                       genMap.FreeGenerators[Position(genMap.Generators, 
                               edge.Edge[2])];
        until false;
    else
        # In this case the tree has height 1, so we are done with one
        # single lookup
        
        edge := MATRIXSS_GetSchreierTreeEdge(schreierTree, point);
        Assert(1, not IsBool(edge), "Point not in orbit!\n");
        
        return [edge.Edge, [genMap.FreeGenerators[Position(genMap.Generators, 
                       edge.Edge[1])],
                       genMap.FreeGenerators[Position(genMap.Generators, 
                               edge.Edge[2])]]];
    fi;
end;

###############################################################################
##
#F MATRIXSS_Membership(ssInfo, element, identity)
##
## Check if an element belongs to a group, using sifting
## \beginitems
## `ssInfo' & Main information structure about our stabiliser chain. The Schreier
##          trees is used during the sifting.
##
## `element' & the element to check for membership 
##
## `identity' & group identity
## \enditems
##
###############################################################################
MATRIXSS_Membership := 
  function(ssInfo, element, identity)
    local level, residue, word, point;
    
    residue := ShallowCopy(element);
    
    # Find an expression of element in terms of the generators of the
    # groups in our stabiliser chain, using the Schreier trees
    for level in [1 .. Length(ssInfo)] do
        MATRIXSS_DebugPrint(9, ["residue: ", residue[1], "\nbase: ", 
                ssInfo[level].partialBase, "\naction", 
                ssInfo[level].action]);
        point := ssInfo[level].action(ssInfo[level].partialBase, 
                         residue[1]);
        
        if not MATRIXSS_IsPointInOrbit(ssInfo[level].schreierTree.Tree, 
                   point) then
            return [Immutable(residue), level];
        fi;
        
        word := MATRIXSS_OrbitElement(ssInfo[level].schreierTree.Tree, point, 
                        ssInfo[level].action, identity);
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
## `ssInfo' & Main information structure about our stabiliser chain. 
##            The Schreier trees is used during the sifting.
##
## `element' & the element to check for membership 
##
## `identity' & group identity
##
## `freeGroup' & the free group in which the sifted element will be expressed
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
        if not MATRIXSS_IsPointInOrbit(ssInfo[level].schreierTree.Tree, 
                   point) then
            return [Immutable(residue), level];
        fi;
        
        representative := 
          MATRIXSS_OrbitElement_ToddCoxeter(ssInfo[level].schreierTree.Tree, 
                  point, ssInfo[level].action, identity, 
                  ssInfo[level].freeGroup, ssInfo[level].genMap);
        
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
  function(schreierTree, generator, point, action, identity)
    local element1, element2, edge, inv_edge;
    
    element1 := MATRIXSS_OrbitElement(schreierTree, point, action,
                        identity);
    element2 := MATRIXSS_OrbitElement(schreierTree, 
                        action(point, generator[1]), action, identity);
    
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
  function(schreierTree, generator, point, action, identity, 
          freeGroup, genMap)
    local element1, element2, edge, inv_edge;
    
    element1 := MATRIXSS_OrbitElement_ToddCoxeter(schreierTree, point, action,
                        identity, freeGroup, genMap);
    element2 := MATRIXSS_OrbitElement_ToddCoxeter(schreierTree, 
                        action(point, generator[1]), action, identity,
                        freeGroup, genMap);
    
    edge := [element1[1][1] * generator[1] * element2[1][2],
             element1[2][1] * genMap.FreeGenerators[Position(genMap.Generators,
                     generator[1])] *
             element2[2][2]];
    inv_edge := [element2[1][1] * generator[2] * element1[1][2],
                 element2[2][1] * 
                 genMap.FreeGenerators[Position(genMap.Generators, 
                         generator[2])] *
                 element1[2][2]];
    
    return [[edge[1], inv_edge[1]], [edge[2], inv_edge[2], freeGroup]];
end;

###############################################################################
##
#F MATRIXSS_RandomSubproduct(elements, identity)
##
## Return a random subproduct of `elements'.
##
###############################################################################
MATRIXSS_RandomSubproduct := function(elements, identity)
    local element, product, perm, list, subProdGroup;
    
    list := [];
    
    subProdGroup := LookupDictionary(MATRIXSS_SubProdGroups, 
                                     Length(elements));
    if IsBool(subProdGroup) then
        subProdGroup := SymmetricGroup(Length(elements));
        AddDictionary(MATRIXSS_SubProdGroups, Length(elements), subProdGroup);
    fi;
    
    # permute elements randomly
    perm := PseudoRandom(subProdGroup);
    for element in [1 .. Length(elements)] do
        list[element] := elements[element^perm];
    od;
    
    # create random subproduct of elements
    product := [identity, identity];
    for element in elements do
        if Random([true, false]) then
            product[1] := product[1] * element[1];
            product[2] := element[2] * product[2];
        fi;
    od;
    
    return product;
end;

###############################################################################
##
#F MATRIXSS_RandomSchreierGenerator(schreierTree, elements, action, identity)
##
## Return a random Schreier generator constructed from the points in 
## `schreierTree' and the generators in `elements'.
##
###############################################################################
MATRIXSS_RandomSchreierGenerator :=
  function(schreierTree, elements, action, identity)
    local subsetLen, randomElements, subProd1, subProd2,
          elementBase, word, point, element;
    
    point := MATRIXSS_RandomOrbitPoint(schreierTree);
    subsetLen := Int(Random([0 .. Length(elements) - 1]) / 2);
    subProd1 := MATRIXSS_RandomSubproduct(elements, identity);
    
    randomElements := [];
    elementBase := ShallowCopy(elements);
    while Length(randomElements) < subsetLen do
        element := Random(elementBase);
        AddSet(randomElements, element);
        RemoveSet(elementBase, element);
    od;
    subProd2 := MATRIXSS_RandomSubproduct(randomElements, identity);
    
    word := Random([subProd1, subProd2]);
    
    return MATRIXSS_GetSchreierGenerator(schreierTree, word, point, 
                   action, identity);
end;

MATRIXSS_AugmentBase := function(ssInfo, newPoint, action, hash, identity)
    local levelStruct;
    
    MATRIXSS_DebugPrint(3, ["Extending base"]);
    
    # Extend base
    if ValueOption("AlternatingActions") <> fail then
        levelStruct := 
          rec(
              partialSGS := [],
              partialBase := NormedRowVector(newPoint),
              action := MATRIXSS_ProjectiveAction,
              points := hash[3],
              hash := hash,
              schreierTree := MATRIXSS_CreateInitialSchreierTree(
                      NormedRowVector(newPoint), hash, identity),
              oldSGS := AsSet([]),
              relations := [],
              IsIdentity := MATRIXSS_ProjectiveIsIdentity);
        Add(ssInfo, levelStruct); 
    fi;
    
    levelStruct := 
      rec(
          partialSGS := [],
          partialBase := newPoint,
          action := action,
          points := hash[3],
          hash := hash,
          schreierTree := MATRIXSS_CreateInitialSchreierTree(newPoint, 
                  hash, identity),
          oldSGS := AsSet([]),
          relations := [],
          IsIdentity := MATRIXSS_IsIdentity);
    Add(ssInfo, levelStruct); 
        
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
    local newPoint, length;
    
    MATRIXSS_DebugPrint(3, ["Finding new base point"]);
    
    length := Length(ssInfo);
    
    # Find new base point
    newPoint := MATRIXSS_NewBasePoint(badElement[1], identity, 
                        ssInfo[length].points);
    
    MATRIXSS_AugmentBase(ssInfo, newPoint, MATRIXSS_PointAction, 
            ssInfo[length].hash, identity);
end;

###############################################################################
##
#F MATRIXSS_GetPartialBaseSGS(generators, identity, field)
##
## Constructs a partial base and a partial SGS given a set of generators
## for a group. Returns the partial SGS and the `ssInfo' structure, 
## see "ssInfo".
## \beginitems
## `generators' & given set of generators
##
## `identity'   & group identity element (the identity matrix)
##
## `field'      & the vector space on which the group acts
## \enditems
##
###############################################################################
MATRIXSS_GetPartialBaseSGS := 
  function(generators, identity, field)
    local newPoint, element, gen, invGen, newSGS, level, dictinfo, 
          levelStruct, point, ssInfo;
    
###############################################################################
##
#V ssInfo
##
## Main structure holding information for the algorithm. This is not a global
## variable, but the same structure is used in all the variants of the 
## algorithm, but all members are not necessarily used.
##
## The structure `ssInfo' is a list of records, with a record for each level 
## in the algorithm, ie one record for each base point. New base points
## may of course be added to the base during the execution of the algorithm,
## and then a new record is added to the end of the list.
##
## The members of the record are:
## \beginitems
## `partialSGS' & the elements in the current partial SGS that fixes all
##              points at lower levels, or the whole partial SGS for the
##              first level
##
## `partialBase' & the base point for this level
##
## `action' & the action (function) at this level
##
## `points' & the field where the base point `partialBase' comes from
##
## `hash' & the hash function for the Schreier tree at this level
##
## `schreierTree' & the Schreier tree for this level, representing the
##                  basic orbit at this level, ie the orbit of `partialBase'
##                  under the action of `partialSGS' at the previous (lower) 
##                  level. Thus, the root of the tree is `partialBase'.
##
## `oldSGS' & the whole partial SGS at the last call of SchreierSims at
##            this level
##
## IsIdentity & the function to check if a point is the identity at this
##                level
##    
###############################################################################
    ssInfo := [];
    
    newSGS := [];
    
    if ValueOption("CleverBasePoints") <> fail then
        # Get a list of possibly good base points for this group
        MATRIXSS_BasePointStore := 
          BasisVectorsForMatrixAction(Group(generators));
    fi;
    
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
        MATRIXSS_DebugPrint(3, ["Length of ssInfo : ", Length(ssInfo)]);
        for level in [1 .. Length(ssInfo)] do
            MATRIXSS_DebugPrint(9, ["ssInfo at level ", level, " is ", 
                    ssInfo[level]]);
            point := ssInfo[level].partialBase;
            if ssInfo[level].action(point, element) = point then
                AddSet(ssInfo[level].partialSGS, gen);
                AddSet(ssInfo[level].partialSGS, invGen);
            else
                level := level - 1;
                break;
            fi;
        od;
        
        if level >= Length(ssInfo) then
            MATRIXSS_DebugPrint(8, ["Matrix ", element, " fixes all points "]);
            
            if Length(ssInfo) > 0 then
                MATRIXSS_ExtendBase(ssInfo, gen, identity);
            else              
                # Get initial point
                newPoint := MATRIXSS_NewBasePoint(gen[1], identity, field);
                dictinfo := [newPoint, true, field];
                
                MATRIXSS_AugmentBase(ssInfo, newPoint, MATRIXSS_PointAction, 
                        dictinfo, identity);
            fi;
            
            MATRIXSS_DebugPrint(7, ["Adding ", gen, " to SGS"]);
        fi;
        
        # Save reference to generator and its inverse
        # Then inverses need not be calculated later
        AddSet(newSGS, gen);
        AddSet(newSGS, invGen);
    od;
    
    return [newSGS, ssInfo];
end;

InstallGlobalFunction(MatrixGroupOrderStabChain, function(ssInfo)
    local order, levelStruct;
    
    order := 1;
    for levelStruct in ssInfo do
        order := order * MATRIXSS_GetOrbitSize(levelStruct.schreierTree.Tree);
    od;
    
    return order;
end);

###############################################################################
##
#M Size(G)
##
## A method for Size for finite matrix groups, that uses the implementation of
## Schreier-Sims algorithm in this package, ie it uses the StabChainMatrixGroup
## attribute to compute the order of `G'.
##
## This method is only installed if `MATRIXSS_TEST' is not defined when the 
## package is loaded.
##
###############################################################################
if IsBoundGlobal("MATRIXSS_TEST") then
    InstallMethod(Size, "for finite matrix groups", 
            [IsMatrixGroup and IsFinite], function(G)
        local ret, orbit, order;
        
        # Compute SGS and base and orbits (ie Schreier trees)
        ret := StabChainMatrixGroup(G);
        
        # Compute order of group using computed orbit sizes
        return MatrixGroupOrderStabChain(ret.SchreierStructure);
    end);
fi;

###############################################################################
#E