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
#MakeReadOnlyGlobal("MATRIXSS_DEBUGLEVEL");

DebugPrint := function(level, message)
    if level <= MATRIXSS_DEBUGLEVEL then
        CallFuncList(Print, message);
        Print("\n");
    fi;
end;

# The action of a group element g (a matrix) on a point p (a row vector)
# The action is from the right
# Used in Matrix Schreier-Sims
MSSAction := function(element, point)
    
    if not IsMatrix(element) or not IsRowVector(point) then
        Error("<element> must be a matrix and <point> must be a row vector");
    fi;
    
    return point * element;
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
ComputeSchreierTree := function(generators, root, action, tree)
    local newElement, childElements, generator, child, elements, element;
    
    elements := [root];
    repeat
        childElements := [];
        for element in elements do
            for generator in generators do
                newElement := action(generator, element);
                
                # check that the element is not already in the Schreier tree
                # (do not create cycles)
                if IsBool(GetSchreierTreeEdge(tree, newElement)) then
                    AddHashEntry(tree, newElement, generator);
                    Add(childElements, newElement);
                fi;
            od;    
        od;
        elements := childElements;
    until IsEmpty(elements);
    
    return tree;
end;      
    
# Computes a Schreier tree 
# root - root of the Schreier tree
# generators - generators for group
# points - point set where root comes from
# action - the action used to create the tree
#
# This is just a computation of a spanning tree for a connected component
SchreierTree := function(generators, points, root, action, identity)
    local tree;
    
    # Create Schreier vector
    tree := SparseHashTable(SparseIntKey(points, root));
    
    # Make the root point to itself 
    AddHashEntry(tree, root, identity);
    
    # Fill Schreier vector
    return ComputeSchreierTree(generators, root, action, tree);
end;

# Compute the group element that connects the root of the Schreier tree to
# a given point
# this function assumes that the point actually is in the orbit described by
# the given Schreier tree
OrbitElement := function(schreierTree, point, action, identity)
    local element, edge;
    
    element := identity;
    
    repeat
        edge := GetSchreierTreeEdge(schreierTree, point);
        
        Assert(1, not IsBool(edge), "Point not in orbit!\n");
        
        if edge = identity then
            return element;
        fi;
        
        point := action(Inverse(edge), point);
        element := edge * element;
    until false;
end;

# return generating set of a point stabiliser
# generators - generating set for group to compute in
# root - point to compute stabiliser of
# schreierTree - should have root as root
# action - the action to use
# identity - the identity element in the group
Stabiliser := function(generators, root, schreierTree, action, identity)
    local stabiliser, orbit, element1, element2, generator, point, edge;
    
    stabiliser := [];
    orbit := GetOrbit(schreierTree);
    
    # Compute Schreier generators of the stabiliser
    for generator in generators do
        for point in orbit do
            if not GetSchreierTreeEdge(schreierTree, point) = generator then
                element1 := OrbitElement(schreierTree, point, 
                                    action, identity);
                element2 := OrbitElement(schreierTree, 
                                    action(generator, point), 
                                    action, identity);
                
                edge := element1 * generator * Inverse(element2);
                
                # Important to make sure no generator is the identity
                if not edge = identity then
                    AddSet(stabiliser, edge);
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
Membership2 := function(schreierTrees, base, sgs, element, action, identity)
    local level, residue, word, point;
    
    residue := element;
    for level in [1 .. Length(base)] do
        point := action(residue, base[level]);
        
        if not IsPointInOrbit(schreierTrees[level], point) then
            return false;
        fi;
        
        word := OrbitElement(schreierTrees[level], point, action, identity);
        residue := residue * Inverse(word);
    od;
    
    if residue = identity then
        return true;
    else
        return false;
    fi;
end; 
    
# check if an element belongs to a group, using sifting
# schreierTrees - trees with the base points as roots, and with generators from
# the stabiliser of the previos point
# base - base for the group under consideration
# sgs - strong generating set for the stabiliser chain corresponding to base
# element - the element to check membership for
# action - the action to use
Membership := function(schreierTrees, base, sgs, element, action, identity)
    local point, generator, numFixedBasePoints, word, childElement;
                     
    numFixedBasePoints := 0;
    for point in base do
        if action(element, point) = point then
            numFixedBasePoints := numFixedBasePoints + 1;
        else
            break;
        fi;
    od;
        
    # basis case
    if numFixedBasePoints = Length(base) then
        if element = identity then
            return identity;
        else
            return fail;
        fi;
    fi;
    
    point := action(element, base[numFixedBasePoints + 1]);
    if not IsPointInOrbit(schreierTrees[numFixedBasePoints + 1], 
               point) then
        return fail;
    else
        word := OrbitElement(schreierTrees[numFixedBasePoints + 1],
                        point, action, identity);
        
        childElement := Membership(schreierTrees, base, sgs,
                                element * Inverse(word), action, identity);
        
        if not IsBool(childElement) then
            return childElement * word;
        else
            return fail;
        fi;
    fi;
end;

# find a point not in base that is moved by element
# (element fixes the base)
NewBasePoint := function(base, element, action, identity, field)
    local basis, point, basePoint;
    
    basePoint := fail;
    basis := BasisVectors(CanonicalBasis(field));
    
    for point in basis do
        if not action(element, point) = point and not point in base then
            basePoint := point;
            break;
        fi;
    od;
    
    if IsBool(basePoint) then
        Error("No new base point found!\n");
    fi;
    
    return basePoint;
end;

# Compute Schreier trees for the stabiliser chain generated by sgs and
# associated to base.
# Tree j is created with base point j as root, using generators in sgs that fix
# base points 1 .. j-1 (or all generators, in case of j = 1)
GetSchreierTrees := function(sgs, base, field, action, identity)
    local generators, point, element, schreierTrees, gens;
    
    # Create mutable generator list
    generators := [];
    UniteSet(generators, sgs);
    
    schreierTrees := [];
    for point in base do
        Add(schreierTrees, 
            SchreierTree(generators, field, point, action, identity));
        
        gens := ShallowCopy(generators);
        for element in generators do
            if not action(element, point) = point then
                RemoveSet(gens, element);
            fi;
        od;
        generators := gens;
    od;
    
    # If sgs is a strong generating set, then no generator fixes all points
    Assert(3, IsEmpty(generators), "Not a real base!\n");
    
    return schreierTrees;
end;

# Recursive Schreier-Sims over matrix group generated by S
# partialBase - a partial base
# partialSGS - a partial strong generating set
SchreierSims := function(partialBase, partialSGS, schreierTrees, action,
                        identity, field)
                local newBase, newSGS, newSchreierTrees, ret, element, point, 
                      orbit, generator, 
                      numFixedBasePoints, stabiliser, generators;
    
    if IsEmpty(partialSGS) then
        DebugPrint(2, ["Basis case, no generators"]);
        return [partialBase, partialSGS, schreierTrees];
    fi;
    
    repeat        
        DebugPrint(3, ["At current level"]);
        DebugPrint(3, ["Partial base : ", partialBase]);
        DebugPrint(3, ["Partial sgs : ", partialSGS]);
        
        DebugPrint(2, ["Getting partial base and sgs for next lower level"]);
        newBase := partialBase{[2 .. Length(partialBase)]};
        newSchreierTrees := schreierTrees{[2 .. Length(schreierTrees)]};
        newSGS  := [];
        for element in partialSGS do
            if action(element, partialBase[1]) = partialBase[1] then
                AddSet(newSGS, element);
            fi;
        od;
        
        DebugPrint(3, ["Partial base : ", newBase]);
        DebugPrint(3, ["Partial sgs : ", newSGS]);
        
        DebugPrint(2, ["Updating to complete base and sgs"]);
        ret              := SchreierSims(newBase, newSGS, newSchreierTrees, 
                                    action, identity, field);
        newBase          := ret[1];
        newSGS           := ret[2];
        newSchreierTrees := ret[3];
        
        DebugPrint(3, ["Base : ", newBase]);
        DebugPrint(3, ["SGS : ", newSGS]);

        partialBase := partialBase{[1]};
        Append(partialBase, newBase);
        partialSGS := UnionSet(partialSGS, newSGS);
        
        # recompute Schreier trees, when base and SGS has changed
        # only recompute first tree? we only need first one in this recursion
        schreierTrees := GetSchreierTrees(partialSGS, partialBase, field, 
                                 action, identity);
        
        # get Schreier generators of stabiliser
        stabiliser := Stabiliser(partialSGS, partialBase[1],
                              schreierTrees[1], action, identity);

        # check each Schreier generators if they are in the group generated by
        # newSGS
        generators := [];
        for element in stabiliser do                
            if Membership2(newSchreierTrees, newBase, newSGS, 
                       element, action, identity) then
                AddSet(generators, element);
            else
                break;
            fi;
        od;
        
        if Length(generators) = Length(stabiliser) then
            DebugPrint(3, ["No extension needed"]);
            return [partialBase, partialSGS, schreierTrees];
        else
            DebugPrint(3, ["Extending SGS with : ", element]);
            
            # element is an element in the stabiliser that is not in the
            # group generated by newSGS
            AddSet(partialSGS, element);
            
            # Check if new generator fixes all points
            
            numFixedBasePoints := 0;
            for point in partialBase do
                if action(element, point) = point then
                    numFixedBasePoints := numFixedBasePoints + 1;
                else
                    break;
                fi;
            od;
            
            if numFixedBasePoints = Length(partialBase) then
                # find point that is moved by the Schreier generator
                # add point and compute its Schreier tree
                
                DebugPrint(6, ["Base : ", partialBase]);
                point := NewBasePoint(partialBase, element,
                                 action, identity, field);
                Add(partialBase, point);
                
                DebugPrint(3, ["Adding new base point : ", point]);           
                DebugPrint(8, ["Base : ", partialBase]);
                
                # recompute Schreier trees when base has changed
                # do we really need to recompute all?
                schreierTrees := GetSchreierTrees(partialSGS, partialBase, 
                                         field, action, identity);             
            fi;
        fi;
        
    until false;
    
end;


# An implementation of the Schreier-Sims algorithm, for matrix groups
MatrixSchreierSims := function(G)
    local generators, base, trees, points, element, gens;

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
            AddSet(generators, element);
            AddSet(generators, Inverse(element));
        fi;
    od;
            
    SetAssertionLevel(0);
    
    # Compute Schreier trees
    trees := GetSchreierTrees(generators, base, points, MSSAction, 
                     Identity(G));
    
    DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    DebugPrint(3, ["Partial base : ", base]);
    DebugPrint(3, ["Partial sgs : ", generators]);
    
    return SchreierSims(base, generators, trees, MSSAction, Identity(G), 
                   points);
end;

InstallGlobalFunction(MatrixGroupOrder, function(G)
    local ret, schreierTrees, tree, order;
    
    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    DebugPrint(2, ["Entering MatrixGroupOrder"]);
    DebugPrint(2, ["Input group : ", G]);
    
    ret := MatrixSchreierSims(G);
    
    DebugPrint(2, ["Base : ", ret[1]]);
    DebugPrint(2, ["SGS : ", ret[2]]);
    
    order := 1;
    for tree in ret[3] do
        DebugPrint(3, ["Orbit size : ", GetOrbitSize(tree)]);
        order := order * GetOrbitSize(tree);
    od;
    
    DebugPrint(2, ["Order : ", order]);

    return order;
end);

InstallGlobalFunction(MatrixSchreierSimsTest, function(maxDegree, maxFieldSize)
    local degree, power, size1, size2, primeNr, prime;
    
    degree := 2;
    while degree <= maxDegree do
        DebugPrint(1, ["Checking orders for degree : ", degree]);
        
        primeNr := 1;
        while primeNr <= 168 and Primes[primeNr] <= maxFieldSize do
            prime := Primes[primeNr];
            DebugPrint(1, ["\tChecking orders for prime : ", prime]);
            
            power := 1;
            while Primes[primeNr]^power <= maxFieldSize do
                DebugPrint(1, ["\t\tChecking order for field size : ", 
                        prime^power]);
                
                size1 := Order(GL(degree, prime^power));
                size2 := MatrixGroupOrder(GL(degree, prime^power));
            
                if not size1 = size2 then
                    Print("Prime = ", prime, " power = ", power, "\n");
                    Print("Degree = ", degree, "\n");
                    Error("Order difference!");
                fi;
                
                power := power + 1;
            od;
            
            primeNr := primeNr + 1;
        od;
        
        degree := degree + 1;
    od;
    
    Print("No order differences\n");
    return true;
end);

#E