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
          ComputeSchreierTree, ExtendSchreierTree, SchreierSims, 
          NewBasePoint, ProjectiveNewBasePoint, Membership, OrbitElement, 
          GetSchreierTreeEdge, GetOrbitSize, GetSchreierGenerator, 
          GetPartialBaseSGS, ExtendBase, IsIdentity, ProjectiveIsIdentity, 
          IsConstantList, CreateInitialSchreierTree, CopySchreierTree,
    # Local variables
          ssInfo, list, generators, level, points, element;

    DebugPrint := function(level, message)
        #Info(MatrixSchreierSimsInfo, level, message);
        #Info(MatrixSchreierSimsInfo, level, message, "\n");
        if level <= MATRIXSS_DEBUGLEVEL then
            CallFuncList(Print, Concatenation(message, ["\n"]));
        fi;
    end;

    # The action of a group element (a matrix) on a point (a row vector)
    # The action is from the right 
    PointAction := OnRight;
    
    # The projective action of a matrix on a row vector
    # The one-dimensional subspace corresponding to the point is represented
    # by the corresponding normed row vector
    ProjectiveAction := OnLines;
    
    # Identity check when using the PointAction
    IsIdentity := function(element, identity)
        return element = identity;
    end;
        
    # Identity check when using projective action (all scalar matrices are
    # considered equal to the identity)
    ProjectiveIsIdentity := function(element, identity)
      return ForAll(identity, i -> i = OnLines(i, element));
    end;
    
    # return all points (as a list) in the orbit of the point 
    # which is root of the schreier tree
    # The list is not necessarily sorted, and it is mutable
    GetOrbit := function(schreierTree)
        return HashKeyEnumerator(schreierTree);
    end;

    # Check if a given point is in the orbit defined by the given Schreier tree
    IsPointInOrbit := function(schreierTree, point)
        if not IsBool(LookupDictionary(schreierTree, point)) then
            return true;
        else
            return false;
        fi;    
    end;

    # Get the label of the edge originating at the given point, and directed 
    # towards the root of the given Schreier tree
    GetSchreierTreeEdge := function(schreierTree, point)
        return LookupDictionary(schreierTree, point);
    end;

    # Get size of orbit defined by the given Schreier tree
    GetOrbitSize := function(schreierTree)
        return Size(schreierTree);
    end;
    
    # Create a Schreier tree containing only the root
    CreateInitialSchreierTree := function(root, dictinfo, identity)
        local tree;
        
        # Create Schreier vector
	tree := NewDictionary(dictinfo[1], dictinfo[2], dictinfo[3]);
        
        # Make the root point to itself 
        AddDictionary(tree, root, Immutable([identity, identity]));
        
        return tree;
    end;
    
    # Creates a copy of a whole Schreier tree
    CopySchreierTree := function(tree, dictinfo)
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
    ExtendSchreierTree := 
      function(oldTree, generators, oldGenerators, action, dictinfo)
      local tree, point, generator, newPoint, newPoints, orbit, list, element;
     
      list := CopySchreierTree(oldTree, dictinfo);
      tree := list[1];
      orbit := ShallowCopy(list[2]);
      
      DebugPrint(4, ["Old orbit: ", orbit]);
      DebugPrint(4, ["Old gens: ", oldGenerators]);
      DebugPrint(4, ["Gens    : ", generators]);
      
      if ValueOption("SimpleSchreierTree") = fail then
          repeat
              newPoints := [];
              for point in orbit do
                  for generator in generators do
                      
                    # Add edges for all new points and new generators
                      if not IsPointInOrbit(oldTree, point) or
                         not generator in oldGenerators then
                          newPoint := action(point, generator[1]);
                          
                          if not IsPointInOrbit(tree, newPoint) then
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
                      if not IsPointInOrbit(oldTree, point) or
                         not generator in oldGenerators then
                          newPoint := action(point, generator[1]);
                          
                          # Make Schreier tree have height 1
                          if not IsPointInOrbit(tree, newPoint) then
                              element := ShallowCopy(GetSchreierTreeEdge(tree, 
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
  ComputeSchreierTree := 
    function(tree, generators, action)
      local point, generator, newPoint, newPoints, orbit, element;
      
      orbit := GetOrbit(tree);
      
      if ValueOption("SimpleSchreierTree") = fail then
          repeat
              newPoints := [];
              for point in orbit do
                  for generator in generators do
                      
                      newPoint := action(point, generator[1]);
                      
                      if not IsPointInOrbit(tree, newPoint) then
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
                      if not IsPointInOrbit(tree, newPoint) then
                          element := ShallowCopy(GetSchreierTreeEdge(tree, 
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
    OrbitElement := 
      function(schreierTree, point, action, identity, IsIdentity)
        local element, edge;
        
        if ValueOption("SimpleSchreierTree") = fail then
            # the group element and its inverse
            element := [identity, identity];
            
            repeat
                edge := GetSchreierTreeEdge(schreierTree, point);
                
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
            
            edge := GetSchreierTreeEdge(schreierTree, point);
                
            Assert(1, not IsBool(edge), "Point not in orbit!\n");
            return edge;
        fi;
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
            DebugPrint(9, ["residue: ", residue[1], "\nbase: ", 
                    ssInfo[level].partialBase, "\naction", 
                    ssInfo[level].action]);
            point := ssInfo[level].action(ssInfo[level].partialBase, 
                             residue[1]);
            
            if not IsPointInOrbit(ssInfo[level].schreierTree, point) then
                return [Immutable(residue), level];
            fi;
            
            word := OrbitElement(ssInfo[level].schreierTree, point, 
                            ssInfo[level].action, identity, 
                            ssInfo[level].IsIdentity);
            residue[1] := residue[1] * word[2];
            residue[2] := word[1] * residue[2];
        od;
        
        return [Immutable(residue), Length(ssInfo) + 1];
    end; 
        
    # Find a point not in base that is moved by element
    # (element fixes the base)
    NewBasePoint := function(element, identity, field)
        local basis, point, basePoint, i, j, length;
        
        DebugPrint(3, ["Matrix that fixes whole base: ", element]);
        DebugPrint(5, ["Point field: ", field]);
        field := AsList(field);
        
        length := Length(element);
        for i in [1 .. length] do
            for j in [1 .. length] do
                DebugPrint(8, ["Checking matrix element: ", element[i][j]]);           
                # If the element is not a diagonal matrix
                if i <> j and not IsZero(element[i][j]) then
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
    GetSchreierGenerator := 
      function(schreierTree, generator, point, action, identity, IsIdentity)
        local element1, element2, edge, inv_edge;
        
        element1 := OrbitElement(schreierTree, point, action,
                            identity, IsIdentity);
        element2 := OrbitElement(schreierTree, 
                            action(point, generator[1]), action, identity,
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
        newPoint := NewBasePoint(badElement[1], identity, 
                            ssInfo[length].points);
        
        DebugPrint(3, ["Extending base"]);
        
        # Extend base
        levelStruct := 
          rec(
              partialSGS := [],
              partialBase := newPoint,
              action := PointAction,
              points := ssInfo[1].points,
              hash := ssInfo[length].hash,
              schreierTree := CreateInitialSchreierTree(newPoint, 
                      ssInfo[length].hash, identity),
              oldSGS := AsSSortedList([]),
              IsIdentity := IsIdentity);
        Add(ssInfo, levelStruct); 

        if ValueOption("AlternatingActions") <> fail then
            levelStruct := 
              rec(
                  partialSGS := [],
                  partialBase := NormedRowVector(newPoint),
                  action := ProjectiveAction,
                  points := ssInfo[length].points,
                  hash := ssInfo[length].hash,
                  schreierTree := CreateInitialSchreierTree(
                          NormedRowVector(newPoint), ssInfo[length].hash, 
                          identity),
                  oldSGS := AsSSortedList([]),
                  IsIdentity := ProjectiveIsIdentity);
            Add(ssInfo, levelStruct); 
        fi;
    end;
    
    # The main Schreier-Sims function
    # ssInfo - main information structure for the current Schreier-Sims run
    # partialSGS - given partial strong generating set
    # level - the level of the call to Schreier-Sims
    # identity - the group identity
    SchreierSims := function(ssInfo, partialSGS, level, identity)
        local generator, point, orbit, strip, schreierGenerator, element, 
              action, recursiveLevel, schreierTree, SGS, oldSGS, points, 
              newPoint, oldSchreierTree, field, newBasePoint, oldOrbit,
              newInverseGenerator, newSGS, cosetTable, list, TC1, TC2;
        
        DebugPrint(2, ["Schreier-Sims at level ", level]);
        
        action := ssInfo[level].action;
        field  := ssInfo[level].points;
                
        # S^(i + 1) = ssInfo[i].partialSGS
        
        # Find S^(i)        
        # Find the generators from the partial SGS that fixes all points at
        # lower levels.
        if level > 1 then
            SGS := ShallowCopy(ssInfo[level - 1].partialSGS);
        else
            SGS := ShallowCopy(partialSGS);
        fi;
        MakeImmutable(SGS);
        
        DebugPrint(3, ["Saved SGS that fixes first ", level - 1, " points ",
                Length(SGS)]);
        
        DebugPrint(4, ["Base point : ", ssInfo[level].partialBase]);
        DebugPrint(9, ["Hash func : ", ssInfo[level].hash]);
        
        # Compute schreier tree for current level
        # Compute \Delta_i = \beta_i^(H^i) = \beta_i^(<S_i>)
        oldSchreierTree := ssInfo[level].schreierTree;
        
        if ValueOption("ExtendSchreierTree") <> fail then
            ssInfo[level].schreierTree := 
              ExtendSchreierTree(ssInfo[level].schreierTree, 
                      SGS, ssInfo[level].oldSGS, action, ssInfo[level].hash);
        else
            DebugPrint(7, ["Creating new empty Schreier Tree"]);
            ssInfo[level].schreierTree := 
              CreateInitialSchreierTree(ssInfo[level].partialBase,
                      ssInfo[level].hash, identity);
            DebugPrint(7, ["Filling new Schreier Tree"]);
            ssInfo[level].schreierTree :=
              ComputeSchreierTree(ssInfo[level].schreierTree, SGS, action);
        fi;
        
        orbit := Immutable(GetOrbit(ssInfo[level].schreierTree));
        
        DebugPrint(6, ["New Schreier Tree : ", ssInfo[level].schreierTree]);
        DebugPrint(6, ["Old Schreier Tree : ", oldSchreierTree]);
        Assert(1, not IsIdenticalObj(ssInfo[level].schreierTree, 
                oldSchreierTree));
        DebugPrint(4, ["Orbit size for level ", level, " is ", 
                Length(orbit)]); 
        
        # We now want to make sure that SGS also fixes the current level
        
        for point in orbit do
            for generator in SGS do
                
                # Avoid rechecking Schreier generators
                if not IsPointInOrbit(oldSchreierTree, point) or 
                   not generator in ssInfo[level].oldSGS then
                                        
                    # Compute Schreier generator g for current level
                    schreierGenerator := 
                      GetSchreierGenerator(ssInfo[level].schreierTree,
                              generator, point, action, identity,
                              ssInfo[level].IsIdentity);
                                        
                    DebugPrint(6, ["Schreier Generator : ", 
                            schreierGenerator]);
                                        
                    if ssInfo[level].IsIdentity(schreierGenerator[1], 
                               identity) then
                        continue;
                    fi;
                    
                    # Check if Schreier generator is in stabiliser at 
                    # the current level
                    # Check if g \in H^(i + 1) = <S^(i + 1)>
                    points := [level + 1 .. Length(ssInfo)];
                    strip := Membership(ssInfo{points},
                                     schreierGenerator, 
                                     identity);
                                        
                    # The drop-out level is in range 
                    # [1 .. Length(ssInfo) + 1 - level] 
                    # but we want the range given by points
                    strip[2] := strip[2] + level;
                    
                    DebugPrint(5, ["Dropout level : ", strip[2]]);
                    
                    if strip[1][1] <> identity then
                        DebugPrint(3, ["Residue found"]);
                        
                        # We have found a Schreier generator which is not in
                        # the stabiliser of the current level, and so we must
                        # add the residue of this generator to our partial SGS
                        # in order to make it into a real SGS
                        
                        newInverseGenerator := Immutable(Reversed(strip[1]));
                        
                        # Add residue to partial SGS
                        # This makes some levels incomplete and so we must
                        # recompute them recursively
                        AddSet(partialSGS, strip[1]);
                        AddSet(partialSGS, newInverseGenerator);
                                                
                        # Possibly extend the base if the Schreier generator
                        # fixes all points in our base
                        if strip[2] = Length(ssInfo) + 1 then
                            ExtendBase(ssInfo, strip[1], identity);
                        fi;
                        
                        # Update partial SGS at each level
                        for recursiveLevel in [level .. strip[2] - 1] do
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   strip[1]);
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   newInverseGenerator);
                        od;
                        
                        if ValueOption("AlternatingActions") <> fail then
                            strip[2] := strip[2] + 1;
                        fi;
                        
                        # We must not recompute all levels downward from the
                        # dropout level
                        for recursiveLevel in Reversed([level + 1 .. 
                                strip[2]]) do
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
            
            DebugPrint(3, ["Considering generator ", element]);
            
            gen := Immutable([element, Inverse(element)]);
            invGen := Immutable(Reversed(gen));
            level := 1;
            while level <= Length(ssInfo) do
                DebugPrint(9, ["ssInfo at level ", level, " is ", 
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
                DebugPrint(8, ["Matrix ", element, " fixes all points "]);
                
                if Length(ssInfo) > 0 then
                    ExtendBase(ssInfo, gen, identity);
                else              
                    # Get initial point
                    newPoint := NewBasePoint(gen[1], identity, field);
                    
                    dictinfo := [newPoint, true, field];
                    
                    levelStruct := 
                      rec(
                          partialSGS := [],
                          partialBase := newPoint,
                          action := PointAction,
                          points := field,
                          hash := dictinfo,
                          schreierTree := 
                          CreateInitialSchreierTree(newPoint, dictinfo, 
                                  identity),
                          oldSGS := AsSSortedList([]),
                          IsIdentity := IsIdentity);
                    Add(ssInfo, levelStruct); 

                    if ValueOption("AlternatingActions") <> fail then
                        levelStruct := 
                          rec(
                              partialSGS := [],
                              partialBase := NormedRowVector(newPoint),
                              action := ProjectiveAction,
                              points := NormedRowVectors(field),
                              hash := dictinfo,
                              schreierTree := 
                              CreateInitialSchreierTree(
                                      NormedRowVector(newPoint), dictinfo, 
                                      identity),
                              oldSGS := AsSSortedList([]),
                              IsIdentity := ProjectiveIsIdentity);
                        
                        Add(ssInfo, levelStruct); 
                    fi;
                fi;
                
                DebugPrint(7, ["Adding ", gen, " to SGS"]);
                
            fi;
            
            # Save reference to generator and its inverse
            # Then inverses need not be calculated later
            AddSet(newSGS, gen);
            AddSet(newSGS, invGen);
        od;
        
        return newSGS;
    end;


    ### MAIN Schreier-Sims 

    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    # Get initial set of generators, to be extended to a partial SGS
    generators := GeneratorsOfGroup(G);
    
    # The vector space on which the group acts
    points := FullRowSpace(FieldOfMatrixGroup(G), DimensionOfMatrixGroup(G));
    
    # Main structure holding information needed by the algorithm
    ssInfo := [];
        
    # ssInfo has a record for each level in the algorithm, and there is one
    # level for each base point. The members of the record are:
    #   partialSGS - the elements in the current partial SGS that fixes all
    #                points at lower levels, or the whole partial SGS for the
    #                first level
    #   partialBase - the base point for this level
    #   action - the action (function) at this level
    #   points - the field where the base points come from
    #   hash - the hash function for the Schreier tree at this level
    #   schreierTree - the Schreier tree for this level, representing the
    #                  basic orbit at this level, ie the orbit of the member
    #                  "partialBase" at this level, under the action of
    #                  "partialSGS" at the previous (lower) level
    #                  Thus, the root of the tree is "partialBase".
    #   oldSGS - the whole partial SGS at the last call of SchreierSims at
    #            this level
    #   IsIdentity - the function to check if a point is the identity at this
    #                level
    
    
    DebugPrint(3, ["Group generators : ", generators]);
    
    # Compute initial partial SGS and base and fill ssInfo
    generators := GetPartialBaseSGS(generators, ssInfo, Identity(G), points);
    
    DebugPrint(3, ["Partial sgs : ", generators]);
    DebugPrint(3, ["Initial base length : ", Length(ssInfo)]);
    DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    
    # Call Schreier-Sims algorithm for each level (starting from top)
    for level in Reversed([1 .. Length(ssInfo)]) do
        SchreierSims(ssInfo, generators, level, Identity(G));
    od;
    
    DebugPrint(2, ["Matrix Schreier-Sims done"]);
    
    # Create output structure
    list := [[], generators, []];
    for level in [1 .. Length(ssInfo)] do
        Add(list[1], ssInfo[level].partialBase);
        Add(list[3], ssInfo[level].schreierTree);
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
              GeneralOrthogonalGroup, SpecialOrthogonalGroup]);
#              GeneralUnitaryGroup, SpecialUnitaryGroup]);
    
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
            
        if size1 <> size2 then
            Print("Correct order: ", size1, "\n");
            Print("Computed order: ", size2, "\n");
            Print("Order difference!\n");
        fi;
    od;
    
    Print("No order differences\n");
    return true;
end);

InstallGlobalFunction(MatrixSchreierSimsBenchmark, function(maxDegree, 
        maxFieldSize, maxReeSize, maxSuzukiSize)
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