###############################################################################
#1
#W    stcs.gi  The Matrix Schreier Sims package 
#W             Schreier-Todd-Coxeter-Sims implementation
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-07-05 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
##
## These are the Schreier-Todd-Coxeter-Sims routines, ie Schreier-Sims 
## algorithm with additional calls to Todd-Coxeter coset enumeration to
## possibly speed up the process. It is known to be fast when the input is
## already a base and SGS, and therefore it is good for verifying a proposed
## base and SGS, for example the output of a probabilistic algorithm.
##
###############################################################################

Revision.("matrixss/lib/stcs_gi") := 
  "@(#)$Id$";

###############################################################################
##
#F MATRIXSS_SchreierToddCoxeterSims(ssInfo, partialSGS, level, identity)
##
## It has the same interface as the deterministic algorithm, see 
## "SchreierSims".
##
###############################################################################
MATRIXSS_SchreierToddCoxeterSims := function(ssInfo, partialSGS, level, 
                                            identity)
        local generator, point, orbit, strip, schreierGenerator, element, 
              action, recursiveLevel, schreierTree, SGS, oldSGS, points, 
              newPoint, oldSchreierTree, newBasePoint, oldOrbit,
              newInverseGenerator, newSGS, freeGroup, cosetTable, word,
              subgroupGens, relation, gens1, gens2, gens3, relations, 
              levelGroup, freeGroupHomo, gen, relHomo;
        
        MATRIXSS_DebugPrint(2, ["Schreier-Sims at level ", level]);
        
        action := ssInfo[level].action;
                
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
        
        # Compute free group for current level
        ssInfo[level].freeGroup := FreeGroup(Length(SGS));
        MATRIXSS_DebugPrint(6, ["Generators : ", List(SGS, i -> i[1])]);
        levelGroup := GroupWithGenerators(List(SGS, i -> i[1]), identity);
        
        # Save mapping between generators of current level group and free group
        gens1 := ShallowCopy(GeneratorsOfGroup(levelGroup));
        gens2 := ShallowCopy(GeneratorsOfGroup(ssInfo[level].freeGroup));
        SortParallel(gens1, gens2);
        ssInfo[level].genMap := Immutable([gens1, gens2]);
        
        # Compute homomorphism for use in coset enumeration
        freeGroupHomo := 
          GroupHomomorphismByImagesNC(ssInfo[level].freeGroup, levelGroup,
                  GeneratorsOfGroup(ssInfo[level].freeGroup), 
                  GeneratorsOfGroup(levelGroup));
        
        MATRIXSS_DebugPrint(3, ["Saved SGS that fixes first ", level - 1, 
                " points ", Length(SGS)]);
        
        MATRIXSS_DebugPrint(4, ["Base point : ", ssInfo[level].partialBase]);
        MATRIXSS_DebugPrint(9, ["Hash func : ", ssInfo[level].hash]);
        
        # Compute schreier tree for current level
        # Compute \Delta_i = \beta_i^(H^i) = \beta_i^(<S_i>)
        oldSchreierTree := ssInfo[level].schreierTree;
        
        if ValueOption("ExtendSchreierTree") <> fail then
            ssInfo[level].schreierTree := 
              MATRIXSS_ExtendSchreierTree(ssInfo[level].schreierTree, 
                      SGS, ssInfo[level].oldSGS, action, ssInfo[level].hash);
        else
            MATRIXSS_DebugPrint(7, ["Creating new empty Schreier Tree"]);
            ssInfo[level].schreierTree := 
              MATRIXSS_CreateInitialSchreierTree(ssInfo[level].partialBase,
                      ssInfo[level].hash, identity);
            MATRIXSS_DebugPrint(7, ["Filling new Schreier Tree"]);
            ssInfo[level].schreierTree :=
              MATRIXSS_ComputeSchreierTree(ssInfo[level].schreierTree, SGS, 
                      action);
        fi;
        
        orbit := Immutable(MATRIXSS_GetOrbit(ssInfo[level].schreierTree));
        
        MATRIXSS_DebugPrint(6, ["New Schreier Tree : ", 
                ssInfo[level].schreierTree]);
        MATRIXSS_DebugPrint(6, ["Old Schreier Tree : ", oldSchreierTree]);
        Assert(1, not IsIdenticalObj(ssInfo[level].schreierTree, 
                oldSchreierTree));
        MATRIXSS_DebugPrint(4, ["Orbit size for level ", level, " is ", 
                Length(orbit)]); 
        
        # We now want to make sure that SGS also fixes the current level
                
        for point in orbit do
            for generator in SGS do
                
                # Avoid rechecking Schreier generators
                if not MATRIXSS_IsPointInOrbit(oldSchreierTree, point) or 
                   not generator in ssInfo[level].oldSGS then
                    
                    # Map relations to current free group
                    relations := [];
                    for element in ssInfo[level].relations do
                        gens1 := GeneratorsOfGroup(element[2]);
                        gens2 := GeneratorsOfGroup(ssInfo[level].freeGroup);
                        
                        MATRIXSS_DebugPrint(3, ["Mapping ", element[1],
                                " from ", element[2], " to ",
                                ssInfo[level].freeGroup]);
                        if Length(gens1) <= Length(gens2) then
                            AddSet(relations, 
                                   MappedWord(element[1], gens1, 
                                           gens2{[1 .. Length(gens1)]}));
                        else
                            MATRIXSS_DebugPrint(2, ["3 : Gens1 : ", gens1]);
                            MATRIXSS_DebugPrint(2, ["3 : Gens2 : ", gens2]);
                            
                        fi;
                    od;
                    
                    # Express subgroup generators as words in generators
                    # of free group
                    subgroupGens := [];
                    for element in List(ssInfo[level].partialSGS, i -> i[1]) do
                        element := 
                          PreImagesRepresentative(freeGroupHomo, element);
                        MATRIXSS_DebugPrint(6, ["Adding ", element,
                                " to subgroup gens"]);
                        AddSet(subgroupGens, element);
                        MATRIXSS_DebugPrint(6, ["SubGroup gens : ", 
                                subgroupGens]);
                    od;
                    
                    MATRIXSS_DebugPrint(3, ["Running coset enum with gens : ", 
                            GeneratorsOfGroup(ssInfo[level].freeGroup), 
                            " relations ", relations, " subgroup gens ", 
                            subgroupGens]);
                    
                    # Perform (interruptible) Todd-Coxeter coset enumeration
                    cosetTable := 
                      CosetTableFromGensAndRels(
                              GeneratorsOfGroup(ssInfo[level].freeGroup), 
                              relations, subgroupGens :
                              max := 1 + MATRIXSS_GetOrbitSize(
                                      ssInfo[level].schreierTree),
                              silent);
                    
                    # If coset enumeration was successful and index of the
                    # subgroup was equal to our orbit size, then we know
                    # that the subgroup is stabiliser and we can exit
                    if cosetTable <> fail then
                        MATRIXSS_DebugPrint(2, ["Nof cosets: ", 
                                Length(cosetTable)]);
                        MATRIXSS_DebugPrint(2, ["Orbit size: ", 
                                MATRIXSS_GetOrbitSize(
                                        ssInfo[level].schreierTree)]);
                        if Length(cosetTable) = MATRIXSS_GetOrbitSize(
                                   ssInfo[level].schreierTree) then
                            ssInfo[level].oldSGS := SGS;
                            return;
                        fi;
                    fi;
                    
                    # Compute Schreier generator g for current level
                    schreierGenerator := 
                      MATRIXSS_GetSchreierGenerator_ToddCoxeter(
                              ssInfo[level].schreierTree,
                              generator, point, action, identity,
                              ssInfo[level].IsIdentity, 
                              ssInfo[level].freeGroup, ssInfo[level].genMap);
                                        
                    MATRIXSS_DebugPrint(6, ["Schreier Generator : ", 
                            schreierGenerator]);
                                        
                    if ssInfo[level].IsIdentity(schreierGenerator[1][1], 
                               identity) then
                        continue;
                    fi;
                    
                    # Check if Schreier generator is in stabiliser at 
                    # the current level
                    # Check if g \in H^(i + 1) = <S^(i + 1)>
                    points := [level + 1 .. Length(ssInfo)];
                    MATRIXSS_DebugPrint(4, ["Sifting element ", 
                            schreierGenerator[1], " on levels ", points]);
                    strip := MATRIXSS_Membership_ToddCoxeter(ssInfo{points},
                                     schreierGenerator[1], 
                                     identity, ssInfo[level].freeGroup);
                    MATRIXSS_DebugPrint(6, ["Got sift: ", strip]);
                    word := ShallowCopy(strip[1][2]);
                    strip := [strip[1][1], strip[2]];
                    
                    MATRIXSS_DebugPrint(6, ["Word : ", word]);
                    MATRIXSS_DebugPrint(6, ["Sift : ", strip]);
                    
                    # The drop-out level is in range
                    # [1 .. Length(ssInfo) + 1 - level]
                    # but we want the range given by points
                    strip[2] := strip[2] + level;
                    
                    MakeImmutable(strip);
                    MakeImmutable(word);
                    
                    MATRIXSS_DebugPrint(4, ["Dropout level : ", strip[2]]);
                    
                    if strip[1][1] <> identity then
                        MATRIXSS_DebugPrint(4, ["Residue found"]);
                        
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
                        if strip[2] > Length(ssInfo) then
                            MATRIXSS_ExtendBase(ssInfo, strip[1], identity);
                        fi;
                        
                        # Update partial SGS at each level
                        for recursiveLevel in [level .. strip[2] - 1] do
                            MATRIXSS_DebugPrint(8, ["Adding ",
                                    strip[1][1], " and ", strip[1][2],
                                    " to generators"]);
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   strip[1]);
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   newInverseGenerator);
                        od;
                    fi;
                    
                    for recursiveLevel in [level .. strip[2] - 1] do
                        freeGroup := 
                          FreeGroup(Length(ssInfo[recursiveLevel].
                                  partialSGS));
                        MATRIXSS_DebugPrint(6, ["Generators : ",
                                List(ssInfo[recursiveLevel].
                                     partialSGS, i -> i[1])]);
                        
                        # Important to use GroupWithGenerators, since we want
                        # all our matrices as generators to create the 
                        # generator mapping, even if some are redundant as
                        # generators (ie some are inverses of each other)
                        levelGroup := 
                          GroupWithGenerators(List(ssInfo[recursiveLevel].
                                  partialSGS, i -> i[1]), identity);
                        
                        # Recompute generator mapping, since we have added a
                        # generator
                        gens1 := ShallowCopy(GeneratorsOfGroup(levelGroup));
                        gens2 := ShallowCopy(GeneratorsOfGroup(freeGroup));
                        SortParallel(gens1, gens2);
                        ssInfo[recursiveLevel + 1].genMap := 
                          Immutable([gens1, gens2]);
                        ssInfo[recursiveLevel + 1].freeGroup := freeGroup;
                        
                        # Retrieve generators so that we can map all words to
                        # the same free group
                        gens1 := GeneratorsOfGroup(word[3]);
                        gens2 := GeneratorsOfGroup(ssInfo[recursiveLevel + 1].
                                         freeGroup);
                        
                        if Length(gens1) >= Length(gens2) then
                            MATRIXSS_DebugPrint(6, ["Looking up ", strip[1][2],
                                    " in ", ssInfo[recursiveLevel + 1].
                                    genMap[1]]);
                            
                            # The identity is not among our generators
                            if strip[1][1] <> identity then
                                element := 
                                  ssInfo[recursiveLevel + 1].genMap[2][
                                          Position(ssInfo[recursiveLevel + 1].
                                                  genMap[1], strip[1][2])];
                                
                                # Map the words to the generators in the same
                                # free group                            
                                relation := 
                                  schreierGenerator[2][1] * word[2] *
                                  MappedWord(element, gens2, 
                                          gens1{[1 .. Length(gens2)]});
                            else
                                relation := schreierGenerator[2][1] * word[2];
                            fi;
                            
                            MATRIXSS_DebugPrint(3, ["Adding relation : ",
                                    relation]);
                            
                            # Save relation along with the its free group
                            # so that we can map its generators later
                            Add(ssInfo[recursiveLevel + 1].relations,
                                Immutable([relation, word[3]]));
                        else
                            MATRIXSS_DebugPrint(2, ["2 : Gens1 : ", gens1]);
                            MATRIXSS_DebugPrint(2, ["2 : Gens2 : ", gens2]);
                        fi;
                    od;
                    
                    if strip[1][1] <> identity then
                        # We must now recompute all levels downward from the
                        # dropout level
                        for recursiveLevel in 
                          Reversed([level + 1 .. strip[2]]) do
                            oldSGS := ssInfo[level].oldSGS;
                            ssInfo[level].oldSGS := SGS;
                            MATRIXSS_SchreierToddCoxeterSims(ssInfo, 
                                    partialSGS,
                                    recursiveLevel, identity);
                            ssInfo[level].oldSGS := oldSGS;
                        od;
                    fi;
                fi;
            od;
        od;
        
        ssInfo[level].oldSGS := SGS;
    end;

###############################################################################
#E
