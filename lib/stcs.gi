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
#F MATRIXSS_SchreierToddCoxeterSims(ssInfo, partialSGS, level, identity, cosetFactor)
##
## The main function for the Schreier-Todd-Coxeter-Sims algorithm. It is very
## similar to ordinary Schreier-Sims algorithm and has a similar interface.
## \beginitems    
## `ssInfo' & main information structure for the current Schreier-Sims run
##    
## `partialSGS' & given partial strong generating set
##    
## `level' & the level of the call to Schreier-Sims
##
## `identity' & the group identity
##
## `cosetFactor' & the quotient of the maximum number of cosets generated 
##                  during coset enumeration and the corresponding orbit size
## \enditems
##
###############################################################################
MATRIXSS_SchreierToddCoxeterSims := 
  function(ssInfo, partialSGS, level, identity, cosetFactor)
    local generator, point, orbit, strip, schreierGenerator, element, 
          action, recursiveLevel, schreierTree, SGS, oldSGS, points, 
          newPoint, oldSchreierTree, newBasePoint, oldOrbit,
          newInverseGenerator, newSGS, freeGroup, cosetTable, word,
          subgroupGens, relation, gens1, gens2, gens3, relations, 
          levelGroup, freeGroupHomo, gen, relHomo, ret, dropoutLevel, residue,
          newGen, wordGen, oldFreeGroup, labels, res, slp, baseFreeGens, i,
          gens, freeGens;
    
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
    
    # Compute schreier tree for current level
    # Compute \Delta_i = \beta_i^(H^i) = \beta_i^(<S_i>)
    oldSchreierTree := ssInfo[level].schreierTree;
    
    ssInfo[level].schreierTree := 
      MATRIXSS_GetSchreierTree(ssInfo[level].schreierTree,
              ssInfo[level].partialBase, SGS, ssInfo[level].oldSGS,
              action, ssInfo[level].hash, identity);
    
    #labels := ssInfo[level].schreierTree.Labels;
    
    # Compute free group for current level
    #ssInfo[level].freeGroup := FreeGroup(Length(SGS));
    #MATRIXSS_DebugPrint(4, ["Generators : ", SGS, labels, 
    #        List(SGS, i -> i[1])]);
    #levelGroup := GroupWithGenerators(List(SGS, i -> i[1]), identity);
    
    # Save mapping between generators of current level group and free group
    #gens1 := ShallowCopy(GeneratorsOfGroup(levelGroup));
    #gens2 := ShallowCopy(GeneratorsOfGroup(ssInfo[level].freeGroup));
    #SortParallel(gens1, gens2);
    #ssInfo[level].genMap := Immutable(rec(Generators := gens1, 
    #                                FreeGenerators := gens2));
    
    # Compute homomorphism for use in coset enumeration
    #freeGroupHomo := GroupHomomorphismByImagesNC(ssInfo[level].freeGroup, 
    #                         levelGroup, gens2, gens1);
    #ssInfo[level].word2elt := freeGroupHomo;
    
    MATRIXSS_DebugPrint(3, ["Saved SGS that fixes first ", level - 1, 
            " points ", SGS]);
    
    MATRIXSS_DebugPrint(4, ["Base point : ", ssInfo[level].partialBase]);
    MATRIXSS_DebugPrint(9, ["Hash func : ", ssInfo[level].hash]);
        
    orbit := Immutable(MATRIXSS_GetOrbit(ssInfo[level].schreierTree.Tree));
    
    MATRIXSS_DebugPrint(4, ["Orbit size for level ", level, " is ", 
            Length(orbit)]); 
    
    # We now want to make sure that SGS also fixes the current level
    
    for point in orbit do
        for generator in SGS do
            
            # Avoid rechecking Schreier generators
            if not MATRIXSS_IsPointInOrbit(oldSchreierTree.Tree, point) or 
               not generator in ssInfo[level].oldSGS then
                
                # Map relations to current free group
                #relations := [];
                gens := List(SGS, i -> i[1]);
                MATRIXSS_DebugPrint(4, ["Get gens as words : ", gens]);

                slp := SLPOfElms(gens);
                MATRIXSS_DebugPrint(4, ["Gens : ", Length(gens),
                        " orig gens : ", NrInputsOfStraightLineProgram(slp)]);
                
                freeGroup := FreeGroup(Length(gens) + 
                                     NrInputsOfStraightLineProgram(slp));
                freeGens := ShallowCopy(GeneratorsOfGroup(freeGroup));
                #baseFreeGens := ssInfo[level].genMap.FreeGenerators{[1 .. 
                #                        NrInputsOfStraightLineProgram(slp)]};
                MATRIXSS_DebugPrint(4, ["Gens : ", freeGens{[1 .. 
                        NrInputsOfStraightLineProgram(slp)]},
                        " slp : ", slp]);
                
                relations := 
                  ResultOfStraightLineProgram(slp, 
                          freeGens{[1 .. NrInputsOfStraightLineProgram(slp)]});
                
                for i in [1 .. Length(relations)] do
                    relations[i] := 
                      relations[i] * Inverse(freeGens[
                              NrInputsOfStraightLineProgram(slp) + i]);
                od;
                
                #gens2 := ssInfo[level].genMap.FreeGenerators;
                
                #for element in ssInfo[level].relations do
                #    gens1 := GeneratorsOfGroup(element[2]);
                    
                #    MATRIXSS_DebugPrint(4, ["Mapping ", element[1],
                #            " from ", element[2], " to ",
                #            ssInfo[level].freeGroup]);
                #    if Length(gens1) <= Length(gens2) then
                #        AddSet(relations, 
                #               MappedWord(element[1], gens1, 
                #                       gens2{[1 .. Length(gens1)]}));
                #    fi;
                #od;
                
                # Express subgroup generators as words in generators
                # of free group
                gens := List(ssInfo[level].partialSGS, i -> i[1]);
                
                MATRIXSS_DebugPrint(4, ["Get subgroup gens as words", gens]);
                if Length(gens) > 0 then
                    slp := SLPOfElms(gens);
                    subgroupGens := 
                      ResultOfStraightLineProgram(slp, 
                              freeGens{[1 .. 
                                      NrInputsOfStraightLineProgram(slp)]});
                else
                    subgroupGens := [];
                fi;
                #for element in List(ssInfo[level].partialSGS, i -> i[1]) do
                #    points := [level + 1 .. Length(ssInfo)];
                    
                #    MATRIXSS_DebugPrint(6, ["Get word for ", element,
                #            " using levels ", points]);
                #    strip := MATRIXSS_Membership_ToddCoxeter(ssInfo{points},
                #                     element, 
                #                     identity, ssInfo[level].freeGroup);
                    #MATRIXSS_DebugPrint(6, ["Got sift: ", strip]);
                #    element := ShallowCopy(strip[1][2]);
                #    res := strip[1][1];
                #    Assert(1, res = identity);
                    
                    #element  := word;
                    #element := 
                    #  PreImagesRepresentative(freeGroupHomo, element);
               #     MATRIXSS_DebugPrint(6, ["Adding ", element[1],
               #             " to subgroup gens"]);
               #     AddSet(subgroupGens, element[1]);
               #     MATRIXSS_DebugPrint(6, ["SubGroup gens : ", 
               #             subgroupGens]);
               # od;
                
                MATRIXSS_DebugPrint(2, ["Running coset enum with gens : ", 
                        freeGens, " relations ", relations, " subgroup gens ", 
                        subgroupGens]);
                
                # Perform (interruptible) Todd-Coxeter coset enumeration
                cosetTable := 
                  CosetTableFromGensAndRels(freeGens, AsSet(relations), 
                          AsSet(subgroupGens) :
                          max := 1 + Int(cosetFactor * 
                                  MATRIXSS_GetOrbitSize(
                                          ssInfo[level].schreierTree.Tree)),
                          silent);
                
                # If coset enumeration was successful and index of the
                # subgroup was equal to our orbit size, then we know
                # that the subgroup is stabiliser and we can exit
                if cosetTable <> fail then
                    MATRIXSS_DebugPrint(2, ["Nof cosets: ", 
                            Length(cosetTable)]);
                    MATRIXSS_DebugPrint(2, ["Orbit size: ", 
                            MATRIXSS_GetOrbitSize(
                                    ssInfo[level].schreierTree.Tree)]);
                    if Length(cosetTable) = MATRIXSS_GetOrbitSize(
                               ssInfo[level].schreierTree.Tree) then
                        ssInfo[level].oldSGS := SGS;
                        return;
                    fi;
                fi;
                
                #MATRIXSS_DebugPrint(6, ["Gen map : ", ssInfo[level].genMap,
                #        " generator ", generator]);
                
                # Compute Schreier generator g for current level
                schreierGenerator := 
                  MATRIXSS_GetSchreierGenerator(
                          ssInfo[level].schreierTree.Tree,
                          generator[1], point, action, identity);
                          #ssInfo[level].freeGroup, ssInfo[level].genMap,
                          #ssInfo[level].word2elt);
                
                MATRIXSS_DebugPrint(6, ["Schreier Generator : ", 
                        schreierGenerator]);
                Assert(1, IsObjWithMemory(schreierGenerator));
                
                if schreierGenerator = identity then
                    continue;
                fi;
                
                # Check if Schreier generator is in stabiliser at 
                # the current level
                # Check if g \in H^(i + 1) = <S^(i + 1)>
                points := [level + 1 .. Length(ssInfo)];
                MATRIXSS_DebugPrint(4, ["Sifting element ", 
                        schreierGenerator, " on levels ", points]);
                strip := MATRIXSS_Membership(ssInfo{points},
                                 schreierGenerator, 
                                 identity);
                MATRIXSS_DebugPrint(6, ["Got sift: ", strip]);
                #word := ShallowCopy(strip[1]);
                residue := strip[1];
                Assert(1, IsObjWithMemory(residue));
                
                #oldFreeGroup := word[2];
                newGen := Immutable([residue, Inverse(residue)]);
                Assert(1, IsObjWithMemory(newGen[1]) and 
                       IsObjWithMemory(newGen[2]));
                #wordGen := Immutable([word[1], Inverse(word[1])]);
                
                # The drop-out level is in range
                # [1 .. Length(ssInfo) + 1 - level]
                # but we want the range given by points
                dropoutLevel := strip[2] + level;
                
                #MATRIXSS_DebugPrint(6, ["Word : ", word]);
                
                #MakeImmutable(residue);
                #MakeImmutable(word);
                
                MATRIXSS_DebugPrint(4, ["Dropout level : ", dropoutLevel]);
                
                if newGen[1] <> identity then
                    MATRIXSS_DebugPrint(4, ["Residue found"]);
                    
                    # We have found a Schreier generator which is not in
                    # the stabiliser of the current level, and so we must
                    # add the residue of this generator to our partial SGS
                    # in order to make it into a real SGS
                    
                    #newInverseGenerator := Immutable(Reversed(residue));
                    
                    # Add residue to partial SGS
                    # This makes some levels incomplete and so we must
                    # recompute them recursively
                    AddSet(partialSGS, newGen);
                    AddSet(partialSGS, 
                           Immutable(Reversed(ShallowCopy(newGen))));
                    #Assert(1, IsEvenInt(Length(partialSGS)));
                    
                    # Possibly extend the base if the Schreier generator
                    # fixes all points in our base
                    if dropoutLevel > Length(ssInfo) then
                        MATRIXSS_ExtendBase(ssInfo, newGen[1], identity);
                    fi;
                    
                    # Update partial SGS at each level
                    for recursiveLevel in [level .. dropoutLevel - 1] do
                        if ssInfo[recursiveLevel].action(
                                   ssInfo[recursiveLevel].partialBase,
                                   newGen[1]) = 
                           ssInfo[recursiveLevel].partialBase then
                            MATRIXSS_DebugPrint(8, ["Adding ",
                                    newGen[1], " to generators"]);
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   newGen);
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   Immutable(Reversed(ShallowCopy(newGen))));
                            
                            #Assert(1, IsEvenInt(Length(ssInfo[recursiveLevel].
                            #        partialSGS)));
                        else
                            break;
                        fi;
                    od;
                fi;
                
                #for recursiveLevel in [level .. dropoutLevel - 1] do
                    #freeGroup := 
                    #  FreeGroup(Length(ssInfo[recursiveLevel].partialSGS));
                    
                    # Important to use GroupWithGenerators, since we want
                    # all our matrices as generators to create the 
                    # generator mapping, even if some are redundant as
                    # generators (ie some are inverses of each other)
                    #levelGroup := 
                    #  GroupWithGenerators(List(ssInfo[recursiveLevel].
                    #          partialSGS, i -> i[1]), identity);
                    #MATRIXSS_DebugPrint(6, ["Generators : ",
                    #        GeneratorsOfGroup(levelGroup)]);
                    
                    # Recompute generator mapping, since we have added a
                    # generator
                    #gens1 := ShallowCopy(GeneratorsOfGroup(levelGroup));
                    #Assert(1, IsEvenInt(Length(gens1)));
                    #gens2 := ShallowCopy(GeneratorsOfGroup(freeGroup));
                    #SortParallel(gens1, gens2);
                    #ssInfo[recursiveLevel + 1].genMap := 
                    #  Immutable(rec(Generators := gens1, 
                    #          FreeGenerators := gens2));
                    #ssInfo[recursiveLevel + 1].freeGroup := freeGroup;
                    
                    # Retrieve generators so that we can map all words to
                    # the same free group
                    #gens1 := GeneratorsOfGroup(oldFreeGroup);
                    
                    #MATRIXSS_DebugPrint(6, ["Generators : ", gens1, gens2]);
                    
                    #if Length(gens2) >= Length(gens1) then
                        
                        # The identity is not among our generators
                    #    if newGen[1] <> identity then
                    #        MATRIXSS_DebugPrint(6, ["Looking up ", newGen[2],
                    #                " in ", ssInfo[recursiveLevel + 1].
                    #                genMap.Generators]);
                    #        element := 
                    #          ssInfo[recursiveLevel + 1].genMap.FreeGenerators[
                    #                  Position(ssInfo[recursiveLevel + 1].
                    #                          genMap.Generators, 
                    #                          newGen[2])];
                            
                            # Map the words to the generators in the same
                            # free group                            
                    #        relation := 
                    #          schreierGenerator[2][1] * wordGen[2] *
                    #          MappedWord(element, gens2, 
                    #                  gens1{[1 .. Length(gens2)]});
                    #    else
                    #        relation := schreierGenerator[2][1] * wordGen[2];
                    #    fi;
                        
                    #    MATRIXSS_DebugPrint(3, ["Adding relation : ",
                    #            Immutable([relation, oldFreeGroup])]);
                    #    MATRIXSS_DebugPrint(3, ["Relations: ",
                    #            ssInfo[recursiveLevel + 1].relations]);
                        
                        # Save relation along with the its free group
                        # so that we can map its generators later
                    #    AddSet(ssInfo[recursiveLevel + 1].relations,
                    #        Immutable([relation, oldFreeGroup]));
                    #fi;
                #od;
                
                if newGen[1] <> identity then
                    # We must now recompute all levels downward from the
                    # dropout level
                    for recursiveLevel in 
                      Reversed([level + 1 .. dropoutLevel]) do
                        oldSGS := ssInfo[level].oldSGS;
                        ssInfo[level].oldSGS := SGS;
                        MATRIXSS_SchreierToddCoxeterSims(ssInfo, 
                                partialSGS,
                                recursiveLevel, identity, cosetFactor);
                        ssInfo[level].oldSGS := oldSGS;
                    od;
                fi;
            fi;
        od;
    od;
    
    ssInfo[level].oldSGS := SGS;
end;

InstallGlobalFunction(SchreierToddCoxeterSims, function(G)
    local ssInfo, list, generators, level, points, element, ret, identity;
    
    # Get initial set of generators, to be extended to a partial SGS
    generators := GeneratorsOfGroup(G);
        
    # The vector space on which the group acts
    points := FullRowSpace(FieldOfMatrixGroup(G), DimensionOfMatrixGroup(G));
    
    MATRIXSS_DebugPrint(3, ["Group generators : ", generators]);
    
    # Compute initial partial SGS and base and fill ssInfo
    ret := MATRIXSS_GetPartialBaseSGS(generators, points);
    generators := ret[1];
    ssInfo     := ret[2];
    identity   := ret[3];
    
    MATRIXSS_DebugPrint(3, ["Partial sgs : ", generators]);
    MATRIXSS_DebugPrint(3, ["Initial base length : ", Length(ssInfo)]);
    MATRIXSS_DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    
    # Call Schreier-Sims algorithm for each level (starting from top)
    for level in Reversed([1 .. Length(ssInfo)]) do
        MATRIXSS_SchreierToddCoxeterSims(ssInfo, generators, level, 
                identity, 6/5);
    od;
    
    MATRIXSS_DebugPrint(2, ["Matrix Schreier-Sims done"]);
    MATRIXSS_DebugPrint(2, ["Order is : ", MatrixGroupOrderStabChain(ssInfo)]);
    
    return Immutable(rec(SchreierStructure := ssInfo, SGS := generators));
end);

###############################################################################
#E
