###############################################################################
##
#W    stcs.gi  The Matrix Schreier Sims package 
#W             Schreier-Todd-Coxeter-Sims implementation
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
##    Dev start : 2004-07-05 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/lib/stcs_gi") := 
  "@(#)$Id$";

# An implementation of the Schreier-Sims algorithm, for matrix groups
InstallGlobalFunction(MatrixSchreierToddCoxeterSims, function(G)
    local ssInfo, list, generators, level, points, element, 
          SchreierToddCoxeterSims;
        
    # The main Schreier-Sims function
    # ssInfo - main information structure for the current Schreier-Sims run
    # partialSGS - given partial strong generating set
    # level - the level of the call to Schreier-Sims
    # identity - the group identity
    SchreierToddCoxeterSims := function(ssInfo, partialSGS, level, identity)
        local generator, point, orbit, strip, schreierGenerator, element, 
              action, recursiveLevel, schreierTree, SGS, oldSGS, points, 
              newPoint, oldSchreierTree, newBasePoint, oldOrbit,
              newInverseGenerator, newSGS, freeGroup, cosetTable, word,
              groupGens, subgroupGens, subgroupRels, tc1, tc2, tc3, tc4,
              relation, gens1, gens2, gens3, relations, homo, G, H,
              sgens, fgens, fsgens, grels, T;
        
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
        freeGroup := FreeGroup(Length(SGS));
        MATRIXSS_DebugPrint(2, ["Generators : ", List(SGS, i -> i[1])]);
        tc1 := Group(List(SGS, i -> i[1]), identity);
        ssInfo[level].freeGroupHomo := 
          GroupHomomorphismByImages(freeGroup, tc1,
                  GeneratorsOfGroup(freeGroup), GeneratorsOfGroup(tc1));
        
        MATRIXSS_DebugPrint(2, ["Free group : ", freeGroup]);
        MATRIXSS_DebugPrint(2, ["Free group homo : ", 
                ssInfo[level].freeGroupHomo]);
        
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
                    
                    if level < Length(ssInfo) then
                        freeGroup := FreeGroup(Length(SGS));
                        MATRIXSS_DebugPrint(2, ["Generators : ", 
                                List(SGS, i -> i[1])]);
                        tc1 := Group(List(SGS, i -> i[1]), identity);
                        ssInfo[level].freeGroupHomo := 
                          GroupHomomorphismByImages(freeGroup, tc1,
                                  GeneratorsOfGroup(freeGroup), 
                                  GeneratorsOfGroup(tc1));
                        
                        groupGens := [];
                        for element in List(SGS, i -> i[1]) do
                            MATRIXSS_DebugPrint(2, ["Adding ", 
                                    PreImagesRepresentative(ssInfo[level].
                                            freeGroupHomo, element), 
                                    " to group gens"]);
                            element := PreImagesRepresentative(ssInfo[level].
                                               freeGroupHomo, element);
                            if not Inverse(element) in groupGens then
                                AddSet(groupGens, element);
                            fi;
                            MATRIXSS_DebugPrint(2, ["Group gens : ", groupGens]);
                        od;
                        
                        relations := [];
                        for element in ssInfo[level].relations do
                            AddSet(relations,
                                MappedWord(element[1], 
                                        GeneratorsOfGroup(element[2]),
                                        GeneratorsOfGroup(freeGroup)));
                        od;
                        
                        tc3 := freeGroup / relations;
                        homo := GroupHomomorphismByImages(tc3, tc1,
                                        GeneratorsOfGroup(tc3),
                                        GeneratorsOfGroup(tc1));
                                        
                        subgroupGens := [];
                        for element in List(ssInfo[level].partialSGS, 
                                i -> i[1]) do
                            element := PreImagesRepresentative(homo, element);
                            MATRIXSS_DebugPrint(2, ["Adding ", element,
                                    " to subgroup gens"]);
                            if not Inverse(element) in subgroupGens then
                                AddSet(subgroupGens, element);
                            fi;
                            MATRIXSS_DebugPrint(2, ["SubGroup gens : ", 
                                    groupGens]);
                        od;
                                               
                        MATRIXSS_DebugPrint(2, ["Group1 : ", tc3]);
                        MATRIXSS_DebugPrint(2, ["Type1 : ", TypeObj(tc3)]);
                        MATRIXSS_DebugPrint(2, ["Type2 : ", 
                                TypeObj(subgroupGens)]);
                        tc2 := Subgroup(tc3, subgroupGens);
                        
                        H := tc2;
                        
                        # Get whole group <G> of <H>.
                        G := FamilyObj( H )!.wholeGroup;
                        
                        # get some variables
                        fgens := FreeGeneratorsOfFpGroup( G );
                        grels := RelatorsOfFpGroup( G );
                        sgens := GeneratorsOfGroup( H );
                        fsgens := List( sgens, gen -> UnderlyingElement( gen ) );
                        
                        # Construct the coset table of <G> by <H>.
                        T := CosetTableFromGensAndRels( fgens, grels, fsgens :
                                     
                                     max := Int(11/10 * MATRIXSS_GetOrbitSize(
                                             ssInfo[level].schreierTree)),
                                     silent);
                        cosetTable := T;
                        
                        #CosetTableDefaultMaxLimit := 
                        #  Int(11/10 * MATRIXSS_GetOrbitSize(
                        #          ssInfo[level].schreierTree));
                        #cosetTable := 
                        #  CosetTable(tc3, tc2 : 
                        #          max := 11/10 * MATRIXSS_GetOrbitSize(
                        #                  ssInfo[level].schreierTree),
                        #          silent := true);
                        
                        MATRIXSS_DebugPrint(2, ["Running coset enum with gens : ", 
                                groupGens, " relations ", 
                                relations,
                                " subgroup gens ", subgroupGens]);
                                                
                        #cosetTable := 
                        #  CosetTableFromGensAndRels(
                        #          groupGens,
                        #          relations,
                        #          subgroupGens : 
                        #          max := 11/10 * 
                        #          MATRIXSS_GetOrbitSize(
                        #                  ssInfo[level].schreierTree),
                        #          silent);
                        
                        MATRIXSS_DebugPrint(2, ["Coset table : ", cosetTable]);
                        if cosetTable <> fail and Length(cosetTable) = 
                           MATRIXSS_GetOrbitSize(ssInfo[level].schreierTree) then
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
                              ssInfo[level].freeGroupHomo);
                                        
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
                    strip := MATRIXSS_Membership_ToddCoxeter(ssInfo{points},
                                     schreierGenerator[1], 
                                     identity);
                    MATRIXSS_DebugPrint(2, ["Got sift: ", strip]);
                    word := ShallowCopy(strip[1][2]);
                    strip := [strip[1][1], strip[2]];
                    
                    MATRIXSS_DebugPrint(2, ["Word : ", word]);
                    MATRIXSS_DebugPrint(2, ["Sift : ", strip]);
                    
                    # The drop-out level is in range 
                    # [1 .. Length(ssInfo) + 1 - level] 
                    # but we want the range given by points
                    strip[2] := strip[2] + level;
                    
                    MATRIXSS_DebugPrint(2, ["Dropout level : ", strip[2]]);
                    
                    if strip[1][1] <> identity then
                        MATRIXSS_DebugPrint(2, ["Residue found"]);
                        
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
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   strip[1]);
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   newInverseGenerator);
                            
                            # Recompute free group homomorphism and add new
                            # relation, which involves residue
                            freeGroup := 
                              FreeGroup(Length(ssInfo[recursiveLevel].
                                      partialSGS));
                            MATRIXSS_DebugPrint(6, ["Generators : " ,
                                    List(ssInfo[recursiveLevel].
                                           partialSGS, i -> i[1])]);
                            tc1 := Group(List(ssInfo[recursiveLevel].
                                           partialSGS, i -> i[1]));
                            MATRIXSS_DebugPrint(6, ["Free group : " ,
                                    freeGroup]);
                            MATRIXSS_DebugPrint(6, ["Group : " , tc1]);
                            MATRIXSS_DebugPrint(6, ["Gens1 : " , 
                                    GeneratorsOfGroup(freeGroup)]);
                            MATRIXSS_DebugPrint(6, ["Gens2 : " , 
                                    GeneratorsOfGroup(tc1)]);

                            ssInfo[recursiveLevel + 1].freeGroupHomo := 
                              GroupHomomorphismByImages(freeGroup, tc1,
                                      GeneratorsOfGroup(freeGroup),
                                      GeneratorsOfGroup(tc1));
                            
                            MATRIXSS_DebugPrint(2, ["Schreier gen : " , 
                                    schreierGenerator[2][1]]);
                            MATRIXSS_DebugPrint(6, ["Sift : " , 
                                    strip[1][2]]);
                            MATRIXSS_DebugPrint(6, ["Type1 : ",
                                    TypeObj(schreierGenerator[2][1])]);
                            MATRIXSS_DebugPrint(6, ["Type2 : ",
                                    TypeObj(word[2])]);
                            MATRIXSS_DebugPrint(6, ["Type3 : ",
                                    TypeObj(PreImagesRepresentative(
                                            ssInfo[recursiveLevel].
                                            freeGroupHomo, strip[1][2]))]);
                            MATRIXSS_DebugPrint(2, ["Word : ", word[2]]);
                            
                            gens1 := GeneratorsOfGroup(
                                             PreImages(ssInfo[level + 
                                                     1].freeGroupHomo));
                            gens2 := GeneratorsOfGroup(
                                             PreImages(ssInfo[recursiveLevel + 1].
                                                     freeGroupHomo));
                            gens3 := GeneratorsOfGroup(
                                             PreImages(ssInfo[level].
                                                     freeGroupHomo));
                            
                            MATRIXSS_DebugPrint(2, ["Gens1 : ", gens1]);
                            MATRIXSS_DebugPrint(2, ["Gens2 : ", gens2]);
                            MATRIXSS_DebugPrint(2, ["Gens3 : ", gens3]);
                            
                            if Length(gens2) >= Length(gens3) and
                               Length(gens2) >= Length(gens1) then
                                
                                tc4 := MappedWord(schreierGenerator[2][1],
                                               gens3,
                                               gens2{[1 .. Length(gens3)]});
                                
                                if Length(ExtRepOfObj(word[2])) > 0 then
                                    tc3 := 
                                      MappedWord(word[2], 
                                              gens1,
                                              gens2{[1 .. Length(gens1)]});
                                else
                                    tc3 := Identity(freeGroup);
                                fi;
                                
                                tc2 := PreImagesRepresentative(
                                               ssInfo[recursiveLevel + 1].
                                               freeGroupHomo, strip[1][2]);

                                MATRIXSS_DebugPrint(2, ["Word1 : ", tc3]);
                                MATRIXSS_DebugPrint(2, ["Word2 : ", tc4]);
                                MATRIXSS_DebugPrint(2, ["Word3 : ", tc2]);
                                
                                relation := tc4 * tc3 * tc2;
                            
                                MATRIXSS_DebugPrint(2, ["Adding relation : ",
                                        relation]);
                                AddSet(ssInfo[recursiveLevel + 1].relations,
                                       [relation, freeGroup]);
                            fi;
                            #Add(ssInfo[recursiveLevel + 1].relations, 
                            #    ObjByExtRep(FamilyObj(tc2),
                            #            schreierGenerator[2][1]) *
                            #    ObjByExtRep(FamilyObj(tc2), word[2]) * tc2);
                            #Add(ssInfo[recursiveLevel + 1].relations,
                            #    schreierGenerator[2][1] * word[2] * 
                            #    PreImagesRepresentative(
                            #            ssInfo[recursiveLevel + 1].
                            #            freeGroupHomo, strip[1][2]));
                        od;
                                                
                        # We must now recompute all levels downward from the
                        # dropout level
                        for recursiveLevel in Reversed([level + 1 .. 
                                strip[2]]) do
                            oldSGS := ssInfo[level].oldSGS;
                            ssInfo[level].oldSGS := SGS;
                            SchreierToddCoxeterSims(ssInfo, partialSGS,
                                    recursiveLevel, identity);
                            ssInfo[level].oldSGS := oldSGS;
                        od;
                    else
                        # Add new relation (without residue) to all affected
                        # levels
                        #for recursiveLevel in [level .. strip[2] - 1] do
                        #    AddSet(ssInfo[recursiveLevel].relations, 
                        #           schreierGenerator[2][1] * word[2]);
                        #od;
                    fi;
                    #for recursiveLevel in [level .. strip[2] - 1] do
                    #    Add(ssInfo[recursiveLevel].relations, 
                    #        PreImagesRepresentative(
                    #                ssInfo[level].freeGroupHomo,
                    #                schreierGenerator[1][1] * strip[1][2]));
                    #od;
                    
                fi;
            od;
        od;
        
        ssInfo[level].oldSGS := SGS;
    end;
    

    ### MAIN Schreier-Sims 

    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    # Get initial set of generators, to be extended to a partial SGS
    generators := GeneratorsOfGroup(G);
    
    if ValueOption("CleverBasePoints") <> fail then
        # Get a list of possibly good base points for this group
        MATRIXSS_BasePointStore := BasisVectorsForMatrixAction(G);
    fi;
    
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
    
    
    MATRIXSS_DebugPrint(3, ["Group generators : ", generators]);
    
    # Compute initial partial SGS and base and fill ssInfo
    generators := MATRIXSS_GetPartialBaseSGS(generators, ssInfo, Identity(G), 
                          points);
    
    MATRIXSS_DebugPrint(3, ["Partial sgs : ", generators]);
    MATRIXSS_DebugPrint(3, ["Initial base length : ", Length(ssInfo)]);
    MATRIXSS_DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    
    # Call Schreier-Sims algorithm for each level (starting from top)
    for level in Reversed([1 .. Length(ssInfo)]) do
        SchreierToddCoxeterSims(ssInfo, generators, level, Identity(G));
    od;
    
    MATRIXSS_DebugPrint(2, ["Matrix Schreier-Sims done"]);
    
    # Create output structure
    list := [[], generators, []];
    for level in [1 .. Length(ssInfo)] do
        Add(list[1], ssInfo[level].partialBase);
        Add(list[3], ssInfo[level].schreierTree);
    od;
    
    return Immutable(list);
end);

#E
