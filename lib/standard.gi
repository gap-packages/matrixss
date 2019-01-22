###############################################################################
#1
#W    standard.gi  The Matrix Schreier Sims package 
#W                 Standard Schreier-Sims implementation
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-07-01 
##
## These are the special routines for the deterministic version of 
## Schreier-Sims algorithm.
###############################################################################

###############################################################################
##
#A StabChainMatrixGroup(G)
##
## An implementation of the Schreier-Sims algorithm, for matrix groups,
## probabilistic version. See "StabChainMatrixGroup!general" for general information
## about the attribute.
##
###############################################################################
InstallMethod(StabChainMatrixGroup, [IsMatrixGroup and IsFinite], 1, 
        function(G)
    local ssInfo, list, generators, level, points, element, SchreierSims, ret,
          identity;
        
###############################################################################
##
#F SchreierSims(ssInfo, partialSGS, level, identity)
##
## The main Schreier-Sims function, which is called for each level.
## \beginitems    
## `ssInfo' & main information structure for the current Schreier-Sims run
##    
## `partialSGS' & given partial strong generating set
##    
## `level' & the level of the call to Schreier-Sims
##
## `identity' & the group identity
## \enditems
##
###############################################################################
    SchreierSims := function(ssInfo, partialSGS, level, identity)
        local generator, point, orbit, strip, schreierGenerator, element, 
              action, recursiveLevel, schreierTree, SGS, oldSGS, points, 
              newPoint, oldSchreierTree, newBasePoint, oldOrbit,
              newSGS, residue, dropoutLevel, newGen;
        
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
        
        MATRIXSS_DebugPrint(3, ["Saved SGS that fixes first ", level - 1, 
                " points ", SGS]);
        
        MATRIXSS_DebugPrint(4, ["Base point : ", ssInfo[level].partialBase]);
        MATRIXSS_DebugPrint(9, ["Hash func : ", ssInfo[level].hash]);
        
        # Compute schreier tree for current level
        # Compute \Delta_i = \beta_i^(H^i) = \beta_i^(<S_i>)
        oldSchreierTree := ssInfo[level].schreierTree;
        
        ssInfo[level].schreierTree := 
          MATRIXSS_GetSchreierTree(ssInfo[level].schreierTree,
                  ssInfo[level].partialBase, SGS, 
                  ssInfo[level].oldSGS,
                  action, ssInfo[level].hash, identity);
        
        orbit := Immutable(MATRIXSS_GetOrbit(ssInfo[level].schreierTree.Tree));
        
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
                
                MATRIXSS_DebugPrint(6, ["Consider", point, 
                        " and ", generators]);
                
                # Avoid rechecking Schreier generators
                if not MATRIXSS_IsPointInOrbit(oldSchreierTree.Tree, point) or 
                   not generator in ssInfo[level].oldSGS then
                                        
                    MATRIXSS_DebugPrint(6, ["Get Schreier Generator"]);
                    
                    # Compute Schreier generator g for current level
                    schreierGenerator := 
                      MATRIXSS_GetSchreierGenerator(
                              ssInfo[level].schreierTree,
                              generator[1], point, action, identity);
                                        
                    MATRIXSS_DebugPrint(3, ["Schreier Generator : ", 
                            schreierGenerator]);
                                        
                    if schreierGenerator = identity then
                        continue;
                    fi;
                    
                    # Check if Schreier generator is in stabiliser at 
                    # the current level
                    # Check if g \in H^(i + 1) = <S^(i + 1)>
                    points := [level + 1 .. Length(ssInfo)];
                    strip := MATRIXSS_Membership(ssInfo{points},
                                     schreierGenerator, 
                                     identity);
                    
                    # The drop-out level is in range 
                    # [1 .. Length(ssInfo) + 1 - level] 
                    # but we want the range given by points
                    dropoutLevel := strip[2] + level;
                    residue := strip[1];
                    
                    MATRIXSS_DebugPrint(5, ["Dropout level : ", dropoutLevel]);
                    
                    if residue <> identity then
                        MATRIXSS_DebugPrint(3, ["Residue found", residue]);
                        
                        # We have found a Schreier generator which is not in
                        # the stabiliser of the current level, and so we must
                        # add the residue of this generator to our partial SGS
                        # in order to make it into a real SGS
                        
                        # Add residue to partial SGS
                        # This makes some levels incomplete and so we must
                        # recompute them recursively
                        newGen := 
                          Immutable([residue, Inverse(residue)]); 
                                  #rec(mat := Inverse(residue.mat),
                                  #    slp := InverseOfStraightLineProgram(
                                  #            residue.slp))]);
                        AddSet(partialSGS, newGen);
                        AddSet(partialSGS, 
                               Immutable(Reversed(ShallowCopy(newGen))));
                        
                        # Possibly extend the base if the Schreier generator
                        # fixes all points in our base
                        if dropoutLevel > Length(ssInfo) then
                            MATRIXSS_ExtendBase(ssInfo, newGen[1], 
                                    identity);
                        fi;
                        
                        # Update partial SGS at each level
                        for recursiveLevel in [level .. dropoutLevel - 1] do
                            if ssInfo[recursiveLevel].action(
                                       ssInfo[recursiveLevel].partialBase,
                                       newGen[1]) = 
                               ssInfo[recursiveLevel].partialBase then
                                AddSet(ssInfo[recursiveLevel].partialSGS, 
                                       newGen);
                                AddSet(partialSGS, 
                                       Immutable(Reversed(ShallowCopy(newGen))));
                            else
                                break;
                            fi;
                        od;

                        # We must not recompute all levels downward from the
                        # dropout level
                        for recursiveLevel in Reversed([level + 1 .. 
                                dropoutLevel]) do
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
    

    ### MAIN Schreier-Sims 
    
    if IsBoundGlobal("MATRIXSS_PROFILE") then
        ProfileFunctions([SchreierSims]);
    fi;
    
    # Get initial set of generators, to be extended to a partial SGS
    generators := GeneratorsOfGroup(G);
        
    # The vector space on which the group acts
    points := FullRowSpace(FieldOfMatrixGroup(G), DimensionOfMatrixGroup(G));
    
    MATRIXSS_DebugPrint(3, ["Group generators : ", generators]);
    
    #idSLP := StraightLineProgramNC([[1, 1, 1, -1]], Length(generators));
    
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
        SchreierSims(ssInfo, generators, level, identity);
    od;
    
    MATRIXSS_DebugPrint(2, ["Matrix Schreier-Sims done"]);
    
    return Immutable(rec(SchreierStructure := ssInfo, SGS := generators));
end);

###############################################################################
#E
