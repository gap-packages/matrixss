###############################################################################
##
#W    standard.gi  The Matrix Schreier Sims package 
#W                 Standard Schreier-Sims implementation
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
##    Dev start : 2004-07-01 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
###############################################################################

Revision.("matrixss/lib/standard_gi") := 
  "@(#)$Id$";

# An implementation of the Schreier-Sims algorithm, for matrix groups
InstallGlobalFunction(MatrixSchreierSims, function(G)
    local ssInfo, list, generators, level, points, element, SchreierSims;
        
    # The main Schreier-Sims function
    # ssInfo - main information structure for the current Schreier-Sims run
    # partialSGS - given partial strong generating set
    # level - the level of the call to Schreier-Sims
    # identity - the group identity
    SchreierSims := function(ssInfo, partialSGS, level, identity)
        local generator, point, orbit, strip, schreierGenerator, element, 
              action, recursiveLevel, schreierTree, SGS, oldSGS, points, 
              newPoint, oldSchreierTree, newBasePoint, oldOrbit,
              newInverseGenerator, newSGS;
        
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
                                        
                    # Compute Schreier generator g for current level
                    schreierGenerator := 
                      MATRIXSS_GetSchreierGenerator(ssInfo[level].schreierTree,
                              generator, point, action, identity,
                              ssInfo[level].IsIdentity);
                                        
                    MATRIXSS_DebugPrint(6, ["Schreier Generator : ", 
                            schreierGenerator]);
                                        
                    if ssInfo[level].IsIdentity(schreierGenerator[1], 
                               identity) then
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
                    strip[2] := strip[2] + level;
                    
                    MATRIXSS_DebugPrint(5, ["Dropout level : ", strip[2]]);
                    
                    if strip[1][1] <> identity then
                        MATRIXSS_DebugPrint(3, ["Residue found"]);
                        
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
                        
                        # Update partial SGS at each level
                        for recursiveLevel in [level .. strip[2] - 1] do
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   strip[1]);
                            AddSet(ssInfo[recursiveLevel].partialSGS, 
                                   newInverseGenerator);
                        od;
                        
                        # Possibly extend the base if the Schreier generator
                        # fixes all points in our base
                        if strip[2] > Length(ssInfo) then
                            MATRIXSS_ExtendBase(ssInfo, strip[1], identity);
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
        SchreierSims(ssInfo, generators, level, Identity(G));
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
