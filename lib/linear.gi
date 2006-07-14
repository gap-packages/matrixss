###############################################################################
#1
#W    linear.gi  The Matrix Schreier Sims package 
#W               Nearly linear-time implementation
##
#H    File      : $RCSfile$
#H    Author    : Henrik B‰‰rnhielm
#H    Dev start : 2004-08-08 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
##
## These are the special routines for the nearly linear-time version, described
## in Babai et al, 1991.
###############################################################################

Revision.("matrixss/lib/linear_gi") := 
  "@(#)$Id$";

###############################################################################
##
#A StabChainMatrixGroup(G)
##
## An implementation of the Schreier-Sims algorithm, for matrix groups.
## This version is inspired by the nearly linear time algorithm, described in 
## \cite{seress03}. See "StabChainMatrixGroup!general" for general information
## about the attribute.
##
###############################################################################
InstallMethod(StabChainMatrixGroup, [IsMatrixGroup and IsFinite], 3, 
        function(G)
    local ssInfo, list, generators, level, points, element, ConstructSGS,
          CompletePointStabiliserSubgroup, ret, iters;
    
###############################################################################
##
#F CompletePointStabiliserSubgroup(ssInfo, element, level, identity, maxIdentitySifts)
##
## The work-horse of the nearly linear-time algorithm, called for each level.
## \beginitems    
## `ssInfo' & main information structure for the current Schreier-Sims run
##    
## `element' & the element which
##
## `partialSGS' & given partial strong generating set
##    
## `identity' & the group identity
## \enditems
##
###############################################################################
    CompletePointStabiliserSubgroup := 
      function(ssInfo, element, level, identity, maxIdentitySifts)
        local orbit, point, elements, recursiveLevel, iterations,
              nIdentitySifts, schreierGen, ret, residue, dropoutLevel,
              subsetSize, randomPointSet, nonTrivialResidue;
        
        MATRIXSS_DebugPrint(2, ["CompletePointStabiliserSubgroup at level ", 
                level]);
        
        elements := [element];
        for recursiveLevel in [level .. Length(ssInfo)] do
            UniteSet(elements, 
                    ssInfo[recursiveLevel].schreierTree.Labels);
        od;
        
        MakeImmutable(elements);
        MATRIXSS_DebugPrint(2, ["Transversals : ", elements]);
        
        # Recalculate Schreier trees if necessary
        orbit := MATRIXSS_GetOrbit(ssInfo[level].schreierTree.Tree);
        for point in orbit do            
            if not MATRIXSS_IsPointInOrbit(ssInfo[level].schreierTree.Tree,
                       ssInfo[level].action(point, element[1])) then
                
                MATRIXSS_DebugPrint(2, ["Rebuilding Schreier tree"]);
                ssInfo[level].schreierTree := 
                  MATRIXSS_GetSchreierTree(ssInfo[level].schreierTree,
                          ssInfo[level].partialBase, elements, 
                          ssInfo[level].oldSGS,
                          ssInfo[level].action, ssInfo[level].hash, identity :
                          ExtendSchreierTree, ShallowSchreierTree);
                MATRIXSS_DebugPrint(2, ["Schreier tree height ", 
                        ssInfo[level].schreierTree.Height]);
                ssInfo[level].oldSGS := elements;
                break;
            fi;
        od;
        
        iterations := 1;
        for recursiveLevel in [level .. Length(ssInfo)] do
            iterations := iterations + ssInfo[recursiveLevel].
                          schreierTree.Height;
        od;
        iterations := iterations * 64;
        nIdentitySifts := 0;
        
        MATRIXSS_DebugPrint(2, ["Nof iterations : ", iterations]);
        MATRIXSS_DebugPrint(2, ["Nof outer iters : ", 
                maxIdentitySifts * LogInt(Length(ssInfo), 2)]);
        
        # Repeat until we reach the required number of consecutive inner 
        # iterations with identity sifts only
        repeat
            nonTrivialResidue := false;
            iters := iterations;
            
            # repeat the above calculated number of iterations
            repeat
                MATRIXSS_DebugPrint(4, ["Nof iterations left : ", iters]);
                schreierGen := 
                  MATRIXSS_RandomSchreierGenerator(
                          ssInfo[level].schreierTree.Tree, elements, 
                          ssInfo[level].action, identity);
                MATRIXSS_DebugPrint(8, ["Random Schreier gen: ", schreierGen]);
                
                if not ssInfo[level].IsIdentity(schreierGen, identity) then
                    ret := MATRIXSS_Membership(ssInfo{[level .. 
                                   Length(ssInfo)]},
                                   schreierGen, identity);
                    residue := ret[1];
                    dropoutLevel := ret[2] + level - 1;
                    
                    if residue <> identity then
                        nonTrivialResidue := true;
                        
                        MATRIXSS_DebugPrint(2, ["Non-trivial residue! ", 
                                residue]);
                        MATRIXSS_DebugPrint(2, ["Dropout level ", 
                                dropoutLevel]);
                        
                        if dropoutLevel > Length(ssInfo) then
                            MATRIXSS_ExtendBase(ssInfo, residue, identity);
                        fi;
                        
                        for recursiveLevel in Reversed([level .. 
                                dropoutLevel]) do
                            CompletePointStabiliserSubgroup(ssInfo, residue,
                                    recursiveLevel, identity, 
                                    maxIdentitySifts);
                        od;
                    fi;
                fi;
                
                iters := iters - 1;
                MATRIXSS_DebugPrint(4, ["Nof iterations left : ", iters]);
            until iters <= 0;
            
            # Check if we only got sifts to identity
            if not nonTrivialResidue then
                nIdentitySifts := nIdentitySifts + 1;
            else
                nIdentitySifts := 0;
            fi;
            
            MATRIXSS_DebugPrint(2, ["Nof iters without residues : ",
                    nIdentitySifts]);
        until nIdentitySifts >= maxIdentitySifts * LogInt(Length(ssInfo), 2);
        
        MATRIXSS_DebugPrint(2, ["CompletePointStabiliserSubgroup done"]);
    end;
            
###############################################################################
##
#F ConstructSGS(ssInfo, partialSGS, identity)
##
## The main Schreier-Sims function for the nearly linear-time algorithm.
## \beginitems    
## `ssInfo' & main information structure for the current Schreier-Sims run
##    
## `partialSGS' & given partial strong generating set
##    
## `identity' & the group identity
## \enditems
##
###############################################################################
    ConstructSGS := function(ssInfo, partialSGS, identity)
        local level, element, ret, residue, dropoutLevel;
                
        MATRIXSS_DebugPrint(2, ["First run"]);
        
        for element in partialSGS do
            ret := MATRIXSS_Membership(ssInfo, element[1], identity);
            residue := ret[1];
            dropoutLevel := ret[2];
            
            if residue <> identity then
                if dropoutLevel > Length(ssInfo) then
                    MATRIXSS_ExtendBase(ssInfo, residue, identity);
                fi;
                for level in Reversed([1 .. dropoutLevel]) do
                    CompletePointStabiliserSubgroup(ssInfo, residue, level, 
                            identity, 4);
                od;
            fi;
        od;
        
        MATRIXSS_DebugPrint(2, ["Second run"]);
        
        for level in Reversed([1 .. Length(ssInfo)]) do
            # According to Seress, the value 4 is enough
            CompletePointStabiliserSubgroup(ssInfo, 
                    Immutable([identity, identity]), level, identity, 4);
        od;
        
        MATRIXSS_DebugPrint(2, ["Building SGS"]);
        
        for level in ssInfo do
            level.partialSGS := level.schreierTree.Labels;
        od;
    end;
    

    ### MAIN Schreier-Sims 
    
    if IsBoundGlobal("MATRIXSS_PROFILE") then
        ProfileFunctions([ConstructSGS, CompletePointStabiliserSubgroup]);
    fi;
    
    # Check if we want to use the random Schreier-Sims
    if ValueOption("Linear") = fail then
        TryNextMethod();
    fi;
    
    # Get initial set of generators, to be extended to a partial SGS
    generators := GeneratorsOfGroup(G);
        
    # The vector space on which the group acts
    points := FullRowSpace(FieldOfMatrixGroup(G), DimensionOfMatrixGroup(G));
    
    MATRIXSS_DebugPrint(3, ["Group generators : ", generators]);
    
    # Compute initial partial SGS and base and fill ssInfo
    ret := MATRIXSS_GetPartialBaseSGS(generators, Identity(G), points);
    generators := ret[1];
    ssInfo     := ret[2];
    
    MATRIXSS_DebugPrint(3, ["Partial sgs : ", generators]);
    MATRIXSS_DebugPrint(3, ["Initial base length : ", Length(ssInfo)]);
    MATRIXSS_DebugPrint(3, ["Calling recursive Schreier-Sims"]);
    
    ConstructSGS(ssInfo, generators, Identity(G));
    
    MATRIXSS_DebugPrint(2, ["Matrix Schreier-Sims done"]);
    
    return Immutable(rec(SchreierStructure := ssInfo, SGS := generators));
end);

###############################################################################
#E
