###############################################################################
#1
#W    random.gi  The Matrix Schreier Sims package 
#W               Standard Schreier-Sims implementation
##
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-07-01 
##
## These are the special routines for the probabilistic implementation of 
## Schreier-Sims algorithm.
##
###############################################################################

###############################################################################
##
#A StabChainMatrixGroup(G)
##
## An implementation of the Schreier-Sims algorithm, for matrix groups,
## probabilistic version. See "StabChainMatrixGroup!general" for general information
## about the attribute.
##
## In addition to the general Options of the attribute `StabChainMatrixGroup',
## the probabilistic algorithm is aware of the following:
## \beginitems
## `Probability' & (lower bound for) probability of correct solution, which
##                    defaults to 3/4
##
## `Verify' & Boolean parameter which signifies if the base and SGS computed
##          using the random Schreir-Sims algorithm should be verified
##          using the Schreier-Todd-Coxeter-Sims algorithm.
##          Defaults to `false.'
##
## `OrderLowerBound' & Lower bound for the order of `G', must be $\geq$ 1.
##                   Defaults to 1.
##
## `OrderUpperBound' & Upper bound for the order of `G', or 0 if unknown.
##                   Defaults to 0.
##
## Note that if the order of `G' is known, so that 
## `OrderLowerBound' = `OrderUpperBound' = `Size(G)'
## then the randomized algorithm always produces a correct base and SGS, so
## there is no need of verification. Also, the verification will extend the
## given base and SGS to a complete base and SGS if needed.
##
###############################################################################
InstallMethod(StabChainMatrixGroup, [IsMatrixGroup and IsFinite], 2, 
        function(G)
    local ssInfo, list, generators, level, points, element, RandomSchreierSims,
          identitySifts, UpdateSchreierTrees, lowOrder, highOrder, p, verify,
          cosetFactor, ret, SchreierSims, RandomElement, residue, dropoutLevel,
          identity;
    
    # Updates the given Schreier trees w.r.t. to the given partial SGS
    UpdateSchreierTrees := function(ssInfo, dropoutLevel, partialSGS, identity)
        local SGS, level, ret;
        
        for level in [1 .. dropoutLevel] do
            if level > 1 then
                SGS := ShallowCopy(ssInfo[level - 1].partialSGS);
            else
                SGS := ShallowCopy(partialSGS);
            fi;
            MakeImmutable(SGS);
            
            if SGS <> ssInfo[level].oldSGS then
                ssInfo[level].schreierTree := 
                  MATRIXSS_GetSchreierTree(ssInfo[level].schreierTree,
                          ssInfo[level].partialBase, SGS, 
                          ssInfo[level].oldSGS,
                          ssInfo[level].action, ssInfo[level].hash, identity);
                ssInfo[level].oldSGS := SGS;
            fi;
        od;
    end;
        
###############################################################################
##
#F RandomSchreierSims(ssInfo, partialSGS, maxIdentitySifts, identity, low_order, high_order)
##    
## The main random Schreier-Sims function.
## \beginitems    
## `ssInfo' & main information structure for the Schreier-Sims 
##    
## `partialSGS' & given partial strong generating set
##
## `maxIdentitySifts' & maximum number of consecutive elements that sifts to 
##                    identity before the algorithm terminates
##    
## `identity' & the group identity
##    
## `lowOrder' & lower bound on the group order (must be $\geq$ 1)
##    
## `highOrder' & upper bound on the group order, or 0 if not available
## \enditems
##
###############################################################################
    RandomSchreierSims := 
      function(ssInfo, partialSGS, maxIdentitySifts, identity, 
              lowOrder, highOrder)
      local nIdentitySifts, element, strip, newInverseGenerator, level, order,
            point, sgsGroup, newGen, invGen, gens;
        
        MATRIXSS_DebugPrint(2, ["Update Schreier trees "]);
        UpdateSchreierTrees(ssInfo, Length(ssInfo), partialSGS, identity);
        MATRIXSS_DebugPrint(2, ["Done updating Schreier trees "]);
        
        # Check if we are already done
        if highOrder > 0 or lowOrder > 1 then
           order := MatrixGroupOrderStabChain(ssInfo);
           MATRIXSS_DebugPrint(3, ["Order is : ", order]);
            if order >= highOrder then
                return;
            fi;
        else
            order := 1;
        fi;
        
        nIdentitySifts := 0;
        
        # Sanity check
        Assert(1, lowOrder >= 1 and highOrder >= 0 and 
               (lowOrder <= highOrder or highOrder = 0));
        
        gens := List(partialSGS, i -> i[1]);
        MATRIXSS_DebugPrint(4, ["Generators1 : ", gens]);
        sgsGroup := Group(gens);
        
        # Loop until our order meets the lower bound and we have sifted the
        # given number of consecutive random elements to identity
        while (nIdentitySifts < maxIdentitySifts or
               order < lowOrder) do
            
            # Get a random element, from a hopefully uniform distribution
            element := PseudoRandom(sgsGroup);
            
            # Our functions expect the elements to be vectors with the element
            # and its inverse
            MATRIXSS_DebugPrint(4, ["Random element to sift : ", element]);
            
            # Sift the random element
            strip := MATRIXSS_Membership(ssInfo, element, identity);
            dropoutLevel := strip[2];
            residue := strip[1];
                             
            if residue <> identity then
                
                # Add residue to our partial SGS
                newGen := Immutable([residue, Inverse(residue)]);
                invGen := Immutable(Reversed(ShallowCopy(newGen)));
                AddSet(partialSGS, newGen);
                AddSet(partialSGS, invGen);
                
                MATRIXSS_DebugPrint(4, ["New gens : ", newGen, invGen]);
                
                gens := List(partialSGS, i -> i[1]);
                MATRIXSS_DebugPrint(4, ["Generators2 : ", gens]);
                sgsGroup := Group(gens);
                
                MATRIXSS_DebugPrint(3, ["Dropout level : ", dropoutLevel]);
                                
                # Extend base if needed
                if dropoutLevel > Length(ssInfo) then
                    MATRIXSS_ExtendBase(ssInfo, newGen[1], identity);
                fi;
                                
                # Update partial SGS at each level
                for level in [1 .. dropoutLevel - 1] do
                    if ssInfo[level].action(ssInfo[level].partialBase,
                               newGen[1]) = ssInfo[level].partialBase then
                        AddSet(ssInfo[level].partialSGS, newGen);
                        AddSet(ssInfo[level].partialSGS, invGen);
                    fi;
                od;
                
                # Recompute Schreier trees
                UpdateSchreierTrees(ssInfo, dropoutLevel, partialSGS, 
                        identity);
                if highOrder > 0 or lowOrder > 1 then
                    order := MatrixGroupOrderStabChain(ssInfo);
                
                    # Check if we are done
                    MATRIXSS_DebugPrint(3, ["Order is : ", order]);
                    if order >= highOrder then
                        return;
                    fi;
                fi;
                
                nIdentitySifts := 0;
            else
                nIdentitySifts := nIdentitySifts + 1;
                
                MATRIXSS_DebugPrint(4, ["Sift to identity! Number : ", 
                        nIdentitySifts]);
            fi;
        od;
    end;
    
    RandomElement := function(ssInfo, identity)
        local levelStruct, element;
        
        element := identity;
        for levelStruct in ssInfo do
            element := element * MATRIXSS_RandomCosetRepresentative(
                               ssInfo[level].schreierTree.Tree, 
                               ssInfo[level].action, identity);
        od;
        
        return element;
    end;
        
    ### MAIN Schreier-Sims 
    
    # Check if we want to use the random Schreier-Sims
    if ValueOption("Random") = fail then
        TryNextMethod();
    fi;
    
    if IsBoundGlobal("MATRIXSS_PROFILE") then
        ProfileFunctions([UpdateSchreierTrees, RandomSchreierSims]);
    fi;
    
    p := ValueOption("Probability");
    verify := ValueOption("Verify");
    lowOrder := ValueOption("OrderLowerBound");
    highOrder := ValueOption("OrderUpperBound");
    cosetFactor := ValueOption("CosetFactor");
    
    MATRIXSS_DebugPrint(2, ["Prob: ", p]);
    MATRIXSS_DebugPrint(2, ["Verify: ", verify]);
    
    if not IsRat(p) or p = 0 or p >= 1 then
        p := 3/4;
    fi;
    
    if verify = fail or not IsBool(verify) then
        verify := false;
    fi;
    
    if not IsPosInt(lowOrder) then
        lowOrder := 1;
    fi;
    
    if not IsPosInt(highOrder) or 
       highOrder < lowOrder then
        highOrder := 0;
    fi;
    
    if cosetFactor = fail or not IsRat(cosetFactor) then
        cosetFactor := 6/5;
    fi;
    
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
    
    # Calculate number of needed identity sifts to meet the required 
    # probability of correctness
    # Note however, that since the algorithm does not use uniformly distributed
    # random elements, this number must be seen as a lower bound and is not a
    # theoretical guarantee for the required probability
    identitySifts := Maximum([LogInt(DenominatorRat(1 - p), 2) - 
                             LogInt(NumeratorRat(1 - p), 2), 30]);
    
    MATRIXSS_DebugPrint(3, ["Partial sgs : ", generators]);
    MATRIXSS_DebugPrint(3, ["Initial base length : ", Length(ssInfo)]);
    MATRIXSS_DebugPrint(2, ["Calling random Schreier-Sims requiring ", 
            identitySifts, " consecutive identity sifts"]);
    
    Assert(1, identitySifts >= 1);
    
    # Call probabilistic Schreier-Sims algorithm 
    RandomSchreierSims(ssInfo, generators, identitySifts, 
            identity, lowOrder, highOrder);
    
    MATRIXSS_DebugPrint(2, ["Random matrix Schreier-Sims done"]);
    MATRIXSS_DebugPrint(2, ["Order is : ", MatrixGroupOrderStabChain(ssInfo)]);
    
    if verify then
        if ValueOption("STCS") <> fail then 
            MATRIXSS_DebugPrint(2, ["Verifying using STCS"]);
            
            for level in [1 .. Length(ssInfo)] do
                ssInfo[level].oldSGS := AsSet([]); 
                ssInfo[level].schreierTree := 
                  MATRIXSS_CreateInitialSchreierTree(
                          ssInfo[level].partialBase, ssInfo[level].hash, 
                          identity);
            od;
            
        
            # Call Schreier-Sims algorithm for each level (starting from top)
            for level in Reversed([1 .. Length(ssInfo)]) do
                MATRIXSS_SchreierToddCoxeterSims(ssInfo, generators, level, 
                        identity, cosetFactor);
            od;
        else
            repeat
                element := MatrixSchreierSimsVerify(ssInfo, generators, 
                                   Identity(G));
                Assert(2, element.Level = 0, "Verification failed\n");
                if element.Level > 0 then
                    # Continue probabilistic Schreier-Sims...
                    RandomSchreierSims(ssInfo, generators, identitySifts, 
                            identity, lowOrder, highOrder);
                fi;
            until element.Level = 0;
        fi;
    fi;
    
    MATRIXSS_DebugPrint(2, ["Order is : ", MatrixGroupOrderStabChain(ssInfo)]);
    
    return Immutable(rec(SchreierStructure := ssInfo, SGS := generators));
end);

###############################################################################
#E