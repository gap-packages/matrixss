###############################################################################
#1
#W    random.gi  The Matrix Schreier Sims package 
#W               Standard Schreier-Sims implementation
##
#H    File      : $RCSfile$
#H    Author    : Henrik Bäärnhielm
#H    Dev start : 2004-07-01 
##
#H    Version   : $Revision$
#H    Date      : $Date$
#H    Last edit : $Author$
##
#H    @(#)$Id$
##
## These are the special routines for the probabilistic implementation of 
## Schreier-Sims algorithm.
##
###############################################################################

Revision.("matrixss/lib/random_gi") := 
  "@(#)$Id$";

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
          cosetFactor, ret, SchreierSims, RandomElement;
    
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
            point, sgsGroup;
        
        UpdateSchreierTrees(ssInfo, Length(ssInfo), partialSGS, identity);
        
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
        
        sgsGroup := Group(List(partialSGS, i -> i[1]), identity);
        
        # Loop until our order meets the lower bound and we have sifted the
        # given number of consecutive random elements to identity
        while (nIdentitySifts < maxIdentitySifts or
               order < lowOrder) do
            
            # Get a random element, from a hopefully uniform distribution
            element := PseudoRandom(sgsGroup);
            element := [element, Inverse(element)];
            
            # Our functions expect the elements to be vectors with the element
            # and its inverse
            MATRIXSS_DebugPrint(4, ["Random element to sift : ", element]);
            
            # Sift the random element
            strip := MATRIXSS_Membership(ssInfo, element, identity);
                             
            if strip[1][1] <> identity then
                
                # Add residue to our partial SGS
                newInverseGenerator := Immutable(Reversed(strip[1]));
                AddSet(partialSGS, strip[1]);
                AddSet(partialSGS, newInverseGenerator);
                sgsGroup := Group(List(partialSGS, i -> i[1]), identity);
                
                MATRIXSS_DebugPrint(3, ["Dropout level : ", strip[2]]);
                
                # Update partial SGS at each level
                for level in [1 .. strip[2] - 1] do
                    AddSet(ssInfo[level].partialSGS, strip[1]);
                    AddSet(ssInfo[level].partialSGS, newInverseGenerator);
                od;
                
                # Extend base if needed
                if strip[2] > Length(ssInfo) then
                    MATRIXSS_ExtendBase(ssInfo, strip[1], identity);
                fi;
                                
                # Recompute Schreier trees
                UpdateSchreierTrees(ssInfo, strip[2], partialSGS, 
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
    ret := MATRIXSS_GetPartialBaseSGS(generators, Identity(G), points);
    generators := ret[1];
    ssInfo     := ret[2];
    
    # Calculate number of needed identity sifts to meet the required 
    # probability of correctness
    # Note however, that since the algorithm does not use uniformly distributed
    # random elements, this number must be seen as a lower bound and is not a
    # theoretical guarantee for the required probability
    identitySifts := LogInt(DenominatorRat(1 - p), 2) - 
                     LogInt(NumeratorRat(1 - p), 2); 
    
    MATRIXSS_DebugPrint(3, ["Partial sgs : ", generators]);
    MATRIXSS_DebugPrint(3, ["Initial base length : ", Length(ssInfo)]);
    MATRIXSS_DebugPrint(2, ["Calling random Schreier-Sims requiring ", 
            identitySifts, " consecutive identity sifts"]);
    
    Assert(1, identitySifts >= 1);
    
    # Call probabilistic Schreier-Sims algorithm 
    RandomSchreierSims(ssInfo, generators, identitySifts, Identity(G),
            lowOrder, highOrder);
    
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
                          Identity(G));
            od;
            
        
            # Call Schreier-Sims algorithm for each level (starting from top)
            for level in Reversed([1 .. Length(ssInfo)]) do
                MATRIXSS_SchreierToddCoxeterSims(ssInfo, generators, level, 
                        Identity(G), cosetFactor);
            od;
        else
            repeat
                element := MatrixSchreierSimsVerify(ssInfo, generators, 
                                   Identity(G));
                Assert(2, element.Level = 0, "Verification failed\n");
                if element.Level > 0 then
                    # Continue probabilistic Schreier-Sims...
                    RandomSchreierSims(ssInfo, generators, identitySifts, 
                            Identity(G), lowOrder, highOrder);
                fi;
            until element.Level = 0;
        fi;
    fi;
    
    MATRIXSS_DebugPrint(2, ["Order is : ", MatrixGroupOrderStabChain(ssInfo)]);
    
    return Immutable(rec(SchreierStructure := ssInfo, SGS := generators));
end);

###############################################################################
#E