###############################################################################
#1
#W    random.gi  The Matrix Schreier Sims package 
#W               Standard Schreier-Sims implementation
##
#H    File      : $RCSfile$
#H    Author    : Henrik B��rnhielm
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
#M StabChainMatrixGroup(G)
##
## An implementation of the Schreier-Sims algorithm, for matrix groups,
## probabilistic version. See "StabChainMatrixGroup!general" for general information
## about the attribute.
##
## In addition to the general Options of the attribute `StabChainMatrixGroup',
## the probabilistic algorithm is aware of the following:
## \beginitems
## Probability & (lower bound for) probability of correct solution, which
##                    defaults to 3/4
##
## Verify & Boolean parameter which signifies if the base and SGS computed
##          using the random Schreir-Sims algorithm should be verified
##          using the Schreier-Todd-Coxeter-Sims algorithm.
##          Defaults to `false.'
##
## OrderLowerBound & Lower bound for the order of `G', must be >= 1.
##                   Defaults to 1.
##
## OrderUpperBound & Upper bound for the order of `G', or 0 if unknown.
##                   Defaults to 0.
##
## Note that if the order of `G' is known, so that 
## `OrderLowerBound' = `OrderUpperBound' = |G|
## then the randomized algorithm always produces a correct base and SGS, so
## there is no need of verification. Also, the verification will extend the
## given base and SGS to a complete base and SGS if needed.
##
###############################################################################
InstallMethod(StabChainMatrixGroup, [IsMatrixGroup and IsFinite], 2, 
        function(G)
    local ssInfo, list, generators, level, points, element, RandomSchreierSims,
          identitySifts, UpdateSchreierTrees, Rattle, InitRattle,
          ScrambleRattle, AddRattleGenerator, low_order, high_order, p, verify;
    
    # Updates the given Schreier trees w.r.t. to the given partial SGS
    UpdateSchreierTrees := function(ssInfo, dropoutLevel, partialSGS, identity)
        local SGS, level;
        
        for level in [1 .. dropoutLevel] do
            if level > 1 then
                SGS := ShallowCopy(ssInfo[level - 1].partialSGS);
            else
                SGS := ShallowCopy(partialSGS);
            fi;
            MakeImmutable(SGS);
            
            if ValueOption("ExtendSchreierTree") <> fail then
                ssInfo[level].schreierTree := 
                  MATRIXSS_ExtendSchreierTree(
                          ssInfo[level].schreierTree, 
                          SGS, ssInfo[level].oldSGS, 
                          ssInfo[level].action, 
                          ssInfo[level].hash);
            else
                ssInfo[level].schreierTree := 
                  MATRIXSS_CreateInitialSchreierTree(
                          ssInfo[level].partialBase,
                          ssInfo[level].hash, identity);
                ssInfo[level].schreierTree :=
                  MATRIXSS_ComputeSchreierTree(
                          ssInfo[level].schreierTree, 
                          SGS, ssInfo[level].action);
            fi;
            
            ssInfo[level].oldSGS := SGS;
        od;
    end;
        
    # Initialise Rattle random element generator
    InitRattle := function(partialSGS, length, identity, nScrambles)
        local i, RattleState;
        
        RattleState := [ShallowCopy(partialSGS), [identity, identity]];
        for i in [Length(partialSGS) + 1 .. length] do
            Add(RattleState[1], Immutable([identity, identity]));
        od;
        
        ScrambleRattle(RattleState, nScrambles);
        return RattleState;
    end;
    
    # Scramble the Rattle state
    ScrambleRattle := function(RattleState, nScrambles)
        while nScrambles > 0 do
            Rattle(RattleState);
            nScrambles := nScrambles - 1;
        od;
    end;
    
    # Generate a random group element using the given Rattle state
    Rattle := function(RattleState)
        local i, j;
        
        i := Random([1 .. Length(RattleState[1])]);
        RattleState[2] := [RattleState[2][1] * 
                           RattleState[1][i][1],
                           RattleState[1][i][2] * 
                           RattleState[2][2]];
        
        i := Random([1 .. Length(RattleState[1])]);
        
        repeat
            j := Random([1 .. Length(RattleState[1])]);
        until i <> j;
        
        if Random([true, false]) then
            RattleState[1][i] := 
              Immutable([RattleState[1][i][1] * RattleState[1][j][1], 
                      RattleState[1][j][2] * RattleState[1][i][2]]);
        else
            RattleState[1][i] := 
              Immutable([RattleState[1][j][1] * RattleState[1][i][1], 
                      RattleState[1][i][2] * RattleState[1][j][2]]);
        fi;
        
        return RattleState[2];
    end;
    
###############################################################################
##
#F RandomSchreierSims(ssInfo, partialSGS, maxIdentitySifts, identity, low_order, high_order, nScrambles, RattleFactor)
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
## `low_order' & lower bound on the group order (must be >= 1)
##    
## `high_order' & upper bound on the group order, or 0 if not available
##
## `nScrambles' & number of initial scrambles of the Rattle pool each time
##              a new pool is created
##    
## `RattleFactor' & the Rattle pool is this many times bigger than the
##                partial SGS
## \enditems
##
###############################################################################
    RandomSchreierSims := 
      function(ssInfo, partialSGS, maxIdentitySifts, identity, 
              low_order, high_order, nScrambles, RattleFactor)
      local nIdentitySifts, element, strip, newInverseGenerator, level, order,
            RattleState;
        
        # Check if we are already done
        if high_order > 0 or low_order > 1 then
           order := MATRIXSS_ComputeOrder(ssInfo);
            if order >= high_order then
                return;
            fi;
        else
            order := 1;
        fi;
        
        nIdentitySifts := 0;
        
        # Sanity check
        Assert(1, low_order >= 1 and high_order >= 0 and 
               (low_order <= high_order or high_order = 0));
        
        RattleState := InitRattle(partialSGS, RattleFactor * 
                               Length(partialSGS), identity, nScrambles);
        
        # Loop until our order meets the lower bound and we have sifted the
        # given number of consecutive random elements to identity
        while (nIdentitySifts <= maxIdentitySifts or
               order < low_order) do
            
            # Get a random element, from a hopefully uniform distribution
            element := Rattle(RattleState);
            
            # Our functions expect the elements to be vectors with the element
            # and its inverse
            MATRIXSS_DebugPrint(8, ["Random element to sift : ", element]);
            
            # Sift the random element
            strip := MATRIXSS_Membership(ssInfo, element, identity);
                             
            if strip[1][1] <> identity then
                
                # Add residue to our partial SGS
                newInverseGenerator := Immutable(Reversed(strip[1]));
                AddSet(partialSGS, strip[1]);
                AddSet(partialSGS, newInverseGenerator);
                
                RattleState := InitRattle(partialSGS, 
                                       RattleFactor * Length(partialSGS),
                                       identity, nScrambles);
                
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
                if high_order > 0 or low_order > 1 then
                    order := MATRIXSS_ComputeOrder(ssInfo);
                
                    # Check if we are done
                    MATRIXSS_DebugPrint(4, ["Order is : ", order]);
                    if order >= high_order then
                        return;
                    fi;
                fi;
                
                nIdentitySifts := 0;
            else
                nIdentitySifts := nIdentitySifts + 1;
                
                MATRIXSS_DebugPrint(6, ["Sift to identity! Number : ", 
                        nIdentitySifts]);
            fi;
        od;
    end;
    
    ### MAIN Schreier-Sims 
    
    # Check if we want to use the random Schreier-Sims
    if ValueOption("Random") = fail then
        TryNextMethod();
    fi;
    
    p := ValueOption("Probability");
    verify := ValueOption("Verify");
    low_order := ValueOption("OrderLowerBound");
    high_order := ValueOption("OrderUpperBound");
    
    if p = fail or not IsRat(p) or p = 0 or p >= 1 then
        p := 3/4;
    fi;
    
    if verify = fail or not IsBool(verify) then
        verify := false;
    fi;
    
    if low_order = fail or not IsPosInt(low_order) then
        low_order := 1;
    fi;
    
    if high_order = fail or not IsPosInt(high_order) or 
       high_order < low_order then
        high_order := 0;
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
    
    MATRIXSS_DebugPrint(3, ["Group generators : ", generators]);
    
    # Compute initial partial SGS and base and fill ssInfo
    generators := MATRIXSS_GetPartialBaseSGS(generators, ssInfo, Identity(G), 
                          points);
    
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
    
    # Call Schreier-Sims algorithm for each level (starting from top)
    RandomSchreierSims(ssInfo, generators, identitySifts, Identity(G), 
            low_order, high_order, 1000, 10);
    
    MATRIXSS_DebugPrint(2, ["Random matrix Schreier-Sims done"]);
    
    if verify then
        MATRIXSS_DebugPrint(2, ["Verifying using STCS"]);
        
        # Call Schreier-Sims algorithm for each level (starting from top)
        for level in Reversed([1 .. Length(ssInfo)]) do
            MATRIXSS_SchreierToddCoxeterSims(ssInfo, generators, level, 
                    Identity(G));
        od;
    fi;
    
    return Immutable(ssInfo);
end);

###############################################################################
#E