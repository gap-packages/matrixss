###############################################################################
##
#W    random.gi  The Matrix Schreier Sims package 
#W               Standard Schreier-Sims implementation
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

Revision.("matrixss/lib/random_gi") := 
  "@(#)$Id$";

# An implementation of the Schreier-Sims algorithm, for matrix groups
InstallGlobalFunction(MatrixRandomSchreierSims, function(G, p)
    local ssInfo, list, generators, level, points, element, RandomSchreierSims,
          identitySifts, UpdateSchreierTrees, ComputeOrder;
    
    # Updates the given Schreier trees w.r.t. to the given partial SGS
    UpdateSchreierTrees := function(ssInfo, partialSGS, identity)
        local SGS, level;
        
        for level in [1 .. Length(ssInfo)] do
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
    
    # Computes the order of the group defined by the given Schreier trees
    ComputeOrder := function(ssInfo)
        local order, levelStruct;
        
        order := 1;
        for levelStruct in ssInfo do
            order := order * MATRIXSS_GetOrbitSize(levelStruct.schreierTree);
        od;
        
        return order;
    end;
    
    # The main random Schreier-Sims function
    # ssInfo - main information structure for the Schreier-Sims 
    # partialSGS - given partial strong generating set
    # maxIdentitySifts - maximum number of consecutive elements that sifts to 
    #                    identity before the algorithm terminates
    # identity - the group identity
    # low_order - lower bound on the group order (must be >= 1)
    # high_order - upper bound on the group order, or 0 if not available
    RandomSchreierSims := 
      function(ssInfo, partialSGS, maxIdentitySifts, identity, 
              low_order, high_order)
      local nIdentitySifts, element, strip, newInverseGenerator, level, order;
        
        # Check if we are already done
        order := ComputeOrder(ssInfo);
        if high_order > 0 and order >= high_order then
            return;
        fi;
        
        nIdentitySifts := 0;
        
        # Sanity check
        Assert(1, low_order >= 1 and high_order >= 0 and 
               (low_order <= high_order or high_order = 0));
        
        # Loop until our order meets the lower bound and we have sifted the
        # given number of consecutive random elements to identity
        while (nIdentitySifts <= maxIdentitySifts or
               order < low_order) do
            
            # Get a random element, from a hopefully uniform distribution
            element := PseudoRandom(Group(List(partialSGS, i -> i[1])));
            
            # Our functions expect the elements to be vectors with the element
            # and its inverse
            element := [element, Inverse(element)];
            MATRIXSS_DebugPrint(8, ["Random element to sift : ", element]);
            
            # Sift the random element
            strip := MATRIXSS_Membership(ssInfo, element, identity);
                             
            if strip[1][1] <> identity then
                
                # Add residue to our partial SGS
                newInverseGenerator := Immutable(Reversed(strip[1]));
                AddSet(partialSGS, strip[1]);
                AddSet(partialSGS, newInverseGenerator);
                
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
                UpdateSchreierTrees(ssInfo, partialSGS, identity);
                order := ComputeOrder(ssInfo);
                
                # Check if we are done
                MATRIXSS_DebugPrint(4, ["Order is : ", order]);
                if high_order > 0 and order >= high_order then
                    return;
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

    if not IsMatrixGroup(G) then
        Error("<G> must be a matrix group");
    fi;
    
    if not (p > 0 and p < 1) then
        Error("<p> must be a probability, 0 < p < 1");
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
            Order(G), Order(G));
    
    MATRIXSS_DebugPrint(2, ["Random matrix Schreier-Sims done"]);
    
    # Create output structure
    list := [[], generators, []];
    for level in [1 .. Length(ssInfo)] do
        Add(list[1], ssInfo[level].partialBase);
        Add(list[3], ssInfo[level].schreierTree);
    od;
    
    return Immutable(list);
end);
