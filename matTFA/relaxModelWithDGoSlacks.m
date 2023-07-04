function [model, relaxedDGoVarsValues, modelwDGoSlackVars] = relaxModelWithDGoSlacks(model, minObjSolVal, rxnNameListNoDGoRelax)
% This function takes as input a TFA model and attempts to find a feasible
% solution. if the model is infeasible, then it adds DGo slack variables
% and identifies a solution (from many) that will enable feasibility. As
% outputs it provides the model with the relaxed DGo-bounds that enable
% this feasible solution, as wel as a model that contains the DGo slack
% variables for further analysis. Additionally it provides a cell table
% with all the reaction DGo values that needed to be relaxed to obtain this
% solution and their corresponding bounds before and after relaxation
%
% INPUTS:
%
% - model: TFA model
% - minObjSolVal: minimum lower bound of the objective (that is maximized)
%   that the desired solution should satisfy
% - rxnNameListNoDGoRelax: Reactions whose DeltaG should not undergo relaxation
%
% OUTPUTS:
%
% - model: output TFA model with relaxed DGo values that allow feasibility
% - relaxedDGoVarsValues: cell table with four columns. In the 1st column
%   are all the reaction DGo values that needed to be relaxed to obtain this
%   solution. In the columns 2-3 are their corresponding lower-upper bounds
%   before relaxation, and in columns 4-5 the lower-upper bounds after
%   relaxation.
% - modelwDGoSlackVars: TFA model that contains the DGo slack variables for
%   further analysis
%
% Georgios Fengos 12 April 2017

if ~exist('rxnNameListNoDGoRelax','var') || isempty(rxnNameListNoDGoRelax)
    rxnNameListNoDGoRelax = [];
end

% Initialize the model as the model we have created up to this point:
modelwDGoSlackVars = model;

% Check if the model gives at all a feasible solution.
solTFA = solveTFAmodelCplex(model);

% Initialize the cell array where we will put the variable names that
% need to be relaxed, and their corresponding and values
relaxedDGoVarsValues = [];

if isempty(solTFA.x) || solTFA.val<minObjSolVal
    % If not investigate which are the DGo that have to be relaxed in order to
    % obtain a feasible solution
    
    % Add slack variables for the DGo in each of the G-constraints:
    % - find the index of the G-constraints:
    id_G_constraints = find(~cellfun(@isempty, regexpi(modelwDGoSlackVars.constraintNames, '^(G_)')) == 1);
    if isempty(id_G_constraints)
        warning('There exist no DG-constraints in the model!!')
    end
    
    % Create a model with slack variables that can be used to relax the DGo
    % of each reaction with thermodynamic constraints
    
    % - First exclude those reactions for which we would rather not allow thermodynamic relaxation
    if ~isempty(rxnNameListNoDGoRelax)
        G_ToExclude = cellfun(@(x) ['G_',x], rxnNameListNoDGoRelax, 'UniformOutput', false);
        id_G_ToExclude = find_cell(G_ToExclude, modelwDGoSlackVars.constraintNames(id_G_constraints));
        id_G_constraints(id_G_ToExclude) = [];
    end
    
    slackVarIDs = cellfun(@(x)regexprep(x, '^(G_)', ''),...
        modelwDGoSlackVars.constraintNames(id_G_constraints), 'un', 0);
    % remove NoX tag from reactions (GECKO model
    if any(startsWith(model.rxns, 'arm_')) && any(startsWith(model.mets, 'prot_'))
        slackVarIDsUnique = unique(cellfun(@(x)regexprep(x, 'No\d+$', ''),...
            slackVarIDs, 'un', 0), 'stable');
    else
        slackVarIDsUnique = slackVarIDs;
    end
    for i = 1:numel(slackVarIDsUnique)
        % From each G-constraint get the corresponding
        % - reaction name
        correspRxnName = slackVarIDs{find(startsWith(slackVarIDs, slackVarIDsUnique{i}), 1)};
        % - and reaction index
        id_correspRxn = find_cell(correspRxnName, modelwDGoSlackVars.rxns);
        % Only add the slack variables if:
        % - DGoError of the corresponding reaction is not flagged with large values
        % - the reaction is NOT hydrogen transport
        isHTransport = isequal(model.metFormulas(model.S(:,id_correspRxn(1))>0),'H') && isequal(model.metFormulas(model.S(:,id_correspRxn(1))<0),'H');
        if all(modelwDGoSlackVars.rxnDeltaGRerr(id_correspRxn, 1) < 1E6) && ~isHTransport
            correspRxnName = regexprep(correspRxnName, 'No\d+$', '');
            modelwDGoSlackVars = addNewVariableInTFA(modelwDGoSlackVars, strcat('DGPosSlack_', correspRxnName), 'C', [0 1000]);
            modelwDGoSlackVars = addNewVariableInTFA(modelwDGoSlackVars, strcat('DGNegSlack_', correspRxnName), 'C', [0 1000]);
        end
        row_idx = id_G_constraints(startsWith(slackVarIDs, correspRxnName));
        modelwDGoSlackVars.A(row_idx, [end-1 end]) = ...
            repmat([1 -1], size(row_idx));
    end
    
    
    if ~exist('minObjSolVal','var') || isempty(minObjSolVal)
        minObjSolVal = 1e-6;
    end
    
    idOfObjective = find(modelwDGoSlackVars.f == 1);
    
    if ~isequal(size(idOfObjective),[1 1])
        warning('Attention: The model has multiple objective functions assigned, i.e:')
        modelwDGoSlackVars.varNames(idOfObjective)
        error('Please select only one!!')
    end
    
    % Set the lower bound of the objective to the desired solution value
    modelwDGoSlackVars.var_lb(idOfObjective) = minObjSolVal;
    
    % Get the indices of all GDo slack variables
    DGS_vi = getAllVar(modelwDGoSlackVars,{'DGPosSlack','DGNegSlack'});
    if isempty(DGS_vi)
        error('There exist no DGo slack variables!')
    end
    
    modelwDGoSlackVars.f = zeros(length(modelwDGoSlackVars.f),1);
    
    % We want to minimize the sum of the slack variables that we have to
    % add to the system in order to get feasible solution. Therefore we
    % - set all indices of the DGo slack variables to 1
    modelwDGoSlackVars.f(DGS_vi) = 1;
    % - set the optimization to "minimization" of this sum
    modelwDGoSlackVars.objtype = 1;
    
    solTFAwDGoSlacks = solveTFAmodelCplex(modelwDGoSlackVars);
    if isempty(solTFAwDGoSlacks.val) || isnan(solTFAwDGoSlacks.val) || solTFAwDGoSlacks.val==0
        error('Although we added slack variables for the DGo, we were not able to get a feasible solution for the objective:\n %s with value greater equal than: %d \n',modelwDGoSlackVars.varNames{idOfObjective} , minObjSolVal)
    end
    
    % Here we find one solution, the solution that minimizes the sum of all
    % slack variables. The question that arizes is: what
    % alternative solutions exist? THIS PART REMAINS TO BE EXTENDED!!!
    
    % identify the slack variables that need to be non-zero
    DGoSlackVars = modelwDGoSlackVars.varNames(DGS_vi(solTFAwDGoSlacks.x(DGS_vi) > 1e-9));
    fprintf('The slack variables that need to be non-zero, are:\n')
    fprintf('%s\n', DGoSlackVars{:})
    DGoSlackValues = solTFAwDGoSlacks.x(find_cell(DGoSlackVars, modelwDGoSlackVars.varNames));
    % Go to each of the DG-constraints that need to be relaxed, and
    % find by how many sigmas this constraint need to be relaced
    rxnNamesToBeRelaxed = regexprep(DGoSlackVars, '(^DGNegSlack_)|(^DGPosSlack_)', '');
    if any(startsWith(model.rxns, 'arm_')) && any(startsWith(model.mets, 'prot_'))
        % if we are dealing with an ecModel, search for corresponding
        % reaction duplicates
        tmp_res = cellfun(@(x)regexp(model.rxns, [x 'No\d+$'], 'match'),...
            rxnNamesToBeRelaxed, 'un', 0);
        tmp_res = [tmp_res{:}];
        rxn_ids = [tmp_res{~cellfun(@isempty, tmp_res)}]';
        
%         new_slack_var_array = zeros(size(rxn_ids));
%         for i = 1:numel(rxnNamesToBeRelaxed)
%             idx = startsWith(rxn_ids, rxnNamesToBeRelaxed(i));
%             new_slack_var_array(idx) = DGoSlackValues(i);
%         end
%         
        rxnNamesToBeRelaxed = rxn_ids;
%         DGoSlackValues = new_slack_var_array;
    end
    % Since suffixes for forward and and reverse directions of reactions
    % were removed, there are duplicate reaction IDs.
    % In this case, the DGo error for both directions will be relaxed,
    % which is reasonable, because the relaxation should affect both
    % directions. This was also the case in the original implementation but
    % is an issue now, because GECKO ecModels contain only irreversible
    % reactions.
    id_rxnNamesToBeRelaxed = find_cell(rxnNamesToBeRelaxed, model.rxns);
    rxnNamesToBeRelaxed = model.rxns(id_rxnNamesToBeRelaxed);
    
    model.rxnDeltaGRerr_updated = model.rxnDeltaGRerr;
    for i = 1:numel(id_rxnNamesToBeRelaxed)
        
        % Get the name of the i-th DG-naught variable that needs to be relaxed
        varNameToBeRelaxed = ['DGo_',rxnNamesToBeRelaxed{i}];
        
        % Find the DGo and the DGo-error for this reaction
        DGo_RxnToBeRelaxed_i      = model.rxnDeltaGR(id_rxnNamesToBeRelaxed(i));
        DGoError_RxnToBeRelaxed_i = model.rxnDeltaGRerr(id_rxnNamesToBeRelaxed(i));
        if DGoError_RxnToBeRelaxed_i == 0
            warning('The DGo-error of the reaction %s appears to be zero!!\nTo express the corresponding DGo-relaxation in terms of sigmas, we will assume that the DGo_%s-error is 10%% of the DGo_%s!', rxnNamesToBeRelaxed{i}, rxnNamesToBeRelaxed{i}, rxnNamesToBeRelaxed{i})
            DGoError_RxnToBeRelaxed_i = 0.1*DGo_RxnToBeRelaxed_i;
        end
        
        % Find the slacks names and variables that are associated with this reaction:
        PosNegSlackVars_i = [{['DGNegSlack_',rxnNamesToBeRelaxed{i}]};{['DGPosSlack_',rxnNamesToBeRelaxed{i}]}];
        if any(startsWith(model.rxns, 'arm_')) && any(startsWith(model.mets, 'prot_'))
            PosNegSlackVars_i = cellfun(@(x)regexprep(x, 'No\d+$', ''),...
                PosNegSlackVars_i, 'un', 0);
        end     
        DGoSlackVar_i   = DGoSlackVars(find_cell(PosNegSlackVars_i,DGoSlackVars));
        DGoSlackValue_i = DGoSlackValues(find_cell(PosNegSlackVars_i,DGoSlackVars));
        
        % Express the relaxation in terms of multiples of sigmas (standard deviation "error" of the DGo estimation):
        relax_XSigmas = 1 + ceil(abs(DGoSlackValue_i/DGoError_RxnToBeRelaxed_i));
        % Relax the values of the corresponding DGo:
        %modelwDGoSlackVars.var_lb(find_cell(varNameToBeRelaxed, modelwDGoSlackVars.varNames)) = DGo_RxnToBeRelaxed_i - relax_XSigmas * DGoError_RxnToBeRelaxed_i;
        %modelwDGoSlackVars.var_ub(find_cell(varNameToBeRelaxed, modelwDGoSlackVars.varNames)) = DGo_RxnToBeRelaxed_i + relax_XSigmas * DGoError_RxnToBeRelaxed_i;
        fprintf('Relaxation of %s : %0.3f +- %0.3f (sigma) ---> %0.3f +- %d*%0.3f \n', varNameToBeRelaxed, DGo_RxnToBeRelaxed_i, DGoError_RxnToBeRelaxed_i, DGo_RxnToBeRelaxed_i, relax_XSigmas, DGoError_RxnToBeRelaxed_i)
        refRangeAndRelaxedRange = [DGo_RxnToBeRelaxed_i + DGoError_RxnToBeRelaxed_i*[- 1 1] DGo_RxnToBeRelaxed_i + relax_XSigmas*DGoError_RxnToBeRelaxed_i*[- 1 1]];
        relaxedDGoVarsValues = [relaxedDGoVarsValues; [varNameToBeRelaxed num2cell(refRangeAndRelaxedRange)]];
        % save updated DGo error
        model.rxnDeltaGRerr_updated(id_rxnNamesToBeRelaxed(i)) = relax_XSigmas*DGoError_RxnToBeRelaxed_i;
    end
    
    fprintf('It was not possible to obtain TFA-feasible solutions with the original DGo values calculated by GCM!\n')
    fprintf('We had to relax the following DGo ranges as follows:\n')
    [[{'varName','LB-before','UB-before','LB-relaxed','UB-relaxed'}];relaxedDGoVarsValues]
    
    % If we set the relaxed bounds to the original model
    for i = 1:size(relaxedDGoVarsValues,1)
        tmpIdx = find_cell(relaxedDGoVarsValues(i,1), model.varNames);
        model.var_lb(tmpIdx) = cell2mat(relaxedDGoVarsValues(i, 4));
        model.var_ub(tmpIdx) = cell2mat(relaxedDGoVarsValues(i, 5));
    end
    
    solRelaxed = solveTFAmodelCplex(model);
    if isempty(solRelaxed.x)
        error(['Although we relaxed the DGo variables it was not possible'...
            'to find a feasible solution!'])
        
    elseif solRelaxed.val<minObjSolVal
        if solRelaxed.val/minObjSolVal > 0.9
            warning('Objective value after relaxation at %2.2g%% of requested value.',...
                100*solRelaxed.val/minObjSolVal)
        else
            error(['The obtained objective value after relaxation was '...
                'lower than 90%% of the requested objective value: %4.4g h^-1 (%2.2f%%)'],...
                solRelaxed.val, 100*solRelaxed.val/minObjSolVal)
        end
    end
else
    fprintf(['It is possible to obtain TFA-feasible solutions with the' ...
        'original DGo values calculated by GCM!\n'])
end

end