function model = addNewVariableInTFA(model, VName, VType, VRange, appendFlag)
% This function adds one new variable to the model. It creates a zero
% vector in the model.A matrix. The constraints will specify how this
% variable affects the system.

if ~exist('appendFlag', 'var') || isempty(appendFlag)
    appendFlag = true;
end
    
[num_constr,num_vars] = size(model.A);
if appendFlag
    model.varNames{num_vars+1} = VName; % append new variable name
    model.var_lb(num_vars+1)   = VRange(1); % lower bound
    model.var_ub(num_vars+1)   = VRange(2); % upper bound
    model.vartypes{num_vars+1} = VType; % 'C' or 'B'
    model.A(num_constr,num_vars+1) = 0; % increase the columns of model.A by one
else
    % find last non-zero column and add variable at position of first
    % all-zero column
    colIdx = find(all(model.A),1,'last') + 1;
    model.varNames{colIdx} = VName; % fill in variable name
    model.var_lb(colIdx)   = VRange(1); % fill in lower bound
    model.var_ub(colIdx)   = VRange(2); % fill in upper bound
    model.vartypes{colIdx} = VType; % fill in variable type
end
end