function model = addNewConstraintInTFA(model, CName, CType, CLHS, CRHS, appendFlag)
% This function adds one new constraint to the model, and specifies how the
% associated variables are involved in this constraint

if ~exist('appendFlag', 'var') || isempty(appendFlag)
    appendFlag = true;
end

[num_constr,~] = size(model.A);
if appendFlag
    model.constraintNames{num_constr+1,1}   = CName; % strcat('DFSEU_',model.mets{i});
    model.constraintType{num_constr+1,1}    = CType; % type of constraint: '<', or '>', or '='
    model.rhs(num_constr+1)                 = CRHS; % value of the right hand side
    model.A(num_constr+1,CLHS.varIDs)       = CLHS.varCoeffs; % coefficients of the involved variables
else
    % find last non-zero row and add variable at position of first
    % all-zero row
    rowIdx = find(all(model.A,2),1,'last') + 1;
    model.constraintNames{rowIdx}   = CName; % strcat('DFSEU_',model.mets{i});
    model.constraintType{rowIdx}    = CType; % type of constraint: '<', or '>', or '='
    model.rhs(rowIdx)               = CRHS; % value of the right hand side
    model.A(rowIdx,CLHS.varIDs)     = CLHS.varCoeffs; % coefficients of the involved variables
end
end