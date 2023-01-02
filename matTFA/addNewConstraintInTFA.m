function model = addNewConstraintInTFA(model, CName, CType, CLHS, CRHS, appendFlag, rowIdx)
% This function adds one new constraint to the model, and specifies how the
% associated variables are involved in this constraint

if ~exist('appendFlag', 'var') || isempty(appendFlag)
    appendFlag = true;
elseif ~exist('rowIdx', 'var') || isempty(rowIdx)
    appendFlag = false;
end

[num_constr,~] = size(model.A);
if appendFlag
    model.constraintNames{num_constr+1,1}   = CName; % strcat('DFSEU_',model.mets{i});
    model.constraintType{num_constr+1,1}    = CType; % type of constraint: '<', or '>', or '='
    model.rhs(num_constr+1)                 = CRHS; % value of the right hand side
    model.A(num_constr+1,CLHS.varIDs)       = CLHS.varCoeffs; % coefficients of the involved variables
else
    model.constraintNames{rowIdx}   = CName; % strcat('DFSEU_',model.mets{i});
    model.constraintType{rowIdx}    = CType; % type of constraint: '<', or '>', or '='
    model.rhs(rowIdx)               = CRHS; % value of the right hand side
    newRows                         = sparse(numel(rowIdx), size(model.A,2));
    newRows(CLHS.varIDs)            = CLHS.varCoeffs;
    model.A(rowIdx,:)               = newRows; % coefficients of the involved variables
end
end