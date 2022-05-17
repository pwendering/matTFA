function minmax = runMinMax(model, rxnNames, verbose, scalPar, feasTol, emphPar)

if ~exist('scalPar','var') || isempty(scalPar)
    scalPar = [];
end
if ~exist('feasTol','var') || isempty(feasTol)
    feasTol = [];
end
if ~exist('emphPar','var') || isempty(emphPar)
    emphPar = [];
end
if ~exist('rxnNames','var') || isempty(rxnNames)
    rxnNames = model.rxns;
end
if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end

minmax = zeros(length(rxnNames),2);
rxn_id = find_cell(rxnNames, model.rxns);

environment = getEnvironment;
solver = getCobraSolver('LP');
parfor i = 1:length(rxn_id)
    restoreEnvironment(environment);
    changeCobraSolver(solver, 'LP', 0, -1);
    tmp_model = model;    
    sol1 = solveFBAmodelCplex(tmp_model, scalPar, feasTol, emphPar, 'min');
    sol2 = solveFBAmodelCplex(tmp_model, scalPar, feasTol, emphPar, 'max');
    minmax(i,:) = [sol1.f sol2.f]; 
end
