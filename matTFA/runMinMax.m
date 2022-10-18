function minmax = runMinMax(varargin)
% 2022/05/20 Philipp Wendering: add parallel option

% parse input
p = inputParser;

addRequired(p,'model')
addRequired(p,'rxnNames')
addParameter(p,'verbose',false)
addParameter(p,'scalPar',[])
addParameter(p,'feasTol',[])
addParameter(p,'emphPar',[])
addParameter(p,'runParallel',0)

p.parse(varargin{:})

model = p.Results.model;
rxnNames = p.Results.rxnNames;
verbose = p.Results.verbose;
scalPar = p.Results.scalPar;
feasTol = p.Results.feasTol;
emphPar = p.Results.emphPar;
runParallel = p.Results.runParallel;

minmax = zeros(length(rxnNames),2);
rxn_id = find_cell(rxnNames, model.rxns);

if runParallel
    
    % get environment and COBRA solver
    environment = getEnvironment;
    solver = getCobraSolver('LP');
    
    parfor i = 1:length(rxn_id)
        
        % restore environment and COBRA solver within parfor loop
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);
        tmp_model = model;
        sol1 = solveFBAmodelCplex(tmp_model, scalPar, feasTol, emphPar, 'min');
        sol2 = solveFBAmodelCplex(tmp_model, scalPar, feasTol, emphPar, 'max');
        minmax(i,:) = [sol1.f sol2.f];
        
    end
else
    for i = 1:length(rxn_id)
        
        if ~isfield(model,'CS_varNames')
            if verbose
                fprintf('minmax for %s\t',model.rxns{rxn_id(i)})
            end
            model = changeObjective(model,model.rxns{rxn_id(i)});
        elseif isfield(model,'CS_varNames')
            if verbose
                fprintf('minmax for %s\t',model.CS_varNames{rxn_id(i)})
            end
            model.c = zeros(size(model.S,2),1);
            model.c(rxn_id(i)) = 1;
        end
        
        sol = solveFBAmodelCplex(model, scalPar, feasTol, emphPar, 'min');
        minmax(i,1) = sol.f;
        sol = solveFBAmodelCplex(model, scalPar, feasTol, emphPar, 'max');
        minmax(i,2) = sol.f;
        
        if verbose
            fprintf('min: %d\t max: %d\n', minmax(i,1), minmax(i,2));
        end
        
    end
end
