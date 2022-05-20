function TMinMax = runTMinMax(varargin)
% This function uses cplex to runs a min-max of the specified tFBA-model
% variables. Many of the cplex-solver parameters have been set to some
% default values, but can also be adjusted here.
% TO DO: Add other solvers as options
% runTMinMax(Model,Model.varNames(NFids),200,10^3)
% 
% 2022/05/20 Philipp Wendering: add parallel option

% if cplex is installed, and in the path
if isempty(which('cplex.m'))
    error('cplex is either not installed or not in the path')
end

% parse input
p = inputParser;

addRequired(p,'tModel')
addRequired(p,'variables')
addParameter(p,'manualScalingFactor',[])
addParameter(p,'mipTolInt',[])
addParameter(p,'emphPar',[])
addParameter(p,'feasTol',[])
addParameter(p,'scalPar',[])
addParameter(p,'TimeInSec',[])
addParameter(p,'mipDisplay',[])
addParameter(p,'runParallel',0)

p.parse(varargin{:})

tModel = p.Results.tModel;
variables = p.Results.variables;
manualScalingFactor = p.Results.manualScalingFactor;
mipTolInt = p.Results.mipTolInt;
emphPar = p.Results.emphPar;
feasTol = p.Results.feasTol;
scalPar = p.Results.scalPar;
TimeInSec = p.Results.TimeInSec;
mipDisplay = p.Results.mipDisplay;
runParallel = p.Results.runParallel;

num_vars=length(tModel.var_lb);
tModel.f = zeros(num_vars,1);
[~,varList] = ismember(variables,tModel.varNames);
varNames = tModel.varNames;

% initialization
TMinMax_LB = zeros(size(varList,1),1);
TMinMax_UB = zeros(size(varList,1),1);

if runParallel
    
    % get environment and COBRA solver
    environment = getEnvironment;
    solver = getCobraSolver('LP');
    parfor k = 1:length(varList)
        
        % restore environment and COBRA solver within parfor loop
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);
        
        i = ismember(varNames,variables{k});
        
        % Prepare cplex structure
        cplex = changeToCPLEX_WithOptions(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay);
        cplex.Model.obj = zeros(num_vars,1);
        cplex.Model.obj(i) = 1;
        cplex.Param.threads.Cur = 1;
        
        %%%% --- Minimization --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cplex.Model.sense = 'minimize';                       %
        cplexSol = cplex.solve();                             %
        solval = cplexSol.objval;                             %
        if isempty(solval)                                    %
            error('Cplex returned an empty solution!')        %
        else                                                  %
            TMinMax_LB(k) = solval;                           %
        end                                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%% --- Maximization --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cplex.Model.sense = 'maximize';                       %
        cplexSol = cplex.solve();                             %
        solval = cplexSol.objval;                             %
        if isempty(solval)                                    %
            error('Cplex returned an empty solution!')        %
        else                                                  %
            TMinMax_UB(k) = solval;                           %
        end                                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(cplex)
    end
    
    TMinMax = [TMinMax_LB TMinMax_UB];
    
else
    
    NCharsToDel = 0;
    for k = 1:length(varList)
        
        i = ismember(varNames,variables{k});
        % Prepare cplex structure
        cplex = changeToCPLEX_WithOptions(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay);
        cplex.Model.obj = zeros(num_vars,1);
        cplex.Model.obj(i) = 1;
        cplex.Param.threads.Cur = 1;
        
        %%%% --- Minimization --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(repmat('\b',1,NCharsToDel))                   %
        fprintf('Minimizing %s\n',tModel.varNames{i});        %
        cplex.Model.sense = 'minimize';                       %
        cplexSol = cplex.solve();                             %
        solval = cplexSol.objval;                             %
        if isempty(solval)                                    %
            error('Cplex returned an empty solution!')        %
        else                                                  %
            TMinMax_LB(k) = solval;                           %
        end                                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        strToDel = ['Minimizing  ',tModel.varNames{i}];
        NCharsToDel = size(strToDel,2);
        %%%% --- Maximization --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(repmat('\b',1,NCharsToDel))                   %
        fprintf('Maximizing %s\n',tModel.varNames{i});        %
        cplex.Model.sense = 'maximize';                       %
        cplexSol = cplex.solve();                             %
        solval = cplexSol.objval;                             %
        if isempty(solval)                                    %
            error('Cplex returned an empty solution!')        %
        else                                                  %
            TMinMax_UB(k) = solval;                           %
        end                                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(cplex)
    end
    
    TMinMax = [TMinMax_LB TMinMax_UB];
    

end

end
