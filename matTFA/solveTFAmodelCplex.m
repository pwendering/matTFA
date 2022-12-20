function sol = solveTFAmodelCplex(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp)
% Solve a model using specific solver settings. More details in
% changeToCPLEX_WithOptions
%
% INPUTS
%   tModel: a TFA-ready model (has a .A matrix)
%   TimeInSec: timelimit for the solver. 
%   manualScalingFactor: manual scaling factor for the solver
%   mipTolInt: Integer tolerance of the solver
%   emphPar: Solver emphasis (trade-offs between speed, feasibility, 
%       optimality, and moving bounds in MIP - see https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/MIPEmphasis.html)
%   feasTol: solver tolerance for feasibility (error on constraints)
%   mipDisplay: verbosity of the MIP info display
%   CPXPARAMdisp: Turn on/off the cplex problem setup display
% 
%% Changelog
% 2017/04/26 - Modified by Pierre on Georgios Fengos's base, to incorporate
% Vikash's Gurobi hooks in a more global fashion, in a similar way COBRA
% solvers are handled.
% Calls the global parameter TFA_MILP_SOLVER, and check if it set to
% something. If not, default to CPLEX using Fengos's code.
% In the long run, we should rename thins function to remove CPLEX from its
% name.

global TFA_MILP_SOLVER

if ~exist('TFA_MILP_SOLVER','var') || isempty(TFA_MILP_SOLVER)
    TFA_MILP_SOLVER = 'cplex_direct';
end
solver = TFA_MILP_SOLVER;
% solver = 'gurobi_direct';
if ~exist('manualScalingFactor','var') || isempty(manualScalingFactor)
    manualScalingFactor = [];
end
if ~exist('mipTolInt','var') || isempty(mipTolInt)
    mipTolInt = [];
end
if ~exist('emphPar','var') || isempty(emphPar)
    emphPar = [];
end
if ~exist('feasTol','var') || isempty(feasTol)
    feasTol = [];
end
if ~exist('scalPar','var') || isempty(scalPar)
    scalPar = [];
end
if ~exist('TimeInSec','var') || isempty(TimeInSec)
    TimeInSec = [];
end
if ~exist('mipDisplay','var') || isempty(mipDisplay)
    mipDisplay = [];
end
if ~exist('CPXPARAMdisp','var') || isempty(CPXPARAMdisp)
    CPXPARAMdisp = [];
end

switch solver
    %% Case CPLEX
    case 'cplex_direct'
        if isempty(which('cplex.p'))
            error('You need to add CPLEX to the Matlab-path!!')
        end
        sol = x_solveCplex(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp);
    %% Case GUROBI
    case 'gurobi_direct'
        if isempty(which('gurobi'))
            error('You need to add Gurobi to the Matlab-path!!')
        end
        sol = x_solveGurobi(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay);
end
end

%% Private function for CPLEX solve
function sol = x_solveCplex(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp)

% this function solves a TFBA problem using CPLEX

% if cplex is installed, and in the path
if isempty(which('cplex.m'))
    error('cplex is either not installed or not in the path')
end

% Convert problem to cplex
cplex = changeToCPLEX_WithOptions(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp);

% Optimize the problem
try
    CplexSol = cplex.solve();
    if isfield(cplex.Solution,'x')
        x = cplex.Solution.x;
        if ~isempty(x)
            sol.x = cplex.Solution.x;
            sol.val = cplex.Solution.objval;
            sol.cplexSolStatus = cplex.Solution.status;
        else
            sol.x = [];
            sol.val = [];
            disp('Empty solution');
            warning('Cplex returned an empty solution!')
            sol.cplexSolStatus = 'Empty solution';
        end
    else
        sol.x = [];
        sol.val = [];
        disp('No field cplex.Solution.x');
        warning('The solver does not return a solution!')
        sol.cplexSolStatus = 'No field cplex.Solution.x';
    end
catch
    sol.x = NaN;
    sol.val = NaN;
    sol.cplexSolStatus = 'Solver crashed';
end

delete(cplex)
end

%% Private function for GUROBI solve
function sol = x_solveGurobi(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay)
    
    if ~exist('manualScalingFactor','var') || isempty(manualScalingFactor)
        % Do not scale the problem manually
    else
        tModel.A = manualScalingFactor*tModel.A;
        tModel.rhs = manualScalingFactor*tModel.rhs;
    end

    % set gurobi parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TimeLimit
    % Limits the total time expended (in seconds). Optimization returns 
    % with a TIME_LIMIT status if the limit is exceeded.
    % Type:             double
    % Default value:    Infinity
    % Minimum value:    0
    % Maximum value:    Infinity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('TimeInSec','var') || isempty(TimeInSec)
        % LCSB default
        TimeInSec = 1e75;
    else
        if TimeInSec < 0 || TimeInSec > 1e+75
            error('Parameter value out of range!')
        end
    end
    params.TimeLimit = TimeInSec;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FeasibilityTol: Primal feasibility tolerance
    % All constraints must be satisfied to a tolerance of FeasibilityTol.
    % Tightening this tolerance can produce smaller constraint violations, 
    % but for numerically challenging models it can sometimes lead to much 
    % larger iteration counts.
    % Type:             double
    % Default value:    1e-6
    % Minimum value:    1e-9
    % Maximum value:    1e-2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('feasTol','var') || isempty(feasTol)
        % LCSB default
        feasTol = 1e-9;
    else
        if feasTol < 1e-9 || feasTol > 1e-2
            error('Parameter value out of range!')
        end
    end
    params.FeasibilityTol = feasTol;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OutputFlag: Controls Gurobi output
    % Enables or disables solver output.
    % Type:             int
    % Default value:    1
    % Minimum value:    0
    % Maximum value:    1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('mipDisplay','var') || isempty(mipDisplay)
        mipDisplay = 0;
    else
        if ~ismember(mipDisplay,[0 1])
            error('Parameter value out of range!')
        end
    end
    params.OutputFlag = mipDisplay;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IntFeasTol: Integer feasibility tolerance
    % An integrality restriction on a variable is considered satisfied when
    % the variable’s value is less than IntFeasTol from the nearest integer
    % value. Tightening this tolerance can produce smaller integrality 
    % violations, but very tight tolerances may significantly increase
    % runtime. Loosening this tolerance rarely reduces runtime.
    % Type:             double
    % Default value:    1e-5
    % Minimum value:    1e-9
    % Maximum value:    1e-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('mipTolInt','var') || isempty(mipTolInt)
        mipTolInt = 1e-9;
    else
        if mipTolInt > 0.1 || mipTolInt < 0
            error('Parameter value out of range!')
        end
    end
    params.IntFeasTol = mipTolInt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NumericFocus: Numerical focus
    % The NumericFocus parameter controls the degree to which the code 
    % attempts to detect and manage numerical issues.
    % The default setting (0) makes an automatic choice, with a slight 
    % preference for speed.
    % Settings 1-3 increasingly shift the focus towards being more careful 
    % in numerical computations. With higher values, the code will spend 
    % more time checking the numerical accuracy of intermediate results,
    % and it will employ more expensive techniques in order to avoid 
    % potential numerical issues.
    % Type:             int
    % Default value:    0
    % Minimum value:    0
    % Maximum value:    3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('emphPar','var') || isempty(emphPar)
        emphPar = 3;
    else
        if ~ismember(emphPar,0:3)
            error('Parameter value out of range!')
        end
    end
    params.NumericFocus = emphPar;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ScaleFlag: Model scaling
    % Controls model scaling. By default, the rows and columns of the model
    % are scaled in order to improve the numerical properties of the 
    % constraint matrix. The scaling is removed before the final solution 
    % is returned. Scaling typically reduces solution times, but it may 
    % lead to larger constraint violations in the original, unscaled model.
    % Turning off scaling (ScaleFlag=0) can sometimes produce smaller 
    % constraint violations. Choosing a different scaling option can 
    % sometimes improve performance for particularly numerically difficult 
    % models.
    % Type:             int
    % Default value:    -1
    % Minimum value:    -1
    % Maximum value:    3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('scalPar','var') || isempty(scalPar)
        scalPar = -1;
    else
        if ~ismember(scalPar, -1:3)
            error('Parameter value out of range!')
        end
    end
    params.ScaleFlag = scalPar;
    
    num_constr = length(tModel.constraintType);
    num_vars = length(tModel.vartypes);

    contypes = '';
    vtypes = '';

    % convert contypes and vtypes into the right format
    for i=1:num_constr
        contypes = strcat(contypes,tModel.constraintType{i,1});
    end
    
    for i=1:num_vars
        vtypes = strcat(vtypes,tModel.vartypes{i,1});
    end
    
    gmodel.A=tModel.A;
    gmodel.obj=tModel.f;
    gmodel.lb=tModel.var_lb;
    gmodel.ub=tModel.var_ub;
    gmodel.rhs=tModel.rhs;
    gmodel.sense=contypes;
    gmodel.vtype=vtypes;
  
    gmodel.varnames=tModel.varNames;
    
    if tModel.objtype==-1
      gmodel.modelsense='max';
    elseif tModel.objtype==1
        gmodel.modelsense='min';
    else
        error(['No objective type specified ' ...
            '(model.objtype should be in {-1,1})']);
    end
    
    try
        result=gurobi(gmodel, params);
        if isfield(result,'x')
            x = result.x;
            x(find(abs(x) < 1E-9))=0;
            
        else
           
            warning('The solver does not return a solution!')
            result.x=[];
            result.objval=[];
        end
    catch
        result.status='0';
        x=NaN;
        result.x=NaN;
        result.objval=NaN;
    end
    
    sol.x=result.x;
    sol.val=result.objval;
    sol.status=result.status;
    % TODO: Add exitflag translation

end