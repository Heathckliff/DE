function [x_parent,fv_Obj,fv_Con,InformToExit,ExitFlag] = DE(DeInfo,datapath,PlotIt)
% Defferential Evolution Optimization Toolbox.
% Problem can be solved:
%   1. Unconstrained sigle/multiple object optimization.
%   2. Constrained single/multiple object optimization.
%   3. Constrain satisfaction problem.
% Input:
%       DeInfo:     1-by-1 struct contains DE options, include the
%                   following fields:
%       FunObjName: Function name which contains mathmatical description of
%                   the problem.
%       MaxGen:     Maximum generation allowed.
%       F:          Differential Mutation scale factor.
%       Cr:         Crossover probability.
%       EqErr:      |fv_Con(Eq)| <= EqErr.
%       BaseVectorStrategy:
%                   'rand' | 'best' | 'target-to-best' | 'n-best'
%       ScaleFactorStrategy:
%                   'constant' | 'jitter' | 'dither-per-vector' |
%                   'dither-per-generation'
%       ScaleFactorRandomizeMagnitude:
%                   must < 2*F
%       EitherOrFactor:
%                   1 means pure differential mutation,
%                   0 means pure arithmatic recombination,
%                   value between (0,1) means the algorithm of 'either-or'.
%       TargetToBestFactor:
%                   If 'ScaleFactorStrategy' set to 'target-to-best', this
%                   option must be set a value between (0,1).
%       NBest:      If 'ScaleFactorStrategy' set to 'n-best', NBest must be
%                   set to an integer value <= NoP.
%       datapath:   If 'datapath' is [], DE will start with a randomly
%                   initialized population; otherwise, DE will proceed with
%                   with the provided population pointed by datapath.
%       PlotIt:     If PlotIt != 0, plot DE process.
% Output:
%       x_parent:   NoD-by-NoP population matrix of the last generation.
%       fv_Obj:     NoObj-by-NoP object function values of the last 
%                   generation. If no Object function, return fv_Obj = [].
%       fv_Con:     (NoIeq+NoEq)-by-NoP constrain function values of the
%                   last generation. If no constrains, return fv_Con = [].
%       ObjInfo:    1-by-1 struct contains problem information.
%       ExitFlag:   Indicates which exit condition is matched.

% Author:
%       Yu XuanFei, Harbin Institute of Technology.
%       E-mail:     14b902037@hit.edu.cn
% Update Info:
% 2015/10/19    v0.0.0.0


fprintf('# DE...\n')
%% Get DE Options
fprintf('   # Get DE Options...\n')
FunObjName = DeInfo.FunObjName;
MaxGen =DeInfo.MaxGen;
F = DeInfo.F;
Cr = DeInfo.Cr;
BaseVectorStrategy =DeInfo.BaseVectorStrategy;
ScaleFactorStrategy =DeInfo.ScaleFactorStrategy;
ScaleFactorRandomizeMagnitude =DeInfo.ScaleFactorRandomizeMagnitude;
EitherOrFactor =DeInfo.EitherOrFactor;
TargetToBestFactor =DeInfo.TargetToBestFactor;
NBest =DeInfo.NBest;
EqErr =DeInfo.EqErr;                  % Equality constrain precision
%% Get Object Problem Information:
fprintf('   # Get Object Information...\n')
[~, ~, ~, ObjInfo] = feval(FunObjName);
NoD = ObjInfo.NoD;              % Number of Dimension
NoP = ObjInfo.NoP;              % Number of Population
NoObj = ObjInfo.NoObj;          % Number of Object functions
NoIeq = ObjInfo.NoIeq;          % Number of Inequality constrains
NoEq = ObjInfo.NoEq;            % Number of Equalilty constrains
x_b_min = ObjInfo.x_b_min;
x_b_max = ObjInfo.x_b_max;
ScaleObjFunVal = ObjInfo.ScaleObjFunVal;
%% Initialization population
if isempty(datapath)
    fprintf('   # Initialization...\n')
    rng('shuffle');
    x_parent = zeros(NoD,NoP);
    for kk = 1:NoD
        x_parent(kk,:) = x_b_min(kk) + (x_b_max(kk) - x_b_min(kk)).*rand(1,NoP);
    end

    fprintf('   # Evaluate Objfun(initialize)...\n')
    [fv_Obj, fv_Con] = feval(FunObjName, x_parent, ObjInfo);
    if NoEq >= 1
        fv_Con(NoIeq+1:NoIeq+NoEq,:) = abs(fv_Con(NoIeq+1:NoIeq+NoEq,:)) - EqErr;
    end
    rc1 = fv_Con <= 0;
    fv_Con(rc1) = 0;
else
    load(datapath,'-mat')
end
%% The Main Cycle
fprintf('   # DE Main Cycle...\n')
jj = 1;
while true
fprintf('   %d\n',jj)
%% Method of Sacle Factor
if strcmp(ScaleFactorStrategy,'jitter')
    F_j = F + ScaleFactorRandomizeMagnitude*(rand(NoD,NoP) - 0.5);
elseif strcmp(ScaleFactorStrategy,'dither-per-vector')
    F_j = ones(NoD,1)*(F + ScaleFactorRandomizeMagnitude*(rand(1,NoP) - 0.5));
elseif strcmp(ScaleFactorStrategy,'dither-per-generation')
    F_j = F + ScaleFactorRandomizeMagnitude*(rand - 0.5);
else
    F_j = F;
end
%% Method of Mutation
% Choose shuffle location r1, r2, r3
r = randi(NoP,1,3);
while r(1)==r(2) || r(2)==r(3) || r(3)==r(1) || r(1)==NoP || r(2)==NoP || r(3)==NoP
    r = randi(NoP,1,3);
end
% Choose base vector
if strcmp(BaseVectorStrategy,'rand')
    x_base = x_parent(:,[r(1)+1:end,1:r(1)]);
elseif strcmp(BaseVectorStrategy,'best') || strcmp(BaseVectorStrategy,'target-to-best')
    fprintf('   # Evaluate ParetoBest...\n')
    rc_best = ParetoBest(fv_Obj,fv_Con,1,ObjInfo);
    x_base = x_parent(:,rc_best)*ones(1,NoP);
    if strcmp(BaseVectorStrategy,'target-to-best')
        x_base = x_parent + TargetToBestFactor.*(x_base - x_parent);
    end
elseif strcmp(BaseVectorStrategy,'n-best')
    fprintf('   # Evaluate ParetoBest...\n')
    rc_best = ParetoBest(fv_Obj,fv_Con,NBest,ObjInfo);
    x_base = x_parent;
    for kk = 1:floor(NoP/NBest)
        x_base(:,1+(kk-1)*NBest:kk*NBest) = x_parent(:,rc_best(randperm(NBest)));
    end
    x_base(:,kk*NBest+1:end) = x_parent(:,rc_best(randperm(NBest,NoP-kk*NBest)));
end
% Way to mutate
if rand < EitherOrFactor    % Defferential Mutation
    x_mut = x_base + ...
    	F_j.*(x_parent(:,[r(2)+1:end,1:r(2)]) - x_parent(:,[r(3)+1:end,1:r(3)]));
else                        % Arithmetic Recombination
    x_mut = x_base + ...
            0.5.*(F_j + 1).*(x_parent(:,[r(2)+1:end,1:r(2)]) + ...
            x_parent(:,[r(3)+1:end,1:r(3)]) - 2.*x_base);
end
% Evaluate Boundary constrains(Bounce-Back)
for kk = 1:NoD
    rc1 = x_mut(kk,:) < x_b_min(kk);
    x_mut(kk,rc1) = 0.5.*(x_base(kk,rc1) + x_b_min(kk));
    
    rc1 = x_mut(kk,:) > x_b_max(kk);
    x_mut(kk,rc1) = 0.5.*(x_base(kk,rc1) + x_b_max(kk));
end
%% Crossover
x_baby = x_parent;
% if rand(0,1) <= Cr
rc1 = rand(NoD,NoP) <= Cr;
x_baby(rc1) = x_mut(rc1);
% or if j = randi(1,Np)
idx1 = sub2ind(size(x_baby),randi(NoD,1,NoP),1:NoP);
x_baby(idx1) = x_mut(idx1);

fprintf('   # Evaluate Objfun(baby)...\n')
[fv_Obj_baby, fv_Con_baby] = feval(FunObjName, x_baby, ObjInfo);
if NoEq >= 1
    fv_Con_baby(NoIeq+1:NoIeq+NoEq,:) = abs(fv_Con_baby(NoIeq+1:NoIeq+NoEq,:)) - EqErr;
end
rc1 = fv_Con_baby <= 0;
fv_Con_baby(rc1) = 0;
%% Selection for next generation
fprintf('   # Evaluate ParetoSelect...\n')
rc1 = ParetoSelect(fv_Obj_baby,fv_Obj,fv_Con_baby,fv_Con,ObjInfo);
x_parent(:,rc1) = x_baby(:,rc1);

sum_rc1 = sum(rc1);

fprintf('   # Evaluate Objfun(nextGen)...\n')
if NoObj >= 1 && (NoIeq+NoEq) >= 1 && (sum_rc1 ~= 0)
    [fv_Obj(:,rc1), fv_Con(:,rc1), InformToExit] = feval(FunObjName, x_parent(:,rc1), ObjInfo);
elseif NoObj >= 1 && (NoIeq+NoEq) == 0 && (sum_rc1 ~= 0)
    [fv_Obj(:,rc1), fv_Con, InformToExit] = feval(FunObjName, x_parent(:,rc1), ObjInfo);
elseif NoObj == 0 && (NoIeq+NoEq) >= 1 && (sum_rc1 ~= 0)
    [fv_Obj, fv_Con(:,rc1), InformToExit] = feval(FunObjName, x_parent(:,rc1), ObjInfo);
end
if NoEq >= 1
    fv_Con(NoIeq+1:NoIeq+NoEq,:) = abs(fv_Con(NoIeq+1:NoIeq+NoEq,:)) - EqErr;
end
rc1 = fv_Con <= 0;
fv_Con(rc1) = 0;
%% Output section
if NoObj >= 1
%     norm_fv_Obj = max(abs(fv_Obj),[],2);
    norm_fv_Obj = mean(fv_Obj,2)./ScaleObjFunVal;
    fprintf('           normInf(Obj_f%02d) = %e\n', [1:NoObj; norm_fv_Obj'])
    %
    if PlotIt
        figure(1); hold on; grid on; box on;
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot(jj,norm_fv_Obj,'o');
        drawnow
    end
else
    Violations = NoP - sum(sum(rc1,1) == (NoIeq + NoEq));
    fprintf('           Violations: %d\n', Violations)
    %
    if PlotIt
        if NoD >= 2
            figure(1); hold on; grid on; box on;
            plot(x_parent(1,:),x_parent(2,:),'.');
            drawnow
        else
            figure(1); hold on; grid on; box on;
            plot(x_parent,'.');
            drawnow
        end
    end
end
%% Update
jj = jj + 1;
%% Termination Criteria
if InformToExit% == true    % Stop by user defined criterion
    ExitFlag = 0;
    return
elseif (jj > MaxGen)        % Maximum generation reached
    ExitFlag = 1;
    return
elseif NoObj == 0 && (sum(sum(fv_Con <= 0,1))/NoP == (NoIeq + NoEq))
                            % For constrain satisfaction problem, all
                            % constrains are met
    ExitFlag = 2;
    return
end
end