function [rc_best] = ParetoBest(fv_Obj,fv_Con,NBest,ObjInfo)
% Find the locations of the n best-so-far solutions based on
% Parero-Dominance.
% Input:
%       fv_Obj:     NoObj-by-NoP object function value.
%                   If no Object function, set fv_Obj = [].
%       fv_Con:     (NoIeq+NoEq)-by-NoP constrain function value,
%                   Eq constrains must put after all the Ieq constrains.
%                   If no constrains, set fv_Con = [];
%       NBest:      Number of the best-so-far solutions required.
%       ObjInfo:    1-by-1 struct contains problem information.
% Output:
%       rc_best:    Location of the best-so-far solution.

% Author:
%       Yu XuanFei, Harbin Institute of Technology.
%       E-mail:     14b902037@hit.edu.cn
% Update Info:
% 2015/10/19    v0.0.0.0


NoP = ObjInfo.NoP;              % Number of Population
NoObj = ObjInfo.NoObj;          % Number of Object functions
NoIeq = ObjInfo.NoIeq;          % Number of Inequality constrains
NoEq = ObjInfo.NoEq;            % Number of Equalilty constrains

if ((NoIeq + NoEq) == 0) && (NoObj == 1) % Unconstrained Single Object Problem
    [~, rc_best] = sort(fv_Obj);
    rc_best = rc_best(1:NBest);
    return
end

rc_best = ones(1,NBest);
temp = 1:NoP;

for jj = 1:NBest
	rc_best(jj) = temp(1);
    for ii = temp
        if (NoIeq + NoEq) == 0  % Unconstrained Problem
            winner = ParetoSelect(fv_Obj(:,ii),fv_Obj(:,rc_best(jj)),[],[],ObjInfo);
        elseif NoObj == 0       % Constrain Satisfaction Problem
            winner = ParetoSelect([],[],fv_Con(:,ii),fv_Con(:,rc_best(jj)),ObjInfo);
        else                    % Constrained Optimization
            winner = ParetoSelect(fv_Obj(:,ii),fv_Obj(:,rc_best(jj)),fv_Con(:,ii),fv_Con(:,rc_best(jj)),ObjInfo);
        end

        if winner
            rc_best(jj) = ii;
        end
    end
    temp = temp(temp ~= rc_best(jj));
end