function [winner] = ParetoSelect(fv_Obj_baby,fv_Obj,fv_Con_baby,fv_Con,ObjInfo)
% Pareto dominance based selection with direct constraint handling.
% Input:
%       fv_Obj_baby:	NoObj-by-NoP baby's object matrix.
%       fv_Obj:         NoObj-by-NoP parent's object matrix.
%       fv_Con_baby:    (NoIeq+NoEq)-by-NoP baby's constrain matrix.
%       fv_Con:         (NoIeq+NoEq)-by-NoP parent's constrain matrix.
%       ObjInfo:        1-by-1 struct contains problem information.
% Output:
%       winner:         1-by-NoP logical matrix indicates who is allowed
%                       to enter the next generation. TRUE means baby win,
%                       and FALSE means parent win.

% Author:
%       Yu XuanFei, Harbin Institute of Technology.
%       E-mail:     14b902037@hit.edu.cn
% Update Info:
% 2015/10/19    v0.0.0.0


NoD = ObjInfo.NoD;              % Number of Dimension
NoObj = ObjInfo.NoObj;          % Number of Object functions
NoIeq = ObjInfo.NoIeq;          % Number of Inequality constrains
NoEq = ObjInfo.NoEq;            % Number of Equalilty constrains

% Default values set to parent win.
if isempty(fv_Obj)
    rc_best = ones(1,length(fv_Con(1,:))) < 0;
elseif isempty(fv_Con)
    rc_best = ones(1,length(fv_Obj(1,:))) < 0;
end
% Conditions baby will win.
if (NoIeq + NoEq) >= 1% Constrained Problem
  % Condition 1
    rc = (sum(fv_Con_baby <= 0,1) == (NoIeq + NoEq)) &...
         (sum(fv_Con <= 0,1) == (NoIeq + NoEq));
       if NoObj >= 1% Objective Optimization Problem
           rc = rc & (sum(fv_Obj_baby <= fv_Obj,1) == NoObj);
           rc_best(rc) = true;
       else% Constraint Satisfaction Problem
           rc_best(rc) = true;
       end
  % Condition 2
	rc =  (sum(fv_Con_baby <= 0,1) == (NoIeq + NoEq)) &...
          (sum(fv_Con <= 0,1) ~= (NoIeq + NoEq));
	rc_best(rc) = true;
  % Condition 3
    rc =  (sum(fv_Con_baby <= 0,1) ~= (NoIeq + NoEq)) &...
          (sum(fv_Con <= 0,1) ~= (NoIeq + NoEq)) &...
          (sum(fv_Con_baby <= fv_Con,1) == (NoIeq + NoEq));
	rc_best(rc) = true;
else% Unconstrained Problem
    rc = sum(fv_Obj_baby <= fv_Obj,1) == NoObj;
    rc_best(rc) = true;
end
%     winner = ones(NoD,1)*rc_best == 1;
    winner = rc_best;
end