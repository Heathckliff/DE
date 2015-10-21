function [fv_Obj, fv_Con, InformToExit, ObjInfo] = ObjfunExample(x_pop)
% Just an example!
% Input:
%       x_pop:      NoD-by-NoP population matrix.
% Output:
%       fv_Obj:     NoObj-by-NoP object function value.
%                   If none, set fv_Obj = [].
%       fv_Con:     (NoIeq+NoEq)-by-NoP constrain function value,
%                   Eq constrains must put after all the Ieq constrains.
%                   If none, set fv_Con = [];
%       ObjInfo:    1-by-1 struct contains problem information.
%       InformToExit:
%                   1-by-1 logical value informs DE when to stop. If
%                   InformToExit = TRUE, DE will stop iteration and return,
%                   else iteration will go on until InformToExit becomes
%                   TRUE or some other exit conditions are met.

% Author:
%       Yu XuanFei, Harbin Institute of Technology.
%       E-mail:     14b902037@hit.edu.cn
% Update Info:
% 2015/10/19    v0.0.0.0


if nargin == 0
    % Object function information
    ObjInfo.NoD = 2;                % Number of Dimension
    ObjInfo.NoP = 200;            	% Number of Population
    ObjInfo.NoObj = 1;              % Number of Object functions
    ObjInfo.NoIeq = 0;              % Number of Inequality constrains
    ObjInfo.NoEq = 0;               % Number of Equalilty constrains
    ObjInfo.x_b_min = -100*ones(ObjInfo.NoD,1);
    ObjInfo.x_b_max = 100*ones(ObjInfo.NoD,1);
    % dummy return values when this section call
    fv_Obj = [];
    fv_Con = [];
    InformToExit = false;
    return
end


% Evaluate Object functions (minimum)--------------------------------------
% Objf_1 = ...
% Objf_2 = ...
% ...
% Objf_NoObj = ...
Objf_1 = x_pop(1,:).^2 + x_pop(2,:).^2;


% Evaluate Inequality constrains (must <= 0 )------------------------------
% Ieqf_1 = ...
% Ieqf_2 = ...
% ...
% Ieqf_NoIeq = ...


% Evaluate Equality constranins (must = 0 )--------------------------------
% Eqf_1 = ...
% Eqf_2 = ...
% ...
% Eqf_NoEq = ...


% Asseble output matrix----------------------------------------------------
% fv_Obj = [Objf_1; Objf_2; ...; Objf_NoObj];
% fv_Con = [Ieqf_1; Ieqf_2; ...; Ieqf_NoIeq; Eqf_1; Eqf_2; ...; Eqf_NoEq];
fv_Obj = Objf_1;
fv_Con = [];


% Evaluate when to stop DE-------------------------------------------------
% If (conditions are met)
%       InformToExit = true;
% else
%       InformToExit = false;
% end
InformToExit = false;


end
