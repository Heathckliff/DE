function [fv_Obj, fv_Con, InformToExit, ObjInfo] = ObjfunModel(x_pop)
% Problem definnning, put your problem mathmatical description here!
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


if nargin == 0 % Output ObjInfo only.
    % Object function information
    ObjInfo.NoD = 2;                % Number of Dimension
    ObjInfo.NoP = 10*ObjInfo.NoD;  	% Number of Population
    ObjInfo.NoObj = 0;              % Number of Object functions
    ObjInfo.NoIeq = 2;              % Number of Inequality constrains
    ObjInfo.NoEq = 0;               % Number of Equalilty constrains
    ObjInfo.x_b_min = [-10;-10];  	% Lower boundary constrain,NoD-by-1
    ObjInfo.x_b_max = [10;10];    	% Upper boundary constrain,NoD-by-1

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


% Evaluate Inequality constrains (must <= 0 )------------------------------
% Ieqf_1 = ...
% Ieqf_2 = ...
% ...
% Ieqf_NoIeq = ...
Ieq_f1 = 7 - sqrt(x_pop(1,:).^2 + x_pop(2,:).^2);
Ieq_f2 = -8 + sqrt(x_pop(1,:).^2 + (2.*x_pop(2,:)).^2);


% Evaluate Equality constranins (must = 0 )--------------------------------
% Eqf_1 = ...
% Eqf_2 = ...
% ...
% Eqf_NoEq = ...


% Asseble output matrix----------------------------------------------------
% fv_Obj = [Objf_1; Objf_2; ...; Objf_NoObj];
% fv_Con = [Ieqf_1; Ieqf_2; ...; Ieqf_NoIeq; Eqf_1; Eqf_2; ...; Eqf_NoEq];
fv_Obj = [];
fv_Con = [Ieq_f1; Ieq_f2];


% Evaluate when to stop DE-------------------------------------------------
% If (conditions are met)
%       InformToExit = true;
% else
%       InformToExit = false;
% end
InformToExit = false;


end