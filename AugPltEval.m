%% Evaluation of a generalized plant by plugging in theta value
function Gplt = AugPltEval(Gplt_sym,Par)
% Gplt_sym is augmented plant with symbol variables;
% Gplt is augmented plant evaluated at selected theta
if size(theta,1) == 1
    Gplt.A = double(subs(Gplt_sym.A,theta));   
    Gplt.B1 = double(subs(Gplt_sym.B1,theta));   
    Gplt.B2 = double(subs(Gplt_sym.B2,theta));   
    Gplt.C1 = double(subs(Gplt_sym.C1,theta));   
    Gplt.D11 = double(subs(Gplt_sym.D11,theta));   
    Gplt.D12 = double(subs(Gplt_sym.D12,theta));   
    Gplt.C2 = double(subs(Gplt_sym.C2,theta));   
    Gplt.D21 = double(subs(Gplt_sym.D21,theta));   
    Gplt.D22 = double(subs(Gplt_sym.D22,theta));   
elseif size(theta,1) == 2
    theta1 = theta(1); theta2 = theta(2);
    Gplt.A = double(subs(Gplt_sym.A));   
    Gplt.B1 = double(subs(Gplt_sym.B1));   
    Gplt.B2 = double(subs(Gplt_sym.B2));   
    Gplt.C1 = double(subs(Gplt_sym.C1));   
    Gplt.D11 = double(subs(Gplt_sym.D11));   
    Gplt.D12 = double(subs(Gplt_sym.D12));   
    Gplt.C2 = double(subs(Gplt_sym.C2));   
    Gplt.D21 = double(subs(Gplt_sym.D21));   
    Gplt.D22 = double(subs(Gplt_sym.D22));   
elseif size(theta,1) == 3
    theta1 = theta(1); theta2 = theta(2); theta3 = theta(3);
    Gplt.A = double(subs(Gplt_sym.A));   
    Gplt.B1 = double(subs(Gplt_sym.B1));   
    Gplt.B2 = double(subs(Gplt_sym.B2));   
    Gplt.C1 = double(subs(Gplt_sym.C1));   
    Gplt.D11 = double(subs(Gplt_sym.D11));   
    Gplt.D12 = double(subs(Gplt_sym.D12));   
    Gplt.C2 = double(subs(Gplt_sym.C2));   
    Gplt.D21 = double(subs(Gplt_sym.D21));   
    Gplt.D22 = double(subs(Gplt_sym.D22));   
end
