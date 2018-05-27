%*************************************************************************
%  
%  Switching LPV controller design under uncertain scheduling parameters 
%             ----------------------------------------
%
% By Pan Zhao, Control Engineering Lab, UBC, BC, Canada
% Under supervision of Prof. Ryozo Nagamune.
% Creation: May 12, 2015.
% Revision: Oct 3, 2015.
% 
% Function: Main function.
%               
% Note: 
%   * Can be used for both switching and non-switching LPV controller design
%   * Can be used for both one and two scheduling parameters,
%   * Questions can be sent to panbrzhao@gmail.com or panzhao@mech.ubc.ca.
clear; clc;
tic;
%% Parameters 
XY_PD = 2           % 0 for parameter-independent X and Y, 1 for parameter-dependent X and constant Y; 2 for parameter-dependent Y and constant X
SwLogic = 2         % 0 for non-switching, 1 for hesteresis switching, 2 for switching with average dwell time
ASWCond_Ceq = 1;    % 1 for using constant terms (X or Y) for equality (default), 0 for PD terms for equality, only valid for average dwell time switching. 
IsGridding = 0;     % 0 for non-gridding method (affine LMIs), 1 for gridding method (Apkarian) method
GSUncertType = 1;   %1 for absolute uncertainty; 2 for proportional uncertainty (not work currently)
ab910type = 1;      % from [Sato13] for determining the LMI form. See [Sato 13] for detail, 1:9(a),2:9(b),3:10(a),4:10(b). 1 is used for Automatica paper: 
% vareps = 1e-7^2;    % tolerance for expressing matrix inequality constraints using LMIs 
% vareps1 = 1e-7;     % used for BRL condition to make the inequality strict

% line search parameter
Epsi= logspace(-5,2,15);
% Paras for average-dwell-time switching
ta = 100;           %average dwell time
Mu = logspace(1e-2, 2, 11); 
Mu = Mu(end-1);

lmivars0 = [];      % struct variable for storing the initial values for all the optimization variables in LMI Lab; can be set empty
IntuitiveSw = 0;    % intuitive switching controller design without any theoretical guarantee. default:0

% Uncertainty  
deltam = 0 
Delta = [-deltam, deltam];

if XY_PD ~= 0 %has to use gridding method for PDLF
    IsGridding = 1; % 1 for using gridding method to solve the PDLMIs
end
for hiding =1
if GSUncertType == 1 
    delta1_l = -deltam; delta1_u = deltam;
    delta2_l = -deltam; delta2_u = deltam; 
elseif GSUncertType == 2 && SwLogic ~= 0
    delta1_l = min(Theta1(1,2)/(1-lambda_u), Theta1(1,2)/(1+lambda_u))-Theta1(1,2);
    delta1_u = max(Theta1(1,2)/(1-lambda_u), Theta1(1,2)/(1+lambda_u))-Theta1(1,2);
    delta2_l = min(Theta1(2,1)/(1-lambda_u), Theta1(2,1)/(1+lambda_u))-Theta1(2,1);
    delta2_u = max(Theta1(2,1)/(1-lambda_u), Theta1(2,1)/(1+lambda_u))-Theta1(2,1);
end

%% plant definition
Theta1 = [-2 0.2;-0.2 2];
Theta2 = [-2 0.2;-0.2 2];
syms theta1 theta2

NumGain = 1e6; % just for avoiding numerical problem
zeta = 50e-3;% 50e-3;
wn = 100;%m = 1;    
A = [0 1+theta1^2; -wn^2*(1+(0.5+theta1)^2-0.5*theta1) -2*zeta*wn+theta2];%-0.5*theta1
%     A = [0 1; -(wn^2+wn^2/2  *theta1) -2*zeta*wn]
B = [0;1];
C = [1 0];
We = makeweight(2,50,0.1); %makeweight(30,50,0.5)    
Wk = ss(tf([1 1],[1 1]));% Integrator
Wu = 1e-7;%1e-7;%2e-6;

Gplt.A = [A zeros(2); -Wk.b*C Wk.a 0; -We.b*C 0 We.a]; Gplt.B1 = [0;0;Wk.b;We.b]; Gplt.B2 = [B;0;0]*NumGain;
Gplt.C1 = [-We.d*C 0 We.c; 0 0 0 0]; Gplt.D11 =[We.d;0]; Gplt.D12 = [0;Wu]*NumGain;
Gplt.C2 = [-Wk.d*C Wk.c 0]; Gplt.D21 = Wk.d; Gplt.D22 = 0;
A_we =  [-Wk.b*C Wk.a 0 ;-We.b*C 0 We.a];    

%     %% Hinf controller design; for debug
%     Ga = AugPltEv(Gplt,0);
%     a = Ga.A;

%     b = [Ga.B1 Ga.B2];
%     c = [Ga.C1;Ga.C2];
%     d = [Ga.D11 Ga.D12;
%         Ga.D21 Ga.D22];
%     ga = ss(a,b,c,d);    
%     [Kinf,CL,Gam,INFO] = hinfsyn(ga,1,1,'method','lmi');
%     Plt = ss([0 1; -1e4 -2*zeta*100],B,C,0);
%     return;           
gs_num = 2;
F_theta = @(x) [1 x(1) x(2)]; %function for PD matrices
d_F_theta = @(x) [0 1 1];  %function for derivative of PD matrices
FthetaNum = [1 1 1]; %affine form       
d_thetah = (1+10)*2;  
theta_min = Theta1(1,1);
theta_max = Theta1(end,2); 
IsGridding = 1;
if SwLogic == 0
    Theta1 = [Theta1(1,1) Theta1(end,end)];
    if gs_num == 2
        Theta2 = [Theta2(1,1) Theta2(end,end)];
    end
end
n =  size(Gplt.A,1);
nw = size(Gplt.B1,2);
nu = size(Gplt.B2,2);
nz = size(Gplt.C1,1);
ny = size(Gplt.C2,1);  
clear theta1 theta2
theta1 = 0; theta2 = 0;    
A_0 = double(subs(Gplt.A));  
theta1 = 1; theta2 = 0;
A_1 = double(subs(Gplt.A))-A_0;
theta1 = 0; theta2 = 1;
A_2 = double(subs(Gplt.A))-A_0;   
B2 = Gplt.B2;
C2 = Gplt.C2;  
I = eye(n);
end % hiding the code
%  epsi
if deltam == 0 
    Epsi = 0
    d_thetah = d_thetah/11;
else
    Epsi = 3.162e-5; % Epsi(1);
end
d_Thetah = {[-d_thetah d_thetah],[-d_thetah d_thetah]};   

% regnum and ssnum determination,
regnum1 = size(Theta1,1); %num. of subsets for theta1, only partition theta_1
regnum2 = size(Theta2,1); % num. of subsets for theta2;
regnum = regnum1 * regnum2; % Total num. of subsets. S shape from the bottom to order them, an example:
% 4 5 6
% 1 2 3
SSNum = (regnum1-1)*2*regnum2+(regnum2-1)*2*regnum1; % switching surface num
% Determination of regid1 & regid2  
REGID = zeros(regnum,3);
for regid = 1:regnum
   if mod(regid,regnum1) == 0  
        regid1 = regnum1;
    else
        regid1 = mod(regid,regnum1);
   end
   if regid > regnum1
        regid2 = 1+ floor(regid/(regnum1+0.1));
    else
        regid2 = 1;
   end  
   REGID(regid,:) = [regid regid1 regid2];   
end
clear regid1 regid2;
%% grid of the (theta,delta) space for solving PDLMIs.
delta_l= delta1_l; delta_u = delta1_u;
for regid1 = 1:regnum1
    [ThetaT{regid1},DeltaT{regid1}] = AdmRegGrid(Theta1(regid1,:),delta_l,delta_u,theta_min,theta_max,IsGridding,IntuitiveSw); 
end
if gs_num == 2
    for regid2 = 1:regnum2
        [Theta2T{regid2},Delta2T{regid2}] = AdmRegGrid(Theta2(regid2,:),delta_l,delta_u,theta_min,theta_max,IsGridding,IntuitiveSw); 
    end
else       
    Theta2T = {0}; Delta2T = {0};
end
plant_paras.ThetaT = ThetaT;  
plant_paras.DeltaT = DeltaT;
plant_paras.Theta2T = Theta2T; 
plant_paras.Delta2T = Delta2T;
plant_paras.Theta1 = Theta1;
plant_paras.Theta2 = Theta2;
plant_paras.d_Thetah = d_Thetah;
plant_paras.REGID = REGID;
design_paras.SwLogic = SwLogic;
design_paras.ab910type = ab910type;
design_paras.XY_PD = XY_PD;
design_paras.F_theta = F_theta;
design_paras.d_F_theta = d_F_theta;
design_paras.FthetaNum = FthetaNum;

%% Line search for epsi or mu 
% Note that line search will be implemented at most for one of epsi or mu. 
Gam_mu_eps = zeros(length(Mu)+length(Epsi)-1,regnum+3);
for mu_index = 1:length(Mu)
    mu = Mu(mu_index);
    lambda0 = log(mu)/ta + 1e-8; 
    ADTPara = [mu,ta,ASWCond_Ceq];
    design_paras.ADTPara = ADTPara;
    for epsi_index = 1:length(Epsi);% length(Epsi)        
        epsi = Epsi(epsi_index)*ones(1,regnum);  
        design_paras.epsi = epsi;
        [Gam,gam,X,Y,Ah,Bh,Ch,Dh,SolverInfo] = SLPV_uncert(Gplt,plant_paras,design_paras,lmivars0);    
        mu_epsi_index = mu_index+ epsi_index -1;
        Gam_mu_eps(mu_epsi_index,1) = Gam;
        Gam_mu_eps(mu_epsi_index,2) = mu;
        Gam_mu_eps(mu_epsi_index,3:length(epsi)+2) = epsi;  
    end
    ElapsedTime = toc;   
end    
[Gam_opt, index] = min(Gam_mu_eps(:,1))

%% Check the stability of the CL systems 
if SwLogic == 0 && XY_PD == 0 % for LPV control and constant Lyapunov functions
    % specify theta and delta value    
    theta1 = 0; theta2 = 0; delta = [0 0]';
    theta1h = theta1+ delta(1); theta2h = theta2+ delta(2); 
    A = A_0+theta1*A_1+theta2*A_2;
    A_t = A_0+theta1h*A_1+theta2h*A_2;
    Ah_t  = Ah(:,:,1)+theta1h*Ah(:,:,2)+theta2h*Ah(:,:,3);
    Bh_t  = Bh(:,:,1)+theta1h*Bh(:,:,2)+theta2h*Bh(:,:,3);
    Ch_t  = Ch(:,:,1)+theta1h*Ch(:,:,2)+theta2h*Ch(:,:,3);
    Dh_t  = Dh(:,:,1)+theta1h*Dh(:,:,2)+theta2h*Dh(:,:,3);
    N  = -Y ;
    M  = X -I/Y ;

    Dk = Dh_t ;
    Ck = (Ch_t - Dk*C2*X)/M';
    Bk = N\(Bh_t -Y*B2*Dk);
    Ak =  N\(Ah_t -Bh_t*C2*X -...
        Y*B2*Ch_t - Y*(A_t-B2*Dh_t*C2)*X)/M'; 
    A_cl = [A+B2*Dk*C2 B2*Ck; Bk*C2 Ak];
    if max(real(eig(A_cl))) >= 0
        disp('CL system is unstable!');
        eig(A_cl)
        A_cl_cond = cond(A_cl)
    else
%         disp('CL system is stable!');
%         eig(A_cl)
%         A_cl_cond = cond(A_cl)
    end
end