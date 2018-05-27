 function [Gam,gam,X,Y,Ah,Bh,Ch,Dh,SolverInfo] = SLPV_uncert(Gplt,plant_paras,design_paras,lmivars0)
%*********************************************************************************************************************************************
%  
%  Switching LPV controller design under uncertain scheduling parameters 
%             ----------------------------------------
%
% By Pan Zhao, Control Engineering Lab, UBC, BC, Canada.
% Under supervision of Prof. Ryozo Nagamune.
% Creation: May 12, 2015.
% Revision: Sep 14, 2017.
% 
% Function: creating and solving LMI problems involved in the controller design.
%               
% Note: 
%   * Can also be used for exact measurement of scheduling parameter case
%   * The default GS parameter number is 2, for one scheduling parameter
%     case, need to set:
%       - Theta2T = {0}, Delta2T = {0}, d_Thetah ={xx,0},
%       - Fcn_theta is only function of theta1, e.g. Fcn_theta = @(x) [1 x(1)]  
% Modified on 10/07/2015 to use same constant X or Y, so that X(i) = X(j) need not be imposed. In this case, under average dwell time switching,
% constant terms are also imposed to be equal at switching surfaces. The contents for ASWCond_Ceq = 0 does not work now. 

% Gplt, generalized plant, created using Matlab symbols
% ThetaT, cell array, gridded points of theta1 for all subsets, ThetaT{1} for subset 1, ThetaT{2} for subset 2,...
% DeltaT, cell arry, gridded points of delta1 for all subsets
% Theta2T, cell arry, gridded points of theta2 for all subsets
% Delta2T, cell arry, gridded points of delta2 for all subsets
% Theta1, matrix, the jth row representing the edge value of theta1 for subset j
% Theta2, matrix, the jth row representing the edge value of theta2 for subset j
% d_Thetah, cell array, bounds for derivatives of theta, d_Thetah{1} for theta1, d_Thetah{2} for theta2
% REGID, [regid, regid1, regid2]
% SwLogic, 0 for non-switching, 1 for hysteresis switching, 2 for average-dwell-time switching
% ab910type, from [Sato13] for determining the LMI form 
% epsi, from [Sato13]
% mu, ta, for average dewll time switching
% ASWCond_Ceq, for average dwell time switching, 1 for equal X at switching surfaces 
% XY_PD, 0 for constant X and Y, 1 for 
% Fcn_theta: a function handle, Ftheta(theta) will give all the scalar functions for
% the matrix variables, for instance in X = X0+ f1(theta)X1+f2(theta)X2,
% Fcn_theta(theta) = [1 f1(theta) f2(theta)]
% d_Fcn_theta: a function handle for the derivative of the scalar functions
% FthetaNum: a vector, each element show the number of constant matrices and matrices as a function of GS parameters theta1, theta2, ... 
%            for instance, [1 1 1] means one constant matrix, one matrix as a function of theta1, one matrix as a function of theta2,
%            while [1 2] means one constant matrix, two matrices as a function of theta1, no theta2. 
% lmivars0: Initial value for optimization variables. 

%% Extract the parameters  
ThetaT = plant_paras.ThetaT;
DeltaT = plant_paras.DeltaT;
Theta2T = plant_paras.Theta2T;
Delta2T = plant_paras.Delta2T;
Theta1 = plant_paras.Theta1;
Theta2 = plant_paras.Theta2;
d_Thetah = plant_paras.d_Thetah;
REGID = plant_paras.REGID;
SwLogic = design_paras.SwLogic;
ab910type = design_paras.ab910type;
epsi = design_paras.epsi;
ADTPara = design_paras.ADTPara;
XY_PD = design_paras.XY_PD;
F_theta = design_paras.F_theta;
d_F_theta = design_paras.d_F_theta;
FthetaNum = design_paras.FthetaNum;
mu = design_paras.ADTPara(1);
ta = design_paras.ADTPara(2);
ASWCond_Ceq = design_paras.ADTPara(3);
% for avoiding numerical problem
vareps = 1e-5^2; % used for expressing equality constraint
vareps1 = 0; % used for BRL condition and positivity of Lyapunov function

lambda0 = log(mu)/ta + 1e-8;
deltam = max(DeltaT{1});
 
if norm(Theta2T{1}) == 0 %% only one GS para theta1
   clear Delta2T  % in case 
   Delta2T = {0}; 
   GSParaNum = 1; 
   d_Thetah{2} = 0;
   if max(REGID(:,3)) > 1
       disp('Error!REGID is not for one GS parameter!');
       return;
   end   
else
    GSParaNum = 2;
end
if XY_PD == 0
    d_Thetah = {0,0};
    disp('Derivative of Thetah is set to 0 due to use of constant Lyapunov function');
end
% Generalized plant parameter
B2 = Gplt.B2; C2 = Gplt.C2;       
n =  size(Gplt.A,1);
nw = size(Gplt.B1,2);
nu = size(Gplt.B2,2);
nz = size(Gplt.C1,1);
ny = size(Gplt.C2,1); 
I = eye(n); 
% subset parameter 
regnum1 = size(ThetaT,2); %num. of subsets for theta1, only partition theta_1
regnum2 = size(Theta2T,2); % num. of subsets for theta2;
regnum = regnum1 * regnum2; % Total num. of subsets. S shape from the bottom to order them, an example:
% 4 5 6
% 1 2 3
% just for positivity of the Lyapunov function 
if GSParaNum == 1
    theta2gridnum = 1;
else
    theta2gridnum = 10;
end   

%% Optimization problem definition
for nomeaning = 1
lminum = 0;
% Defining opt variables
setlmis([]);
Gam = lmivar(1,[1 1]);
for regid = 1:regnum
    for Id_Ftheta=1:sum(FthetaNum)%
        if XY_PD == 1 %Y0 is constatnt, so define X
           X(Id_Ftheta,regid)=lmivar(1,[n 1]); 
        elseif XY_PD == 2    %X0 is constatnt
           Y(Id_Ftheta,regid)=lmivar(1,[n 1]);
        end
        Ah(Id_Ftheta,regid)=lmivar(2,[n n]);
        Bh(Id_Ftheta,regid)=lmivar(2,[n ny]);
        Ch(Id_Ftheta,regid)=lmivar(2,[nu n]);
        Dh(Id_Ftheta,regid)=lmivar(2,[nu ny]);
    end 
    gam(regid) = lmivar(1,[1 1]);
end   

switch XY_PD 
    case 1%Y0 is constatnt, so define X
        Y =lmivar(1,[n 1]); 
    case 2 %X0 is constatnt
        X =lmivar(1,[n 1]);
    case 0
        X =lmivar(1,[n 1]);
        Y =lmivar(1,[n 1]);
end

%% subset conditons
for regid = 1:regnum           
   % Get the gridding points for the admissible region  
   regid1 = REGID(regid,2); regid2 = REGID(regid,3);
   thetaT = ThetaT{regid1};deltaT = DeltaT{regid1}; theta2T = Theta2T{regid2}; delta2T = Delta2T{regid2};        

   %% Positivity of Lyapunov function
    if XY_PD == 0 % if using PDLF, then there is only one LMI for positivity of LF.
        lminum = lminum + 1;                
        lmiterm([lminum 1 1 X],1,-1);
        lmiterm([lminum 1 1 0],vareps1); % improving numerical property for simulation
        lmiterm([lminum 2 1 0],-1);
        lmiterm([lminum 2 2 Y],1,-1);
        lmiterm([lminum 2 2 0],vareps1);
    else  
        for theta1h  = unique(thetaT) % linspace(min(thetaT),max(thetaT),5) % Gridding may not be necessary because use of affine Laypunov function
            for theta2h = unique(theta2T) % unique(linspace(min(theta2T),max(theta2T),theta2gridnum))%linspace(min(theta2T),max(theta2T),theta2gridnum)     
                thetah = [theta1h;theta2h];
                Fthetah = F_theta(thetah);                    
                lminum = lminum +1;
                 if XY_PD == 1 % PD X                            
                     for Id_Ftheta = 1:length(Fthetah)
                         lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta),-1);
%                                  X = X+X_theta(:,:,Id_Ftheta,regid)*Fthetah(Id_Ftheta); 
                     end
                     lmiterm([lminum 2 2 Y],1,-1);
                elseif XY_PD == 2 % PD Y
                     for Id_Ftheta = 1:length(Fthetah)
                         lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta),-1);
%                                  Y = Y+Y_theta(:,:,Id_Ftheta,regid)*Fthetah(Id_Ftheta); 
                     end   
                     lmiterm([lminum 1 1 X],1,-1);
                 end                         
                 lmiterm([lminum 1 1 0],vareps1); % improving numerical property for simulation
                 lmiterm([lminum 2 1 0],-1);    
                 lmiterm([lminum 2 2 0],vareps1);
%                          [X I; I Y] >= eye(2*n)*1e-6; %% This is fairly important for successful simulation in Simulink. If using 0, simulation may diverge.            
            end
        end
    end 

    for Id_theta1 = 1:length(thetaT)
        theta1 = thetaT(Id_theta1);
        for Id_theta2 = 1:length(theta2T)
            theta2 = theta2T(Id_theta2);  
            delta1 = deltaT(Id_theta1); delta2 = delta2T(Id_theta2);
            theta1h = theta1+delta1; theta2h = theta2+delta2;   
            theta = [theta1;theta2]; thetah = [theta1h;theta2h];

            Ga = AugPltEval(Gplt, theta);
            A = Ga.A;
            B1 = Ga.B1;B2 = Ga.B2;
            C1 = Ga.C1;C2 = Ga.C2;
            D11 = Ga.D11; D12 = Ga.D12;
            D21 = Ga.D21; D22 = Ga.D22;                    
            Ga =  AugPltEval(Gplt, thetah);
            A_h = Ga.A;  
            Fthetah = F_theta(thetah);   
            d_Fthetah = d_F_theta(thetah);             
            for Id_d_thetah1 = 1:length(d_Thetah{1})
                d_thetah1 = d_Thetah{1}(Id_d_thetah1);
                for Id_d_thetah2 = 1:length(d_Thetah{2}) % d_Thetah{2} = 0, for one GS para case
                    d_thetah2 = d_Thetah{2}(Id_d_thetah2);              
                    lminum = lminum+1;                            
                    switch XY_PD
                        case 0
                            lmiterm([lminum 1 1 X],A,1,'s'); %XA+A'X
                            lmiterm([lminum 4 1 X],C1,1);
                            lmiterm([lminum 2 2 Y],1,A,'s');
                            lmiterm([lminum 3 2 Y],B1',1);                                         
                        case 1
                            for Id_Ftheta=1:sum(FthetaNum)
                                lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*A,1,'s');
                                lmiterm([lminum 4 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*C1,1);                                        
                            end
                            lmiterm([lminum 2 2 Y],1,A,'s');
                            lmiterm([lminum 3 2 Y],B1',1);                                        
                            % - d_X
                            for Id_Ftheta = 2:1+FthetaNum(2)
                                lmiterm([lminum 1 1 X(Id_Ftheta,regid)],-d_Fthetah(Id_Ftheta)*d_thetah1,1);
%                                         -d_X = -(d_X+X_theta(:,:,Id_Ftheta,regid)*d_Fthetah(Id_Ftheta)*d_thetah1); 
                            end                                 
                            for Id_Ftheta = 2+FthetaNum(2):sum(FthetaNum)%
                                lmiterm([lminum 1 1 X(Id_Ftheta,regid)],-d_Fthetah(Id_Ftheta)*d_thetah2,1);
                            end
                        case 2
                            for Id_Ftheta = 1:sum(FthetaNum)
                                lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta),A,'s');
                                lmiterm([lminum 3 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*B1',1);  
                            end
                            lmiterm([lminum 1 1 X ],A,1,'s'); %XA+A'X
                            lmiterm([lminum 4 1 X ],C1,1);   
                            % d_Y
                            for Id_Ftheta = 2:1+FthetaNum(2)
                                lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],d_Fthetah(Id_Ftheta)*d_thetah1,1);
%                                         d_X = d_X+X_theta(:,:,Id_Ftheta,regid)*d_Fthetah(Id_Ftheta)*d_thetah1; 
                            end                                 
                            for Id_Ftheta = 2+FthetaNum(2):sum(FthetaNum)% if 
                                lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],d_Fthetah(Id_Ftheta)*d_thetah2,1);
                            end                                    
                    end   
                    lmiterm([lminum 2 1 0],A');
                    lmiterm([lminum 3 1 0],B1');
                    lmiterm([lminum 3 3 gam(regid)],-1,1);
                    lmiterm([lminum 4 2 0],C1);
                    lmiterm([lminum 4 3 0],D11);
                    lmiterm([lminum 4 4 gam(regid)],-1,1);   

                    for Id_Ftheta=1:sum(FthetaNum)   
%                                 regidtemp = regid;
%                                 regid = 1; 
                        lmiterm([lminum 2 1 Ah(Id_Ftheta,regid)],Fthetah(Id_Ftheta),1);                                
                        lmiterm([lminum 2 2 Bh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*1,C2,'s'); 
                        lmiterm([lminum 3 2 -Bh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D21',1);
                        lmiterm([lminum 2 1 -Dh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*C2',B2');                                
                        lmiterm([lminum 3 1 -Dh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D21',B2'); 
                        lmiterm([lminum 4 2 Dh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D12,C2); 
                        lmiterm([lminum 4 3 Dh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D12,D21);     
                        lmiterm([lminum 1 1 Ch(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*B2,1,'s');
                        lmiterm([lminum 4 1 Ch(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D12,1);
%                                 regid = regidtemp;
                    end    
                 % average-dwell-time switching
                    if SwLogic == 2 %average-dwell-time switching, BRL condition should be modified                                
                        switch XY_PD
                            case 0
                                lmiterm([lminum 1 1 X],lambda0,1); %XA+A'X
                                lmiterm([lminum 2 2 Y],lambda0,1); %XA+A'X
                            case 1
                                for Id_Ftheta=1:sum(FthetaNum)
                                    lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*lambda0,1); %XA+A'X
                                end 
                                lmiterm([lminum 2 2 Y],lambda0,1); %XA+A'X                                        
                            case 2
                                lmiterm([lminum 1 1 X],lambda0,1); %XA+A'X
                                for Id_Ftheta=1:sum(FthetaNum)
                                    lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*lambda0,1); %XA+A'X
                                end
                        end
                        lmiterm([lminum 2 1 0],lambda0/2); %XA+A'X
                    end 
                    % for uncertainty in scheduling parameter measurement
                    if deltam ~= 0
                        for nomeaning = 1
                        switch XY_PD
                            case 0
                                switch ab910type
                                    case 1
                                        lmiterm([lminum 5 1 X ],epsi(regid),1);
                                        lmiterm([lminum 5 2 Y ],(A-A_h)',1);
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                    case 2
                                        lmiterm([lminum 5 1 X ],epsi(regid)*(A-A_h),1);                                                
                                        lmiterm([lminum 5 2 Y ],1,1); 
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                    case 3     
                                        lmiterm([lminum 5 1 X ],epsi(regid),1);
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                        lmiterm([lminum 6 2 Y ],(A-A_h)',1);
                                        lmiterm([lminum 6 6 0],-epsi(regid));
                                    case 4
                                        lmiterm([lminum 5 1 X ],epsi(regid)*(A-A_h),1);    
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                        lmiterm([lminum 6 2 Y ],1,1); 
                                        lmiterm([lminum 6 6 0],-epsi(regid)); 
                                end    
                            case 1
                                switch ab910type
                                    case 1
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 5 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*epsi(regid),1); %XA+A'X
                                        end
                                        lmiterm([lminum 5 2 Y ],(A-A_h)',1);
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                    case 2
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 5 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*epsi(regid)*(A-A_h),1); %XA+A'X
                                        end                                             
                                        lmiterm([lminum 5 2 Y ],1,1); 
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                    case 3     
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 5 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*epsi(regid),1); %XA+A'X
                                        end
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                        lmiterm([lminum 6 2 Y ],(A-A_h)',1);
                                        lmiterm([lminum 6 6 0],-epsi(regid));
                                    case 4
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 5 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*epsi(regid)*(A-A_h),1); %XA+A'X
                                        end 
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                        lmiterm([lminum 6 2 Y ],1,1); 
                                        lmiterm([lminum 6 6 0],-epsi(regid)); 
                                end  
                            case 2
                                switch ab910type
                                    case 1
                                        lmiterm([lminum 5 1 X ],epsi(regid),1);
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 5 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*(A-A_h)',1); %XA+A'X
                                        end
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                    case 2
                                        lmiterm([lminum 5 1 X ],epsi(regid)*(A-A_h),1);    
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 5 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta),1); %XA+A'X
                                        end
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                    case 3     
                                        lmiterm([lminum 5 1 X ],epsi(regid),1);
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 6 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*(A-A_h)',1); %XA+A'X
                                        end
                                        lmiterm([lminum 6 6 0],-epsi(regid));
                                    case 4
                                        lmiterm([lminum 5 1 X ],epsi(regid)*(A-A_h),1);    
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 6 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta),1); %XA+A'X
                                        end
                                        lmiterm([lminum 6 6 0],-epsi(regid)); 
                                end   
                        end   
                        end
                    end %deltam                             
                end %Id_d_thetah2
            end %Id_d_thetah2
        end % Id_theta2
    end % Id_theta1            
    lminum = lminum + 1;
    lmiterm([lminum 1 1 gam(regid)],1,1);
    lmiterm([lminum 1 1 Gam],-1,1);
%             regid = regidtemp;
end %regid  
    % An alternative cost function
%         lminum = lminum+1;
%         for regid = 1: regnum            
%             lmiterm([lminum 1 1 gam(regid)],1/regnum,1);            
%         end
%         lmiterm([lminum 1 1 Gam],-1,1);

%% Switching surface conditions, ie. LMIs for the monotonic property of the Lyapunov function   
lmisys = getlmis;
lminum_sub = lminum; 
%         disp(['The   number of lmis for subsets is ', num2str(lminum_sub)]);
% switching surface index increases in the theta1 direction first
% in the way of 
% 4 5 6
% 1 2 3
% and then in theta2 direction in the way of 
% 2 4
% 1 3   
if SwLogic ~= 0
    ss_para.XY_PD = XY_PD;
    ss_para.ASWCond_Ceq = ASWCond_Ceq;
    ss_para.SwLogic = SwLogic;
    ss_para.mu =  mu;
    % for switching surfaces in theta1 direction     
    %  Note that to express X == Y in LMI lab, -[varepsi*E=Id_Ftheta X-Y; X-Y Id_Ftheta]<=0 is used
    for regid2 = 1:regnum2
        RegId1 = REGID(find(REGID(:,3)==regid2),2);
        RegId1 = RegId1'; RegId1 = sort(RegId1);     
        RegId = (regid2-1)*regnum1+RegId1;
        for id1 = RegId1(1:end-1)  
            for theta2 = unique(Theta2(regid2,:))%becasue inequalities is affine w.r.t theta1 & theta2, so check of vertices is enough
                    theta1 = Theta1(id1,end);                                                      
                    Fthetah = F_theta([theta1;theta2]);   
                    jk = [RegId(id1) RegId(id1+1)];
                    [lmisys,lminum] = LMI_SwSurf(lmisys,X,Y,Fthetah,jk,ss_para,lminum);

                    theta1 = Theta1(id1+1,1);
                    Fthetah = F_theta([theta1;theta2]);                                   
                    jk = [RegId(id1+1) RegId(id1)];                            
                    [lmisys,lminum] = LMI_SwSurf(lmisys,X,Y,Fthetah,jk,ss_para,lminum);
            end   
        end %id1
    end % regid2       

    % for switching surfaces in theta2 direction 
    for regid1 = 1:regnum1
        RegId2 = REGID(find(REGID(:,2)==regid1),3);
        RegId2 = RegId2';RegId2 = sort(RegId2);       
        RegId = (RegId2-1)*regnum1+regid1;
        for id2 = RegId2(1:end-1)
            for theta1 = Theta1(regid1,:) %becasue inequalities is affine w.r.t theta1 & theta2, so check of vertices are enough
                theta2 = Theta2(id2,end);
                Fthetah = F_theta([theta1;theta2]);   
                jk = [RegId(id2) RegId(id2+1)];
                [lmisys,lminum] = LMI_SwSurf(lmisys,X,Y,Fthetah,jk,ss_para,lminum);   

                theta2 = Theta2(id2+1,1);
                Fthetah = F_theta([theta1;theta2]);  
                jk = [RegId(id2+1) RegId(id2)];      
                [lmisys,lminum] = LMI_SwSurf(lmisys,X,Y,Fthetah,jk,ss_para,lminum);  
            end    
        end %id2
    end       
end %if 

% formulating the cost function
nvar = decnbr(lmisys); %number of decision variable.
c = zeros(nvar,1);
c(1)=1; 
% LMI Lab solver setting
options(1)= 1e-4; % relative accuary on the optimal value
options(2)= 1000; %Number of Iteration
options(4) =  10; % J, the code terminates when the objective has not decreased by more than the desired relative accuracy during the last J iterations
options(5)= 0; % 1 for not showing the process
% solve the optimization problem
if ~isempty(lmivars0)  
    % initial value setting (optional)
    for hide = 1
    if regnum == 1 && XY_PD == 1 && sum(FthetaNum) == 3                 
         xinit = mat2dec(lmisys,lmivars0.Gam,...
            lmivars0.X(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
            lmivars0.X(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
            lmivars0.X(:,:,3),lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
            lmivars0.Gam, ...
            lmivars0.Y);
    elseif regnum == 2
        switch XY_PD
            case 1
                if sum(FthetaNum) == 2
                    xinit = mat2dec(lmisys,lmivars0.Gam,....
                        lmivars0.X(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.X(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Gam,...
                        lmivars0.X(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.X(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Gam,...
                        lmivars0.Y);
                elseif sum(FthetaNum) == 3
                     xinit = mat2dec(lmisys,lmivars0.Gam,...
                        lmivars0.X(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.X(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.X(:,:,3),lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
                        lmivars0.Gam, ...
                        lmivars0.X(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.X(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.X(:,:,3),lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
                        lmivars0.Gam,...
                        lmivars0.Y);
                end                    
            case 2
                if sum(FthetaNum) == 2
                    xinit = mat2dec(lmisys,lmivars0.Gam,...
                        lmivars0.Y(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.Y(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Gam,...
                        lmivars0.Y(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.Y(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Gam,...
                        lmivars0.X);                        
                elseif sum(FthetaNum) == 3
                    xinit = mat2dec(lmisys,lmivars0.Gam,...
                        lmivars0.Y(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.Y(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Y(:,:,3),lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
                        lmivars0.Gam, ...
                        lmivars0.Y(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.Y(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Y(:,:,3),lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
                        lmivars0.Gam,...
                        lmivars0.X);
                end
            case 0
                if sum(FthetaNum) == 2
                    xinit = mat2dec(lmisys,lmivars0.Gam,...
                        lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Gam,...
                        lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Gam,...
                        lmivars0.Y,lmivars0.Y);
                elseif sum(FthetaNum) == 3
                     xinit = mat2dec(lmisys,lmivars0.Gam,...
                        lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
                        lmivars0.Gam,...
                        lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
                        lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
                        lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
                        lmivars0.Gam,...
                        lmivars0.X,lmivars0.Y);
                end
        end
    elseif regnum == 4
        if sum(FthetaNum) == 3 && XY_PD == 1
        xinit = mat2dec(lmisys,lmivars0.Gam,...
        lmivars0.X(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
        lmivars0.X(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
        lmivars0.X(:,:,3),lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
        lmivars0.Gam, ...
        lmivars0.X(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
        lmivars0.X(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
        lmivars0.X(:,:,3),lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
        lmivars0.Gam, ...
        lmivars0.X(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
        lmivars0.X(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
        lmivars0.X(:,:,3),lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
        lmivars0.Gam, ...
        lmivars0.X(:,:,1),lmivars0.Ah(:,:,1),lmivars0.Bh(:,:,1),lmivars0.Ch(:,:,1),lmivars0.Dh(:,:,1),...
        lmivars0.X(:,:,2),lmivars0.Ah(:,:,2),lmivars0.Bh(:,:,2),lmivars0.Ch(:,:,2),lmivars0.Dh(:,:,2),...
        lmivars0.X(:,:,3),lmivars0.Ah(:,:,3),lmivars0.Bh(:,:,3),lmivars0.Ch(:,:,3),lmivars0.Dh(:,:,3),...
        lmivars0.Gam,...
        lmivars0.Y);
        end  
    end
    [Gam,xopt] = mincx(lmisys,c,options,xinit);
    end % hide
else
    [Gam,xopt] = mincx(lmisys,c,options);
end

% check whether all constraints are satisfied
% lmifail = 0;
%         evals = evallmi(lmisys,xopt);        
%         for i = 1: lminum
%             [lhs,rhs] = showlmi(evals,i);
%             if max(real(eig(lhs-rhs))) > 5e-7 
%                 eig_max = max(real(eig(lhs-rhs)))
%                 disp ('Not all LMI constranits are satisfied')
%                 lmifail = 1;
% %                 break;
%             end
%         end 
%         lmifail
if ~isempty(Gam)        
for regid = 1:regnum
    for Id_Ftheta = 1:sum(FthetaNum)
        Akh(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Ah(Id_Ftheta,regid));
        Bkh(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Bh(Id_Ftheta,regid));
        Ckh(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Ch(Id_Ftheta,regid));
        Dkh(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Dh(Id_Ftheta,regid));
        switch XY_PD
            case 1
                Xk(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,X(Id_Ftheta,regid));
            case 2
                Yk(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Y(Id_Ftheta,regid)); 
            case 3
                Xk(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,X(Id_Ftheta,regid));
                Yk(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Y(Id_Ftheta,regid)); 
        end
    end
    switch XY_PD 
        case 1
            Yk(:,:,regid) = dec2mat(lmisys,xopt,Y);
        case 2
            Xk(:,:,regid) = dec2mat(lmisys,xopt,X);
        case 0 
            Xk(:,:,regid) = dec2mat(lmisys,xopt,X);
            Yk(:,:,regid) = dec2mat(lmisys,xopt,Y);
    end
    gam_reg(regid)  = dec2mat(lmisys,xopt,gam(regid));
end   
Ah = Akh;
Bh = Bkh;
Ch = Ckh;
Dh = Dkh;        
X = Xk;
Y = Yk;    
gam = gam_reg;
clear Akh Bkh Ckh Dkh Xk Yk  
else
    Ah = [];
    Bh = [];
    Ch = [];
    Dh = [];
    Y = [];
    X = [];
    gam = Inf;
    Gam = Inf;
    lmifail = 1;
end
end % nomeaning, just for hiding the code
SolverInfo.lmifail = lmifail;
