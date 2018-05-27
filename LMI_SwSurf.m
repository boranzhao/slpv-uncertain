function [lmisys,lminum] = LMI_SwSurf(lmisys0,X,Y,Fthetah,jk,ss_para,lminum)
%% Generating LMIs at switching surfaces 
XY_PD = ss_para.XY_PD;
ASWCond_Ceq = ss_para.ASWCond_Ceq;
SwLogic = ss_para.SwLogic;
mu = ss_para.mu;

j = jk(1); k = jk(2);
setlmis(lmisys0)
switch XY_PD 
    case 0 %for constant LF, Lyapunov functions have to be same for all subsets. 
        if SwLogic == 1
        elseif SwLogic == 2   %% not exactly right.
            if ASWCond_Ceq == 1 
               lminum = lminum + 1;
               lmiterm([lminum 1 1 Y(j)],-mu, 1);
               lmiterm([lminum 1 1 Y(r)],1, 1);
               lmiterm([lminum 2 1 0],-sqrt(mu-1));
               lmiterm([lminum 2 2 X(j)],-1, 1);  
%                [mu*Y(j)-Y(r) sqrt(mu-1)*I;
%                                   sqrt(mu-1)*I X(j)] >= 0; 
               
            elseif ASWCond_Ceq == 0 % Not applicable at this moment
%                 lminum = lminum + 1;
%                 lmiterm([lminum 1 1 0],-vareps);
%                 lmiterm([lminum 2 1 Y(j)],-1,1);
%                 lmiterm([lminum 2 1 Y(k)],1,1);
%                 lmiterm([lminum 2 2 0],-1);                              
            end                                               
        end
    case 1  % PD X, so constant Y should  
        if SwLogic == 1                                          
                %Xj > = Xk 
                lminum = lminum + 1;
                for Id_Ftheta = 1:length(Fthetah)
                    lmiterm([lminum 1 1 X(Id_Ftheta,j)],Fthetah(Id_Ftheta), 1);
                    lmiterm([lminum 1 1 X(Id_Ftheta,k)],-Fthetah(Id_Ftheta), 1);
                end                                  
        elseif SwLogic == 2                               
            if ASWCond_Ceq == 0 % not applicable at this moment   
                lminum = lminum + 1;
%                    lmiterm([lminum 1 1 0],-vareps);
%                    for Id_Ftheta = 1:length(Fthetah)
%                         lmiterm([lminum 1 1 X(Id_Ftheta,j)],-Fthetah(Id_Ftheta), 1);
%                         lmiterm([lminum 1 1 X(Id_Ftheta,k)],Fthetah(Id_Ftheta), 1);
%                    end  
%                    lmiterm([lminum 2 2 0],-1);        
% 
%                    lminum = lminum + 1;                        
%                    lmiterm([lminum 1 1 Y(j)],-mu, 1);
%                    lmiterm([lminum 1 1 Y(k)],1, 1);
%                    lmiterm([lminum 2 1 0],-sqrt(mu-1));
%                    for Id_Ftheta = 1:length(Fthetah)
%                         lmiterm([lminum 2 2 X(Id_Ftheta,j)],-Fthetah(Id_Ftheta), 1);
%                    end  
            elseif ASWCond_Ceq == 1 %  
%                Yj = Yk implicitly imposed due to use of the same lmivar.
                lminum = lminum + 1;   
                for Id_Ftheta = 1:length(Fthetah)
                    lmiterm([lminum 1 1 X(Id_Ftheta,k)],-mu*Fthetah(Id_Ftheta), 1);
                    lmiterm([lminum 1 1 X(Id_Ftheta,j)],Fthetah(Id_Ftheta), 1);
                end  
                lmiterm([lminum 2 1 0],-sqrt(mu-1));
                lmiterm([lminum 2 2 Y],-1, 1);   %Y(j)                                   
%                                     [mu*Xk-Xj sqrt(mu-1)*I;
%                                       sqrt(mu-1)*I Yj] >= 0;    
            end  
        end 
    case 2 % PD Y,
        if SwLogic == 1  % constant X should be equal at SS           
            lminum = lminum + 1;
            for Id_Ftheta = 1:length(Fthetah)
                lmiterm([lminum 1 1 Y(Id_Ftheta,j)],-Fthetah(Id_Ftheta), 1);
                lmiterm([lminum 1 1 Y(Id_Ftheta,k)],Fthetah(Id_Ftheta), 1);
            end                        
        elseif SwLogic == 2                                                      
            if ASWCond_Ceq == 1
%              Xj ==  Xk implicitly imposed due to use of same lmivar
                lminum = lminum + 1;               
                for Id_Ftheta = 1:length(Fthetah)
                    lmiterm([lminum 1 1 Y(Id_Ftheta,j)],-mu*Fthetah(Id_Ftheta), 1);
                    lmiterm([lminum 1 1 Y(Id_Ftheta,k)],Fthetah(Id_Ftheta), 1);
                end  
                lmiterm([lminum 2 1 0],-sqrt(mu-1));
                lmiterm([lminum 2 2 X],-1, 1); % instead of X(j), note that X is the same for subsets
%                                    [mu*Y-Y1 sqrt(mu-1)*I;
%                                       sqrt(mu-1)*I X_temp(j)] >= 0;    
            elseif ASWCond_Ceq == 0 % not work
%                 lminum = lminum + 1;
%                 lmiterm([lminum 1 1 0],-vareps);
%                 for Id_Ftheta = 1:length(Fthetah)
%                     lmiterm([lminum 1 1 Y(Id_Ftheta,j)],-Fthetah(Id_Ftheta), 1);
%                     lmiterm([lminum 1 1 Y(Id_Ftheta,k)],Fthetah(Id_Ftheta), 1);
%                 end  
%                 lmiterm([lminum 2 2 0],-1);  
% %                                     Y == Y1;
%                 lminum = lminum + 1;                        
%                 lmiterm([lminum 1 1 X(k)],-mu, 1);
%                 lmiterm([lminum 1 1 X(j)],1, 1);
%                 lmiterm([lminum 2 1 0],-sqrt(mu-1));
%                 for Id_Ftheta = 1:length(Fthetah)
%                     lmiterm([lminum 2 2 Y(Id_Ftheta,j)],-Fthetah(Id_Ftheta), 1);
%                 end                                   
            end              
        end
end % switch    
lmisys = getlmis;
end


