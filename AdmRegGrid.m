%*************************************************************************
%  
%  Switching LPV controller design under uncertain scheduling parameters 
%             ----------------------------------------
%
% By Pan Zhao, Control Engineering Lab, UBC, BC, Canada.
% Under supervision of Prof. Ryozo Nagamune.
% Creation: May 12, 2015.
 
% Function: grid the admissible region Phi_i^[j];    

function [thetaT,deltaT]  = AdmRegGrid(Theta, delta_l,delta_u,theta_min,theta_max,IsGridding,IntuitiveSw)
if length(Theta) < 2
    thetaT = 0;
    deltaT = 0;
    return;
end

if abs(delta_l)+ abs(delta_u) ~= 0       
    if Theta(1) == theta_min
        thetaT = [Theta(1) Theta(1) Theta(1)-delta_l];
        deltaT = [delta_u 0 delta_l];
    else
        if IntuitiveSw == 0
            thetaT = [Theta(1)-delta_u  Theta(1)-delta_l];
            deltaT = [delta_u delta_l];
         elseif  IntuitiveSw == 1
             thetaT = [Theta(1) Theta(1)];
             deltaT = [delta_u delta_l];
        end                     
    end     
    if Theta(2) == theta_max            
        thetaT = [thetaT Theta(2) Theta(2) Theta(2)-delta_u];
        deltaT = [deltaT delta_l 0 delta_u];
    else
        if IntuitiveSw == 0
            thetaT = [thetaT Theta(2)-delta_l  Theta(2)-delta_u];
            deltaT = [deltaT delta_l delta_u];
        elseif  IntuitiveSw == 1
            thetaT = [thetaT Theta(2)  Theta(2)];
            deltaT = [deltaT delta_l delta_u];
        end
    end 
    %% for gridding
    if IsGridding == 1
%         %% gridding
%         for i = 1:length(thetaT)
%             i_1 = i+1;
%             if i_1 > length(thetaT)
%                 i_1 = 1;
%             end
%             theta_new(i) = (thetaT(i)+thetaT(i_1))/2; %add a point in the middle of any two vertices
%             delta_new(i) = (deltaT(i)+deltaT(i_1))/2;                
%         end
%         theta_new(i+1) = sum(thetaT)/length(thetaT); % add a point in the middle of the region
%         delta_new(i+1) = sum(deltaT)/length(thetaT);         
        theta_new = sum(thetaT)/length(thetaT); % add a point in the middle of the region
        delta_new = sum(deltaT)/length(thetaT); 
        thetaT = [thetaT theta_new];deltaT = [deltaT delta_new];
        clear theta_new delta_new
    end    
else
    if  IsGridding == 1 
        thetaT = linspace(Theta(1),Theta(2),5);
        deltaT = thetaT*0;
    else 
        thetaT = [Theta(1)  Theta(2)];
        deltaT =[0 0];
    end
end        

