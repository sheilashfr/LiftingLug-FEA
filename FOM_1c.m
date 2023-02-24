 clc
clear
%material properties and given values
E=70*10^3;
v=0.3;
yield_stress=400;
K_c=20;
theta_a=0;
theta_b=45;
theta_c=90;
theta_p=30;
strain_a=fn_strain(1944637,theta_a);
strain_b=fn_strain(1944637,theta_b);
strain_c=fn_strain(1944637,theta_c);
strain_p=fn_strain(1944637,theta_p);
p=-400:1:400; %range of sigma p

%calculate stress
sig_xx=(E/(1-v^2))*(strain_a+v*strain_c);
sig_yy=(E/(1-v^2))*(strain_c+v*strain_a);
sig_xy=strain_b*(E/(1-v));

%original stress tensor
sig_ij=[sig_xx sig_xy; sig_xy sig_yy];

%rotate stress tensor
Q=[cos(theta_p) sin(theta_p);-sin(theta_p) cos(theta_p)];
sig_30=Q*sig_ij*transpose(Q);

%find principal stresses
for j=1:length(p);
    sig_p1(j)=((sig_30(1,1)+sig_30(2,2)+p(j))/2)+sqrt(((sig_30(1,1)+p(j)-sig_30(2,2))/2)^2+sig_30(1,2)^2);
    sig_p2(j)=((sig_30(1,1)+sig_30(2,2)+p(j))/2)-sqrt(((sig_30(1,1)+p(j)-sig_30(2,2))/2)^2+sig_30(1,2)^2);
end

%von mises' stress
% for j=1:length(p);
% VM(j)=(1/sqrt(2))*sqrt(((sig_p1(j)-sig_p2(j))^2)+((sig_p2(j))^2)+((sig_p1(j))^2))
% end

for j=1:801;
VM(j)=(1/sqrt(2))*sqrt(((sig_p1(j)-sig_p2(j))^2)+((sig_p2(j))^2)+((sig_p1(j))^2))
end 
% 
vm_400=find(VM<400);
safe_vm1= vm_400 - 401

%vm yield criterion 
safe_vm=[]
for j=1:801
if VM(j)<400
    safe_vm(j,1)=VM(j);
    safe_vm(j,2)=j-401;
else
end
end

%tresca stress
for j=1:length(p);
    av1=abs(sig_p1(j)-sig_p2(j));
    av2=abs(sig_p2(j));
    av3=abs(-sig_p1(j));
   tres=[av1 av2 av3];
   tresca(j)=max(tres);
end

%tresca yield criterion
safe_tres=[]
for j=1:801
if tresca(j)<400
    safe_tres(j,1)=tresca(j);
    safe_tres(j,2)=j-401;
else
end
end

SF_VM=[];
%safety factor
for j=1:length(p)
SF_VM(j)=yield_stress/VM(j);
SF_Tres(j)=yield_stress/tresca(j);
end
maxtres=max(SF_Tres);
maxvm=max(SF_VM);

% %plot safety factor
% plot(p,SF_VM);
% xlabel('Stress (MPa)');
% ylabel("Von Mises' Safety Factor")
% figure
% plot(p,SF_Tres);
% xlabel('Stress (MPa)');
% ylabel("Tresca Safety Factor")


%% crack analysis

%rotate stress tensor
sig_new = [];
a=0.0025;
sqr=sqrt(pi*a);

for p = -400:1:400
    i = p + 401;    
    sig_p{i} = [sig_30(1,1) + p, sig_30(1,2); sig_30(2,1), sig_30(2,2)];
    
    for t = 1:1:181
         rotate{t}= [cosd(t-1), sind(t-1); -sind(t-1), cosd(t-1)];
         rotatetensor{i,t} = rotate{t}*sig_p{i}*transpose(rotate{t}); %%transformed sensor in every orientation for all sigma p values
         % sy(i,t) = rotatetensor{i,t}(2,2);
    end 
end

k=1;
% %obtain max k for every sigma p values
% for i = 1:801
%     for j = 1:181
%         max_k(i,j) = rotatetensor{i,j}(2,2)*sqr; %max k values for each orientation (xi-1) each sigma p (yj-1)
%     end
% end
%% 
%LEFM K Values
for i = 1:801
    for j = 1:181

        max_k(i,j) = rotatetensor{i,j}(2,2)*sqr; %max k values for each orientation (xi-1) each sigma p (yj-1)        
     %Validity: if not valid, then save coordinate
    end
end

ktens=[];
for i = 1:801
    for j = 1:181
        if max_k(i,j) <= 0
            ktens(i,j) = 0;
        elseif max_k(i,j) > 0
            ktens(i,j) = max_k(i,j);
        end
    end
end

% 2D contour plot
sig_p_plot = [-400:400];
orien_angle = [0:1:180];
[X,Y] = meshgrid(orien_angle,sig_p_plot);
figure
c = contourf(X,Y,ktens,[20 20],'ShowText',true);
hold on
a = colorbar
xlabel('Crack orientation (degrees)');
ylabel('Additional stress "P" (MPa)');
a.Label.String = 'Stress Intensity Factor (MPa m^0.5)';
hold off
% 

%% LEFM Validity

LEFM_validity = [];
for i = 1:801
    for j = 1:181
        LEFM_validity(i,j) = 4/pi.*(max_k(i,j)/400)^2;
    end
end

%tensile k value
ktens=[];
for i = 1:801
    for j = 1:181
        if LEFM_validity(i,j) <= 0
            ktens(i,j) = 0;
        elseif max_k(i,j) > 0
            ktens(i,j) = LEFM_validity(i,j);
        end
    end
end

% % 2D contour plot
% sig_p_plot = [-400:400];
% orien_angle = [0:1:180];
% [X,Y] = meshgrid(orien_angle,sig_p_plot);
% figure
% c = contourf(X,Y,LEFM_validity,'ShowText',true);
% hold on
% a = colorbar
% xlabel('Crack orientation (degrees)');
% ylabel('Additional stress "P" (MPa)');
% a.Label.String = 'Stress Intensity Factor (MPa m^0.5)';
% hold off
% % 



%% EPFM Check
a = 0.0025;
for i = 1:801
    for j = 1:181
        if a < LEFM_validity(i,j)
            max_k(i,j) = max_k(i,j)/sqrt(1-(1/2)*((rotatetensor{i,j}(2,2))/400)^2);
        elseif a >= LEFM_validity(i,j)
            max_k(i,j) = max_k(i,j);
        end
    end
end

%%
        

%tensile k value
ktens=[];
for i = 1:801
    for j = 1:181
        if max_k(i,j) <= 0
            ktens(i,j) = 0;
        elseif max_k(i,j) > 0
            ktens(i,j) = max_k(i,j);
        end
    end
end


%tensile k value
ktens=[];
for i = 1:801
    for j = 1:181
        if max_k(i,j) <= 0
            ktens(i,j) = 0;
        elseif max_k(i,j) > 0
            ktens(i,j) = max_k(i,j);
        end
    end
end

% 2D contour plot
sig_p_plot = [-400:400];
orien_angle = [0:1:180];
[X,Y] = meshgrid(orien_angle,sig_p_plot);
figure
c = contourf(X,Y,ktens,[20 20],'ShowText',true);
hold on
a = colorbar
xlabel('Crack orientation (degrees)');
ylabel('Additional stress "P" (MPa)');
a.Label.String = 'Stress Intensity Factor (MPa m^0.5)';
hold off
% 


% % 2D contour plot
% sig_p_plot = [-400:400];
% orien_angle = [0:1:180];
% [X,Y] = meshgrid(orien_angle,sig_p_plot);
% figure
% c = contourf(X,Y,ktens,[20 20],'ShowText',true);
% hold on
% a = colorbar
% xlabel('Crack orientation (degrees)');
% ylabel('Additional stress "P" (MPa)');
% a.Label.String = 'Stress Intensity Factor (MPa m^0.5)';
% hold off
% % 


% %% Plots
% mesh(max_k)
% xlabel('Angle(degrees)')
% ylabel('sigma p (MPa)')
% zlabel ('K')

