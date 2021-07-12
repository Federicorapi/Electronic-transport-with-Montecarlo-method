clear all
close all
clc
%% Compute the acoustic deformation potential scattering in the Gamma valley of GaAs (inelastic case)
Q = 1.6021766208e-19; % elementary charge, C

nE = 200; % number of energy points
vE = linspace(0,1,nE)*Q; % energy axis, J
Waco_abs = zeros(2,nE);
Waco_emi = zeros(2,nE);
T = [77 300]; % temperature, K
for i = 1:2
    for ie = 1:nE, E = vE(ie); 
    [Waco_emi(i,ie), Waco_abs(i,ie)] = aco_scat_inel(E,T(i));
    end
end

figure(1), hold on
grid on
plot(vE/Q,Waco_emi(1,:),'r--','linewidth',2)
plot(vE/Q,Waco_abs(1,:),'b--','linewidth',2)
plot(vE/Q,Waco_emi(1,:)+Waco_abs(1,:),'k--','linewidth',2)

plot(vE/Q,Waco_emi(2,:),'r.-','linewidth',2)
plot(vE/Q,Waco_abs(2,:),'b.-','linewidth',2)
plot(vE/Q,Waco_emi(2,:)+Waco_abs(2,:),'k.-','linewidth',2)

set(gca,'FontSize',14,'FontName','Arial','box','on','YScale','log')
ylabel('Acoustic scattering rate, 1/s'), xlabel('Energy, eV')
legend('Emission 77K', 'Absorption 77K', 'total 77K','Emission 300K' ,'Absorption 300K' ,'total 300K')
hold off

%% Compute the acoustic deformation potential scattering rate in the Gamma valley (elastic approximation) 
Q = 1.6021766208e-19; % elementary charge, C
iv = 1;
nE = 200; % number of energy points
vE = linspace(0,1,nE)*Q; % energy axis, J
Waco_el_par = zeros(1,nE);
Waco_el_nonpar = zeros(1,nE);
%T = [77 300]; % temperature, K
 %for i = 1:2
    for ie = 1:nE, E = vE(ie); 
    [Waco_el_par(ie),Waco_el_nonpar(ie)] = aco_scat_el(E,300,iv);   %%%300K
    end
 %end
figure(2), hold on 
grid on
plot(vE/Q,Waco_el_par,vE/Q,Waco_el_nonpar,'linewidth',2);
plot(vE/Q,Waco_emi(2,:)+Waco_abs(2,:),'k.-','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on','YScale','log')
ylabel('Acoustic scattering rate, 1/s'), xlabel('Energy, eV')
legend('Parabolic approx', 'Non parabolic approx','Scattering rate inelastic case')
title('Elastic and inelastic approximation @300K')
hold off

%% Compute the acoustic intervalley derfomation potential scattering from Gamma-L
Q = 1.6021766208e-19; % elementary charge, C
nE = 200; % number of energy points
vE = linspace(0,1,nE)*Q; % energy axis, J
T = 300; % temperature, K
iv = 1;
Waco_inter_em = zeros(1,nE);
Waco_inter_abs = zeros(1,nE);

for ie = 1:nE, E = vE(ie); 
    [Waco_inter_em(ie),Waco_inter_abs(ie)] = aco_inter(E,T,iv);   
end

figure(3), hold on
grid on
plot(vE/Q,Waco_inter_em,vE/Q,Waco_inter_abs,'linewidth',2);
set(gca,'FontSize',14,'FontName','Arial','box','on','YScale','log')
ylabel('Intervalley scattering rate, 1/s'), xlabel('Energy, eV')
legend('\Gamma-L emission','\Gamma-L absorption');
title('Acousting scattering intervalley rate @300K')
hold off

%% Compute the polar optical scattering rate without the parabolic approximation
% In order to run this section, the inputs in the pol_scat function must be
% changed according to what is present in line 92

Q = 1.6021766208e-19; % elementary charge, C
nE = 200; % number of energy points
vE = linspace(0,1,nE)*Q; % energy axis, J
iv = 1;
T = 300; % temperature, K
Wpop_emi = zeros(2,nE);
Wpop_abs = zeros(2,nE);
alpha = [0.64/Q,0];
for i = 1:2
    for ie = 1:nE, E = vE(ie); 
        [Wpop_emi(i,ie),Wpop_abs(i,ie)] = pol_scat(E,T,iv,alpha(i)); 
    end
end
figure(4), hold on
grid on
plot(vE/Q,Wpop_emi(1,:),vE/Q,Wpop_abs(1,:),'linewidth',2);
plot(vE/Q,Wpop_emi(1,:)+Wpop_abs(1,:),'linewidth',2);
set(gca,'FontSize',14,'FontName','Arial','box','on','YScale','log')
ylabel('Scattering rate, 1/s'), xlabel('Energy, eV')
legend('emission','absorption','total');
title('    Polar optical scattering rate, non parabolic approx')
hold off

figure(5), hold on
grid on
plot(vE/Q,Wpop_emi(1,:)+Wpop_abs(1,:),'linewidth',2);
plot(vE/Q,Wpop_emi(2,:)+Wpop_abs(2,:),'linewidth',2);
set(gca,'FontSize',14,'FontName','Arial','box','on','YScale','log')
ylabel('Scattering rate, 1/s'), xlabel('Energy, eV')
legend('Non parabolic approx','Parabolic');
title('Polar optical scattering rate')
hold off

%% Compute the impurity scattering rate
Q = 1.6021766208e-19; % elementary charge, C
nE = 200; % number of energy points
vE = linspace(0,0.6,nE)*Q; % energy axis, J
Wimp = zeros(2,nE);
NI = [1e22 1e24];
iv = 1;
T = 300;
for i = 1:2
    for ie = 1:nE, E = vE(ie); 
        [Wimp(i,ie)] = imp_scat(E,NI(i),T,iv);   
    end
end

figure(6), hold on
grid on
plot(vE/Q,Wimp(1,:),'linewidth',2);
plot(vE/Q,Wimp(2,:),'linewidth',2);
set(gca,'FontSize',14,'FontName','Arial','box','on','YScale','log')
ylabel('Scattering rate, 1/s'), xlabel('Energy, eV')
legend('nI = 1e16/cm^3','nI = 1e18/cm^3');
title('Impurity scattering rate')
hold off



