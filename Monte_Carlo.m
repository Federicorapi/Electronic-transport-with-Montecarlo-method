%% Montecarlo simulation
clear all
close all
clc
tic
Q = 1.6021766208e-19; % elementary charge, C
T = 300; % temperature, K
tau_sim = 1e-9; % simulation time, s
tau = 1e-16; % flight time, s
Fx = [1e4 1e5 3e5 4e5 5e5 7e5 9e5] ; % electric field along the x-axis, V/m (0.1 kV/cm)
time = 0:tau:tau_sim; % time axis, s
nt = length(time); % # of time steps
E0 = [0 0.32]*Q; % energy reference of G and L valleys, J 
nI = 10^20; % impurity concentration, 1/m^3
HBAR = 1.054e-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
eps0 = 8.854e-12; % vacuum permittivity constant, F/m
eps_s = 12.9*eps0; % static dielectric constant, F/m
q0 = sqrt(Q^2*nI/(eps_s*kB*T));
M0 = 9.1095e-31; % electron mass, kg
meff_G = 0.067*M0; %kg
ml = 1.2*M0; % kg
mt = 0.2*M0; %kg
meff_L = (ml*mt^2)^(1/3); %kg
Delta = 0.3*Q; % difference of energy minima, J
int_phon_en = 28e-3*Q; % intervalley phonon energy, J
N_field = length(Fx); % number of electric fields used


% energy dispersion relation, (kx,ky,kz) -> E
k2energy = @(kx,ky,kz,iv) (iv==1)*HBAR^2*(kx^2+ky^2+kz^2)/(2*meff_G) + (iv==2)*HBAR^2*(kx^2+ky^2+kz^2)/(2*meff_L);

%
% inverse dispersion relation E -> norm(kx,ky,kz)
energy2k = @(E,iv) (iv==1)*sqrt(2*meff_G*E/HBAR^2) + (iv==2)*sqrt(2*meff_L*E/HBAR^2);

% initial state 
iv = 1; % valley index (1 -> G, 2 -> L)
kx = 0; ky = 0; kz = 0; % wavevector, 1/m
knorm = 0;
%
% statistical data 
vEkini = zeros(N_field,nt); % array of kinetic energies before drift 
vEkinf = zeros(N_field,nt); % array of kinetic energies after drift
vEkins = zeros(N_field,nt); % array of kinetic energy right after scattering
vEtoti = zeros(N_field,nt); % array of total energies before drift
vEtotf = zeros(N_field,nt); % array of total energies after drift
E_av = zeros(N_field,nt); % array of avarage kinetic energy at time t
Vel_av = zeros(N_field,nt); % array of mean velocity during drift
vEkin_tau = zeros(N_field,nt);
vDeltaEkin = zeros(N_field,nt);
gamma_count = zeros(1,N_field); % count number of time in which the particle is in Gamma
l_count = zeros(1,N_field); % count number of time in which the particle is in L
vel = zeros(N_field,nt);
%
% Loop for different fields------------------------------------------
for j = 1 : N_field
% Monte Carlo main loop ---------------------------------------------
for ii = 1:nt % time loop 
t = time(ii); % current time
%
% update energy
Ekini = k2energy(kx,ky,kz,iv); % kinetic energy before drift, J
Etoti = Ekini + E0(iv); % total energy (kinetic + potential), J 
%
% drift particle in momentum space
kx = kx - Q*Fx(j)*tau/HBAR;
Ekinf = k2energy(kx,ky,kz,iv); % kinetic energy after drift, J
Etotf = Ekinf + E0(iv); % total energy (kinetic + potential), J 
% drift particle in real space
vel(j,ii) = -(Ekinf-Ekini)/(Q*Fx(j)*tau); % instantaneous velocity, m/s 
% collect statistics
vEkini(j,ii) = Ekini;
vEkinf(j,ii) = Ekinf;
vEtoti(j,ii) = Etoti;
vEtotf(j,ii) = Etotf;
% compute scattering rates
[Wpop_em,Wpop_ab] = pol_scat(Ekinf,T,iv); 
[Waco_el_par, Waco_el_nonpar] = aco_scat_el(Ekinf,T,iv);
[Wimp]            = imp_scat(Ekinf,nI,T,iv);
[Wiv_em, Wiv_ab] = aco_inter(Ekinf,T,iv);
W = [Wpop_em Wpop_ab Waco_el_par Wimp Wiv_em Wiv_ab]; % scattering rates, 1/s
% choose scattering mechanism and final state
L = cumsum(W);
r = rand; % pick a random number
%
if(r < L(1)*tau) % select polar optical, emission
    
Ekins = Ekinf - hwpop; % kinetic energy after scattering, J
knorm = energy2k(Ekins,iv);
% compute scattering angles theta and phi
r_theta = rand;
r_phi = rand;
phi = 2*pi*r_phi;
f = (2*sqrt(Ekins*Ekinf))/(sqrt(Ekinf)-sqrt(Ekins))^2;
theta = acos(((1+f)-(1+2*f)^r_theta)/f);
% rotate back to laboratory frame (anisotropic scattering only)
Rx = [ ky/sqrt(kx^2+ky^2),kx*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)), kx/sqrt(kx^2+ky^2+kz^2)];
Ry = [-kx/sqrt(kx^2+ky^2),ky*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)),ky/sqrt(kx^2+ky^2+kz^2)];
Rz = [ 0,-sqrt(kx^2+ky^2)/sqrt(kx^2+ky^2+kz^2),kz/sqrt(kx^2+ky^2+kz^2)];          
kve = knorm*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];      
kx = Rx * kve; ky = Ry * kve; kz = Rz * kve;              
%      
elseif(r < L(2)*tau) % select polar optical, absorption   
Ekins = Ekinf + hwpop; % kinetic energy after scattering, J
knorm = energy2k(Ekins,iv);
% scattering angles theta and phi
r_theta = rand;
r_phi = rand;
phi = 2*pi*r_phi;
f = (2*sqrt(Ekins*Ekinf))/(sqrt(Ekinf)-sqrt(Ekins))^2;
theta = acos(((1+f)-(1+2*f)^r_theta)/f);
% rotate back to laboratory frame (anisotropic scattering only)
Rx = [ ky/sqrt(kx^2+ky^2),kx*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)),kx/sqrt(kx^2+ky^2+kz^2)];
Ry = [-kx/sqrt(kx^2+ky^2), ky*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)), ky/sqrt(kx^2+ky^2+kz^2)];
Rz = [ 0, -sqrt(kx^2+ky^2)/sqrt(kx^2+ky^2+kz^2), kz/sqrt(kx^2+ky^2+kz^2)];          
kve = knorm*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];      
kx = Rx * kve; ky = Ry * kve; kz = Rz * kve;              
%      
elseif(r < L(3)*tau) % select acoustic scattering 
Ekins = Ekinf;
knorm = energy2k(Ekins,iv);
% scattering angles theta and phi
r_theta = rand;
r_phi = rand;
phi = 2*pi*r_phi;
theta = acos(1-2*r_theta);
kve = knorm*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];      
kx = kve(1); ky = kve(2); kz = kve(3);

elseif(r < L(4)*tau) % select impurity scattering 
Ekins = Ekinf;
knorm = energy2k(Ekins,iv);
% scattering angles
r_theta = rand;
r_phi = rand;
phi = 2*pi*r_phi;
theta = acos(1-(2*r_theta)/(1+(1-r_theta)*((2*knorm)/(q0))^2));
% rotate back to laboratory frame (anisotropic scattering only)
Rx = [ ky/sqrt(kx^2+ky^2), kx*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)),kx/sqrt(kx^2+ky^2+kz^2)];
Ry = [-kx/sqrt(kx^2+ky^2),ky*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)),ky/sqrt(kx^2+ky^2+kz^2)];
Rz = [ 0,-sqrt(kx^2+ky^2)/sqrt(kx^2+ky^2+kz^2),kz/sqrt(kx^2+ky^2+kz^2)];          
kve = knorm*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];      
kx = Rx * kve; ky = Ry * kve; kz = Rz * kve; 

elseif(r < L(5)*tau) % select intervalley scattering, emission 
    iv = 1.*(iv==2) + 2.*(iv==1); % update valley index

    if iv == 2
Ekins = Ekinf-int_phon_en-Delta;
    elseif iv == 1
         Ekins = Ekinf-int_phon_en+Delta;
    end
   
knorm = energy2k(Ekins,iv);
r_theta = rand;
r_phi = rand;
phi = 2*pi*r_phi;
theta = acos(1-2*r_theta);
kve = knorm*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];      
kx = kve(1); ky = kve(2); kz = kve(3);


elseif(r < L(6)*tau) % select intervalley scattering, absorption 
    iv = 1.*(iv==2) + 2.*(iv==1); % update valley index

    if iv == 2  
Ekins = Ekinf+int_phon_en-Delta;
    elseif iv == 1
         Ekins = Ekinf+int_phon_en+Delta;
    end
knorm = energy2k(Ekins,iv);
r_theta = rand;
r_phi = rand;
phi = 2*pi*r_phi;
theta = acos(1-2*r_theta);
kve = knorm*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];      
kx = kve(1); ky = kve(2); kz = kve(3);
else; end % no scattering has occurred (do nothing)

% count the occupation of the Gamma and L valley
if iv == 1
    gamma_count(j) = gamma_count(j)+1;
elseif iv == 2
        l_count(j) = l_count(j)+1;
end
end  % Monte Carlo main loop -----------------------------------------
%
%
end % end of the loop for different electric fields-------------------
% post-processing: drift velocity, diffusivity, et alia --------------
toc
%% Maxwellian fit
N_bins = 100;
Te = 300; % K
for i = 1:N_field
mesh = linspace(min(vEtoti(i,:)),max(vEtoti(i,:)),N_bins); % select the energy mesh
N = hist(vEtoti(i,:),mesh);
f_m = sqrt(mesh).*exp(-(mesh)./(kB*Te)); % Maxwell distribution
norm_fm = trapz(mesh,f_m);
f_m = f_m./norm_fm;
norm_N = trapz(mesh,N);
N = N./norm_N;
figure 
grid on
hold on
plot(mesh/Q,f_m*Q,'--r','LineWidth',2);
plot(mesh/Q,N*Q,'b','LineWidth',2);
end
ylabel('Electron distribution');
xlabel('Energy, eV');
legend('Maxwellian fit, T_e = 300 K','MC');
title(['F =',num2str(Fx(i)/(10^5)),'kV/cm, T = 300 K']);
hold off
grid off
%% Plot average energy
for i = 1:N_field
vEkin_tau(i,:) = (vEkini(i,:) + vEkinf(i,:))/2; % average energy during drift, J
E_av(i,:) = ((cumsum(vEkin_tau(i,:)))*tau)./time; % average energy, J
E_av(i,:) = (E_av(i,:)*10^3)/Q; % average energy, meV
reference = ones(1,nt)*1.5*kB*T; % reference for convergence, J
reference = reference*10^3/Q;  % reference for convergence, meV
figure(2)
hold on
grid on
plot(time/10^-9,E_av(i,:),'Linewidth',2);
end
plot(time/10^-9,reference,'--r','LineWidth',2);
grid off
hold off
xlabel('Time, ns');
ylabel('Average energy, meV');
legend('Fx = 0.1 kV/cm','Fx = 1 kV/cm','Fx = 3 kV/cm','Fx = 4 kV/cm','Fx = 5 kV/cm','Fx = 7 kV/cm','Fx = 9 kV/cm','3/2 kBT');
ylim([0 150]);


%% Plot average velocity
for i = 1:N_field
vDeltaEkin(i,:) =  vEkinf(i,:) - vEkini(i,:); % energy variation during drift, J 
Vel_av(i,:) = (cumsum(vDeltaEkin(i,:))./(time*Q*Fx(i)))*10^-5; % average velocity, 10^7 cm/s
figure(3)
hold on
grid on
plot(time/10^-9,Vel_av(i,:),'LineWidth',2);
xlabel('Time, ns');
ylabel('Average velocity, 10^7 cm/s');
end
grid off
hold off
legend('Fx = 0.1 kV/cm','Fx = 1 kV/cm','Fx = 3 kV/cm','Fx = 4 kV/cm','Fx = 5 kV/cm','Fx = 7 kV/cm','Fx = 9 kV/cm');
ylim([-2 2]);

%% Plot velocity vs Field
figure(4)
plot(Fx*10^-5,Vel_av(:,nt),'-o');
xlabel('Electric field, kV/cm');
ylabel('Velocity, 10^7 cm/s');

%% Plot average occupation of Gamma and L valleys
l_count = l_count/nt;
gamma_count = gamma_count/nt;
figure(5)
hold on
grid on
plot(Fx*10^-5,l_count,'b');
plot(Fx*10^-5,gamma_count,'r');
grid off
hold off
ylabel('Electron occupancy');
xlabel('Electric field, kV/cm');
legend('L-valley','\Gamma-valley');

%% Plot autocorrelation function
Dif_coeff = zeros(1,N_field);
time_step = 1e-15; % s
for i = 1:N_field
MAXLAG = 2000;
deltav = vel(i,1:10:end) - mean(vel(i,:));
[Cv,lags] = xcorr(deltav,deltav,MAXLAG);
Dif_coeff(i) = trapz(lags(MAXLAG:end)*time_step,Cv(MAXLAG:end)/length(deltav));
figure(6)
hold on
grid on
plot(lags*time_step*1e12,Cv/length(deltav))
end
set(gca,'FontSize',14,'FontName','Arial','box','on')
xlabel('time, ps')
ylabel('Autocorrelation function')
grid off
hold off
xlim([0 2])
legend('Fx = 0.1 kV/cm','Fx = 1 kV/cm','Fx = 3 kV/cm','Fx = 4 kV/cm','Fx = 5 kV/cm','Fx = 7 kV/cm','Fx = 9 kV/cm');

%% Plot diffusion coefficient as a function of the electric field
figure(7)
hold on
grid on
plot(Fx*10^-5,Dif_coeff*10^2);
ylabel('Diffusion coefficient (10^2 cm^2/s)');
xlabel('Electric field, kV/cm');
grid off
hold off

%% Compare with the relaxation time approxiation
mobility;
D_ein = mu_0_tot(296)*kB*300/Q;
disp(D_ein*10^2);
