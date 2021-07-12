function [Wimp] = imp_scat(E,nI,T,iv)
%%%%%%%%    Constants    %%%%%%%%%%%%%
c_light = 2.99792458e+8; % light velocity, m/s
H = 6.626070040e-34; % Planck constant, J*s
HBAR = 1.054e-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
%T = 300; % temperature, K
eps0 = 8.854e-12; % vacuum permittivity constant, F/m
M0 = 9.1095e-31; % electron mass, kg
%
rho = 5320; % mass density, kg/m^3
eps_s = 12.9*eps0; % static dielectric constant, F/m
eps_infty = 10.9*eps0; % high-frequency dielectric constant, F/m
v_l = 5240; % longitudinal sound velocity, m/s
meff_G = 0.067*M0; %kg
ml = 1.2*M0;
mt = 0.2*M0;
meff_L = (ml*mt^2)^(1/3);
masses = [meff_G meff_L];
meff = masses(iv); % effective mass, kg
%alpha = 0.64/Q; % nonparabolicity factor, 1/J
alpha = 0;
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
Daco = 7*Q; % acoustic deformation potential, J
egap = 1.424*Q; % energy gap, J 
gamma =@(x) x*(1+alpha*x); %energy, J
dgamma_dE =@(x) 1+2*alpha*x;
DOS_par = ((2*meff)^1.5*sqrt(E))/(4*pi^2*HBAR^3);
DOS_nonpar = DOS_par*dgamma_dE(E);
q0 = sqrt((Q^2*nI)/(eps_s*kB*T));  % reciprocal screening length
k2 = (2*meff*gamma(E))/HBAR^2;

Wimp = (pi*nI*1*Q^4)/(HBAR*eps_s^2)*DOS_nonpar*(2/(q0^2*(4*k2+q0^2)));

