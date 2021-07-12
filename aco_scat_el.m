function [Waco_el_par,Waco_el_nonpar] = aco_scat_el(E,T,iv)
%% Compute the acoustic scattering rate in the elastic approximation 

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
ml = 1.2*M0; % longitudinal mass in the L valley, kg
mt = 0.2*M0; % transversal mass in the L valley, kg
meff_L = (ml*mt^2)^(1/3);
masses = [meff_G meff_L];
meff = masses(iv); % effective mass, kg

alpha = 0.64/Q; % nonparabolicity factor, 1/J
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
Daco = 7*Q; % acoustic deformation potential, J
egap = 1.424*Q; % energy gap, J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dGamma_dE = 1+2*alpha*E;

Waco = (2*pi*Daco^2*kB*T)/(HBAR*v_l^2*rho);

DOS_par = ((2*meff)^1.5*sqrt(E))/(4*pi^2*HBAR^3);
DOS_nonpar = DOS_par*dGamma_dE;

Waco_el_par = Waco*DOS_par;
Waco_el_nonpar = Waco*DOS_nonpar;