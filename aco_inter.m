function [Waco_inter_em,Waco_inter_abs] = aco_inter(E,T,iv)
%%%%%%%%    Constants    %%%%%%%%%%%%%
c_light = 2.99792458e+8; % light velocity, m/s
H = 6.626070040e-34; % Planck constant, J*s
HBAR = 1.054e-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
%T = 300; % temperature, K
eps0 = 8.854e-12; % vacuum permittivity constant, F/m
M0 = 9.1095e-31; % electron mass, kg
meff_G = 0.067*M0; %kg
ml = 1.2*M0; % longitudinal mass in the L valley, kg
mt = 0.2*M0; % transversal mass in the L valley, kg
meff_L = (ml*mt^2)^(1/3); % kg
masses = [meff_L meff_G];
meff = masses(iv); % effective mass, kg
Zj = [4 1]; %number of equivalent final valleys (4 for L and 1 for Gamma)
%
rho = 5320; % mass density, kg/m^3
eps_s = 12.9*eps0; % static dielectric constant, F/m
eps_infty = 10.9*eps0; % high-frequency dielectric constant, F/m
v_l = 5240; % longitudinal sound velocity, m/s

%alpha = 0.64/Q; % nonparabolicity factor, 1/J
alpha = 0;
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
Daco = 7*Q; % acoustic deformation potential, J
egap = 1.424*Q; % energy gap, J

Dij = 1e11*Q; % Intervalley coupling constant, J/m
int_phon_en = 28e-3*Q; % intervalley phonon energy, J
Delta = 0.3*Q; % difference of energy minima, J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if iv == 1
Ef_em = E-int_phon_en-Delta;
Ef_ab = E+int_phon_en-Delta;
elseif iv == 2
    Ef_em = E-int_phon_en+Delta;
    Ef_ab = E+int_phon_en+Delta;
end

Waco_inter = (pi*Dij^2*Zj(iv))/(rho*(int_phon_en/HBAR));
Nij = 1/(exp(int_phon_en/(kB*T))-1);
DOS_em = ((2*meff)^1.5*sqrt(Ef_em))/(4*pi^2*HBAR^3);
DOS_ab = ((2*meff)^1.5*sqrt(Ef_ab))/(4*pi^2*HBAR^3);

DOS_em = real(DOS_em);
DOS_ab = real(DOS_ab);

Waco_inter_em = Waco_inter*(Nij+1)*DOS_em;
Waco_inter_abs = Waco_inter*Nij*DOS_ab;