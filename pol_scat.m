function [Wpop_emi,Wpop_abs] = pol_scat(E,T,i,alpha) % insert alpha in input to run scattering.m
%%%%%%%%    Constants    %%%%%%%%%%%%%
c_light = 2.99792458e+8; % light velocity, m/s
H = 6.626070040e-34; % Planck constant, J*s
HBAR = 1.054e-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
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
meff = masses(iv); % effective mass, k

%alpha = 0; %%%% Uncomment this line for the MC simulation

hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
Daco = 7*Q; % acoustic deformation potential, J
egap = 1.424*Q; % energy gap, J

eps_p = 1/((1/eps_infty)-(1/eps_s)); %F/m
gamma =@(x) x*(1+alpha*x); % energy, J
E_i = E;
Wpop = (Q^2*sqrt(meff)*(hwpop/HBAR))/(4*pi*eps_p*sqrt(2)*HBAR*sqrt(gamma(E_i)));
dgamma_dE =@(x) 1+2*alpha*x;
Nq = 1/(exp(hwpop/(kB*T))-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alpha == 0.64/Q
%% Compute the polar optical scattering rate (non parabolic approximation)
%% Emission case
E_f = E_i-hwpop;

A = (2*(1+alpha*E_i)*(1+alpha*E_f)+alpha*(gamma(E_i)+gamma(E_f)))^2;
B = -2*alpha*sqrt(gamma(E_i))*sqrt(gamma(E_f))*(4*(1+alpha*E_i)*(1+alpha*E_f)+alpha*(gamma(E_i)+gamma(E_f)));
C = 4*(1+2*alpha*E_i)*(1+2*alpha*E_f)*(1+alpha*E_i)*(1+alpha*E_f);

if E_i<hwpop
    Wpop_emi = 0;
else
Wpop_emi = Wpop*dgamma_dE(E_f)*(Nq+1)*(1/C)*(A*log(abs(((gamma(E_i))^0.5+(gamma(E_f))^0.5)/((gamma(E_i))^0.5-(gamma(E_f))^0.5)))+B);
end

%% Absorption case
E_f = E_i+hwpop;
A = (2*(1+alpha*E_i)*(1+alpha*E_f)+alpha*(gamma(E_i)+gamma(E_f)))^2;
B = -2*alpha*sqrt(gamma(E_i))*sqrt(gamma(E_f))*(4*(1+alpha*E_i)*(1+alpha*E_f)+alpha*(gamma(E_i)+gamma(E_f)));
C = 4*(1+2*alpha*E_i)*(1+2*alpha*E_f)*(1+alpha*E_i)*(1+alpha*E_f);

Wpop_abs = Wpop*dgamma_dE(E_f)*Nq*(1/C)*(A*log(abs(((gamma(E_i))^0.5+(gamma(E_f))^0.5)/((gamma(E_i))^0.5-(gamma(E_f))^0.5)))+B);

end

if alpha == 0
%% Compute the polar optical scattering rate (parabolic approximation)
k = sqrt(2*meff*E)*(1/HBAR);
Wpop = ((Q^2*(hwpop/HBAR))/(8*pi*eps_p))*(k/E);

%% Emission case
if E < hwpop
    Wpop_emi = 0;
else
qmin = k*abs(1-sqrt(1-hwpop/E));
qmax = k*abs(1+sqrt(1-hwpop/E));
Wpop_emi = Wpop*(Nq+1)*log(qmax/qmin);
end
%% Absorption case
qmin = k*abs(1-sqrt(1+hwpop/E));
qmax = k*abs(1+sqrt(1+hwpop/E));
Wpop_abs = Wpop*Nq*log(qmax/qmin);
end
