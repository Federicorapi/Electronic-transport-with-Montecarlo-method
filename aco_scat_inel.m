function [Waco_emi,Waco_abs] = aco_scat_inel(E,T)
%% Compute the emission and absorption scattering rates depending on E and T

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
meff = 0.067*M0; % effective mass, kg

alpha = 0.64/Q; % nonparabolicity factor, 1/J
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
Daco = 7*Q; % acoustic deformation potential, J
egap = 1.424*Q; % energy gap, J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Gamma = E*(1+alpha*E);
Es = (meff*v_l^2)/2;
C = (4*sqrt(Es))/(kB*T*(1-4*alpha*Es));
control = Es/(1-4*alpha*Es);

F1 =@(x) (x^12/574801920 - x^10/12096000 + x^8/241920 - x^6/4320 + x^4/48 - x^3/6 + x^2/2);
F2 =@(x) ((x^3*(5*x^10 - 234*x^8 + 11440*x^6 - 617760*x^4 + 51891840*x^2 - 389188800*x + 1037836800))/3113510400);
G1 =@(x) (x^12/574801920 - x^10/12096000 + x^8/241920 - x^6/4320 + x^4/48 + x^3/6 + x^2/2);
G2 =@(x) ((x^3*(5*x^10 - 234*x^8 + 11440*x^6 - 617760*x^4 + 51891840*x^2 + 389188800*x + 1037836800))/3113510400);

Waco = (meff^(0.5)*(kB*T)^3*Daco^2) / ((2^(2.5))*pi*HBAR^4*v_l^4*rho*sqrt(Gamma));

%%Integration limits
if Gamma < control
    
    x1_a = C*(sqrt(Es)*(1+2*alpha*E)-sqrt(Gamma));
    x2_a = C*(sqrt(Es)*(1+2*alpha*E)+sqrt(Gamma));
   
   Waco_emi = 0;
   Waco_abs = Waco*((1+2*alpha*E)*(F1(x2_a)-F1(x1_a))+2*alpha*kB*T*(F2(x2_a)-F2(x1_a)));
end
 if Gamma > control
     
        x1_a = 0;
        x2_a = C*(sqrt(Gamma)+sqrt(Es)*(1+2*alpha*E));
        x1_e = 0;
        x2_e = C*(sqrt(Gamma)-sqrt(Es)*(1+2*alpha*E));
       
        Waco_emi = Waco*((1+2*alpha*E)*(G1(x2_e)-G1(x1_e))-2*alpha*kB*T*(G2(x2_e)-G2(x1_e)));
        Waco_abs = Waco*((1+2*alpha*E)*(F1(x2_a)-F1(x1_a))+2*alpha*kB*T*(F2(x2_a)-F2(x1_a)));
       
  end
end




