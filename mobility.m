%% Plot Waco_emi and Waco_abs
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
eps_p = 1/((1/eps_infty)-(1/eps_s)); %F/m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iv = 1;


nT = 1000;  % number of temperature points
T = linspace(7,1000,nT); % temperature axis, K

nE = 1000; % number of energy points
vE = linspace(0.0001,1,nE)*Q; % energy axis, J

mu_0_tot = zeros(1,nT);
mu_0_tot_inv = zeros(1,nT);
mu_0_aco = zeros(1,nT); % acoustic scattering mobility, m^2/(Vs)
mu_0_imp = zeros(1,nT); % impurity scattering mobility, m^2/(Vs)
mu_0_pop = zeros(1,nT); % polar optical scattering mobility, m^2/(Vs)

tau_pop = zeros(1,nE); % relaxation time, s
tau_aco = zeros(1,nE); % relaxation time, s
tau_imp = zeros(1,nE); % relaxation time, s

fM = zeros(1,nE);



for i = 1:nT    
for ie = 1:nE, E = vE(ie); 
    %% Acoustic scattering tau
    [Waco_par, Waco_nonpar] = aco_scat_el(E,T(i),iv);
    tau_aco(ie) = 1/Waco_par;
    
%% Impurity scattering tau  
    nI = 10^20;
    q0 = sqrt((Q^2*nI)/(eps_s*kB*T(i)));
    k2 = (2*meff*E)/(HBAR^2);
    DOS_par = ((2*meff)^1.5*sqrt(E))/(4*pi^2*HBAR^3);
    tau_imp_inv = ((pi*nI*1*Q^4)/(HBAR*eps_s^2))*((1)/(4*k2^2))*(log(1+(4*k2)/(q0^2))-1)*DOS_par;
    tau_imp(ie) = 1/tau_imp_inv;
    
 %% Polar optical scattering tau   
    Nq = 1/(exp(hwpop/(kB*T(i)))-1);
    
    pop_abs = Nq*(sqrt((E+hwpop)/E)-(hwpop/E)*asinh(sqrt(E/hwpop)));
    
    if E<hwpop
        pop_emi = 0; 
    else
        pop_emi = (Nq+1)*(sqrt((E-hwpop)/E)+(hwpop/E)*asinh(sqrt((E-hwpop)/hwpop)));
    end
    
    tau_pop_inv = (Q^2*(hwpop/HBAR))/(8*pi*eps_p)*(sqrt(k2))/(E)*(pop_emi + pop_abs);
    tau_pop(ie) = 1/tau_pop_inv;
    
 %% Maxwell ditribution
 beta = 1/(kB*T(i));
 fM(ie) = exp(-beta*E);
    
end

mu_0_aco(i) = (Q/meff)*2/(3*kB*T(i))*(trapz(vE,fM.*vE.*tau_aco.*DOS_par))/(trapz(vE,fM.*DOS_par));
mu_0_imp(i) = (Q/meff)*2/(3*kB*T(i))*(trapz(vE,fM.*vE.*tau_imp.*DOS_par))/(trapz(vE,fM.*DOS_par));
mu_0_pop(i) = (Q/meff)*2/(3*kB*T(i))*(trapz(vE,fM.*vE.*tau_pop.*DOS_par))/(trapz(vE,fM.*DOS_par));

mu_0_tot_inv(i) = (1/mu_0_aco(i)) + (1/mu_0_imp(i)) + (1/mu_0_pop(i));
mu_0_tot(i) = 1/ mu_0_tot_inv(i);
end


figure(7), hold on
grid on
plot(T,mu_0_aco,'--','LineWidth',2);
plot(T,mu_0_imp,'--','LineWidth',2);
plot(T,mu_0_pop,'--','LineWidth',2);
plot(T,mu_0_tot,'LineWidth',2);
ylim([0 10^5]);
set(gca,'FontSize',14,'FontName','Arial','box','on','XScale','log','YScale','log')
ylabel('Mobility (m^2 V^-1 s^-1)'), xlabel('Temperature (K)')
legend('\mu_{aco}','\mu_{imp}','\mu_{pop}','\mu_{tot}');
% title('Impurity scattering rate')
 hold off


