% pn_depletion.m: 1D pn junction model for carrier-depletion phase modulation
% Wei Shi, UBC, Nov. 2012
%
% usage, e.g.:
%   [n, p, x, xn, xp, Rj, Cj]=pn_depletion(500e-9, 50e-9, 1e-6, 1e-6, 25, -1, 100)
%
function [n, p, x, xn, xp, Rj, Cj]=pn_depletion(wg_width, pn_offset, ds_npp, ds_ppp, T, V, pts)
%
% N_D, N_A: doping densities
% V: applied voltage; positive for forward bias; negative for reverse bias
% ds_npp: distance of the n++ boundary to the pn junction centre
% ds_ppp: distance of the p++ boundary to the pn junction centre
% Rj: junction resistance in ohms
% Cj: junction capacitance in F/m

epsilon0 = 8.854187817620e-12; % [F/m] 
epsilon_s = 11.8; % relative dielectric constant for Si 
q = 1.60217646e-19; % electronic charge [Coulumbs]
kB = 1.3806503e-23; % Boltzmann constant in J/K
T=T+273.15; % Temperature [K]
VT=kB*T/q;

%material constants
NA_plus=4.4e20*1e6;% cm^-3*1e6
ND_plus=4.4e20*1e6;
NA=5e17*1e6;% cm^-3*1e6
ND=3e17*1e6;

Rs_rib_n=2.5e3;
Rs_rib_p=4.0e3;
Rs_slab_n=0.6e4;
Rs_slab_p=1e4;

% waveguide height
h_rib=220e-9;   h_slab=90e-9;

h=4.135e-15; % Plank's constant [eV-s]
m_0=9.11e-31; % electron mass [kg]
m_n=1.08*m_0; % Density-of-states effective mass for electrons
m_p=1.15*m_0; % Density-of-states effective mass for holes
Nc=2*(2*pi*m_n*(kB/q)*T/h^2)^(3/2)/(q)^(3/2); % Effective Density of states for Conduction Band
Nv=2*(2*pi*m_p*(kB/q)*T/h^2)^(3/2)/(q)^(3/2); % Effective Density of states for Valence Band
Eg=1.1242; % band gap for Si [eV]
% ni=1e10*1e6;
ni=sqrt(Nc*Nv).*exp(-Eg/(2*(kB/q)*T)); % intrinsict charge carriers in m^-3

Vbi=VT*log(NA*ND/ni^2); % built-in or diffusion potential
Wd=sqrt(2*epsilon0*epsilon_s*(NA+ND) / (q*NA*ND) *(Vbi-V)); % depletion width
xp=-Wd/(1+NA/ND)+pn_offset;
xn=Wd/(1+ND/NA)+pn_offset;

del_x=wg_width/(pts-1);
x_ppp=-ds_ppp+pn_offset;  x_npp=ds_npp+pn_offset;
x_min=x_ppp-500e-9;  x_max=x_npp+500e-9;
%
x_NA_plus=x_min:del_x:x_ppp-del_x;
x_NA=x_ppp:del_x:xp-del_x;
x_dep=xp:del_x:xn;% for the depletion region
x_ND=xn+del_x:del_x:x_npp;
x_ND_plus=x_npp+del_x:del_x:x_max;
x=[x_NA_plus, x_NA, x_dep, x_ND, x_ND_plus];

n0_NA=ni^2/NA;  p0_ND=ni^2/ND;
n0_NA_plus=ni^2/NA_plus;  p0_ND_plus=ni^2/ND_plus;

% Long-base assumption
% Lp=sqrt(Dp*tau_p);
% Ln=sqrt(Dn*tau_n);
% del_n_NA=n0_NA*(exp(q*V/(kB*T))-1)* exp(-abs(x_NA-xp)/Ln); % minority electron density in p(NA) region
% del_p_ND=p0_ND*(exp(q*V/(kB*T))-1)* exp(-abs(x_ND-xn)/Lp); % minority hole density in n(ND) region

% Short-base assumption
del_n_NA=n0_NA*(exp(q*V/(kB*T))-1)* (1-abs((x_NA-xp)/(xp-x_ppp))); % minority electron density in p(NA) region
del_p_ND=p0_ND*(exp(q*V/(kB*T))-1)* (1-abs((x_ND-xn)/(x_npp-xn))); % minority hole density in n(ND) region


n_NA=n0_NA+del_n_NA;  p_ND=p0_ND+del_p_ND;
p_dep=zeros(1, length(x_dep));  n_dep=zeros(1, length(x_dep));
p_NA=ones(1, length(x_NA))*NA; % majority holes in p(NA) region
n_ND=ones(1, length(x_ND))*ND; % majority electrons in n(ND) region
n_NA_plus=ones(1, length(x_NA_plus))*n0_NA_plus;% assumption of uniform electrons in p++  
p_ND_plus=ones(1, length(x_ND_plus))*p0_ND_plus;% assumption of uniform holes in n++
p_NA_plus=ones(1, length(x_NA_plus))*NA_plus; % majority holes in p++ region
n_ND_plus=ones(1, length(x_ND_plus))*ND_plus; % majority electrons in n++ region

n=[n_NA_plus, n_NA, n_dep, n_ND, n_ND_plus]; p=[p_NA_plus, p_NA, p_dep, p_ND, p_ND_plus];

Rj=(wg_width/2-xn)* Rs_rib_n+(wg_width/2+xp)* Rs_rib_p+(-wg_width/2-x_ppp)* Rs_slab_p+(x_npp-wg_width/2)* Rs_slab_n;
Cj=sqrt(q*epsilon0*epsilon_s/2/ (1/ND+1/NA)/(Vbi-V))*h_rib;