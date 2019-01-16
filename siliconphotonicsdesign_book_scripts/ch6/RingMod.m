% RingMod.m: Ring modulator 1D model
% Usage, e.g.,
%   [Ethru Edrop Qi Qc Rj Cj]=RingMod(1.55e-6, 'all-pass', 10e-6, 0, 2*pi*10e-6, 500e-9, 0, 1e-6, 1e-6, 25, 0);
% 
% Wei Shi, UBC, 2012
% weis@ece.ubc.ca
%
function [Ethru Edrop Qi Qc tau_rt Rj Cj]=RingMod(lambda, Filter_type, r, Lc, L_pn, w, pn_offset, ds_n_plus, ds_p_plus, T, V);
%
% type: all-pass or add-drop
% r: radius
% Lc: coupler length
% Lh: heater length
%
% neff_pn, alpha_pn: effective index and free-carrier obsorption of the phase modulator
% Rj, Cj: junction resistance and capacitance of the phase modulator
%
% predetermined parameters
t=220e-9; t_slab=90e-9; n_core=3.47; n_clad=1.44; n_oxide=1.44; pts=200;
%
[neff_pn alpha_pn Rj Cj]=neff_V(lambda, t, w, t_slab, n_core, n_clad, n_oxide, pn_offset, ds_n_plus, ds_p_plus, pts, T, V);
%
% undoped waveguide mode and effective index
[xwg0 TM_E_TEwg0 neff0]=wg_TElike_1Dprofile_neff(lambda, t, w, t_slab, n_core, n_clad, n_oxide, pts, 2);
neff_exc=neff0;
del_lambda=0.1e-9;
[xwg1 TM_E_TEwg1 neff0_1]=wg_TElike_1Dprofile_neff(lambda+del_lambda, t, w, t_slab, n_core, n_clad, n_oxide, pts, 2);
ng=neff0-(neff0_1-neff0)/del_lambda*lambda;

alpha_wg_dB=5; % optical loss of intrinsic optical waveguide, in dB/cm
alpha_wg=-log(10^(-alpha_wg_dB/10));% converted to /cm
alpha_pn=alpha_wg+alpha_pn;
alpha_exc=alpha_wg; % optical loss of the ring cavity excluding the phase modulator

L_rt=Lc*2+2*pi*r;
L_exc=L_rt-L_pn;
phi_pn=(2*pi/lambda)*neff_pn*L_pn;
phi_exc=(2*pi/lambda)*neff_exc*L_exc;
phi_rt=phi_pn+phi_exc;

c=299792458;
vg=c/ng;
tau_rt=L_rt/vg;% round-trip time

A_pn=exp(-alpha_pn*100*L_pn); % attunation due to pn junciton
A_exc=exp(-alpha_exc*100*L_exc); % attunation over L_exc
A=A_pn*A_exc; % round-trip optical power attenuation

alpha_av=-log(A)/L_rt;% average loss of the cavity
Qi=2*pi*ng/lambda/alpha_av;

%coupling coefficients
k=0.2;
if (Filter_type=='all-pass')
    t=sqrt(1-k^2);
    Ethru=(-sqrt(A)+t*exp(-1i*phi_rt))/(-sqrt(A)*conj(t)+exp(-1i*phi_rt));
    Edrop=0;
    Qc=-(pi*L_rt*ng)/(lambda*log(abs(t)));
elseif (Filter_type=='add-drop')
    k1=k; k2=k1;
    t1=sqrt(1-k1^2);  t2=sqrt(1-k2^2);
    Ethru=(t1-conj(t2)*sqrt(A)*exp(1i*phi_rt))/(1-sqrt(A)*conj(t1)*conj(t2)*exp(1i*phi_rt));
    Edrop=-conj(k1)*k2*sqrt(sqrt(A))*exp(1i*phi_rt/2)/(1-sqrt(A)*conj(t1)*conj(t2)*exp(1i*phi_rt));
    Qc1=-(pi*L_rt*ng)/(lambda*log(abs(t1)));
    Qc2=-(pi*L_rt*ng)/(lambda*log(abs(t2)));
    Qc=1/(1/Qc1+1/Qc2);
else
    error(1, 'The''Filter_type'' has to be ''all-pass'' or ''add-drop''.\n');
end