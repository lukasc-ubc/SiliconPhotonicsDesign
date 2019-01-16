% example:
% [neff alpha delta_neff delta_phi fc]=neff_V_plot(1.5e-6, 220e-9, 500e-9, 90e-9, 3.47, 1.44, 1.44, 50e-9, 10e-6, 10e-6, 500, 25, -(1:5));

function [neff alpha delta_neff delta_phi fc] = neff_V_plot(lambda, t, w, t_slab, n_core, n_clad, n_oxide, pn_offset, ds_n_plus, ds_p_plus, pts, T, V);

neff=zeros(1, length(V)); alpha=zeros(1, length(V));
Rj=zeros(1, length(V)); Cj=zeros(1, length(V));
for i=1:length(V);
    [neff(i) alpha(i) Rj(i) Cj(i)]=neff_V(lambda, t, w, t_slab, n_core, n_clad, n_oxide, pn_offset, ds_n_plus, ds_p_plus, pts, T, V(i))
end

[neff_v0 alpha_v0]=neff_V(lambda, t, w, t_slab, n_core, n_clad, n_oxide, pn_offset, ds_n_plus, ds_p_plus, pts, T, 0);

delta_neff=neff-neff_v0;
alpha_dB=-10*log10(exp(-alpha));

figure; plot(-V, delta_neff)
figure; plot(-V, alpha_dB);

% Phase shift per cm
delta_phi=2*pi/lambda*delta_neff*1e-2/pi;% per cm
figure; plot(-V, delta_phi, 'linewidth', 2);

% Cut-off frequency
fc=1./(2*pi*Rj.*Cj)*1e-9;% in GHz
figure; plot(-V, fc, 'linewidth', 2);
