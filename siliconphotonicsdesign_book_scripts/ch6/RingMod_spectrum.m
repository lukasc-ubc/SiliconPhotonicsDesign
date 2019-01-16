% calculate the ring modulator spectrum
% Wei Shi UBC, 2012 
% weis@ece.ucb.ca

function [Ethru Edrop Qi Qc tau_rt Rj Cj]=RingMod_spectrum(lambda, Filter_type, r, Lc, L_pn, w, pn_offset, ds_n_plus, ds_p_plus, T, V);
%
Ethru=zeros(1, length(lambda));
Edrop=zeros(1, length(lambda));
Qi=zeros(1, length(lambda));
Qc=zeros(1, length(lambda));
tau_rt=zeros(1, length(lambda));
%
for i=1:length(lambda)
    [Ethru(i) Edrop(i) Qi(i) Qc(i) tau_rt(i) Rj Cj]=RingMod(lambda(i), Filter_type, r, Lc, L_pn, w, pn_offset, ds_n_plus, ds_p_plus, T, V);
end