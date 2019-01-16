% neff_V.m: effective index as a function of voltage for carrier-depletion phase modulation 
% Wei Shi, UBC, 2012
% Usage, e.g.:
%    [del_neff alpha Rj Cj]=neff_V(1.55e-6, 220e-9, 500e-9, 90e-9, 3.47, 1.44, 1.44, 50e-9, 1e-6, 1e-6, 500, 25, -1)
function [neff alpha Rj Cj]=neff_V(lambda, t, w, t_slab, n_core, n_clad, n_oxide, pn_offset, ds_n_plus, ds_p_plus, pts, T, V)
[n, p, xdoping, xn, xp, Rj, Cj]=pn_depletion(w, pn_offset, ds_n_plus, ds_p_plus, T, V, pts);

M=min(ds_n_plus-pn_offset+0.5e-6, ds_p_plus+pn_offset+0.5e-6)/w-0.5;

[xwg, TM_E_TEwg, neff0]=wg_TElike_1Dprofile_neff(lambda, t, w, t_slab, n_core, n_clad, n_oxide, pts, M);
Ewg=TM_E_TEwg(:,1)';

pts_x=length(xwg);
dxwg=zeros(1, pts_x); 
dxwg(1)=xwg(2)-xwg(1); dxwg(pts_x)=xwg(pts_x)-xwg(pts_x-1);
for i=2:pts_x-1
    dxwg(i)=xwg(i+1)/2-xwg(i-1)/2;
end
    
n_wg=interp1(xdoping, n, xwg);
p_wg=interp1(xdoping, p, xwg);

del_ne=-3.64e-10*lambda^2*sum(conj(Ewg).*(n_wg*1e-6).*Ewg.*dxwg)/sum(conj(Ewg).*Ewg.*dxwg);
del_nh=-3.51e-6*lambda^2*sum(conj(Ewg).*(p_wg*1e-6).^0.8.*Ewg.*dxwg)/sum(conj(Ewg).*Ewg.*dxwg);
del_neff=del_ne+del_nh;
neff=neff0+del_neff;

del_alpha_e=3.52e-6*lambda^2*sum(conj(Ewg).*(n_wg*1e-6).*Ewg.*dxwg)/sum(conj(Ewg).*Ewg.*dxwg);
del_alpha_h=2.4e-6*lambda^2*sum(conj(Ewg).*(p_wg*1e-6).*Ewg.*dxwg)/sum(conj(Ewg).*Ewg.*dxwg);
alpha=del_alpha_e+del_alpha_h;
