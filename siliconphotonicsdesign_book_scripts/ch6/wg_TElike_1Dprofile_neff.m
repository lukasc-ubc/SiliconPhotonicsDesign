% wg_TElike_1Dprofile.m - Effective Index Method - 1D mode profile
% Lukas Chrostowski, 2012
% modified by Wei Shi, 2012

% usage, e.g.:
%  [xwg, TM_E_TEwg]=wg_TElike_1Dprofile_neff (1.55e-6, 0.22e-6, 0.5e-6, 90e-9,3.47, 1, 1.44, 100, 2);
%  figure; plot(xwg, TM_E_TEwg(:,1))

function [xwg, TM_E_TEwg, neff_TEwg_1st]=wg_TElike_1Dprofile_neff (lambda, t, w, t_slab, n_core, n_clad, n_oxide, pts, M)

% TE (TM) modes of slab waveguide (core and slab portions):
[nTE,nTM]=wg_1D_analytic (lambda, t, n_oxide, n_core, n_clad);
if t_slab>0
    [nTE_slab,nTM_slab]=wg_1D_analytic (lambda, t_slab, n_oxide, n_core, n_clad);
else
    nTE_slab=n_clad; nTM_slab=n_clad;
end
[xslab, TE_Eslab, TE_Hslab, TM_Eslab, TM_Hslab]= wg_1D_mode_profile (lambda, t, n_oxide, n_core, n_clad, pts, M);

% TE-like modes of the etched waveguide (for fundamental slab mode):
[nTE,nTM]=wg_1D_analytic (lambda, w, nTE_slab(1), nTE(1), nTE_slab(1));
neff_TEwg_1st=nTM(1);
[xwg, TE_E_TEwg, TE_H_TEwg, TM_E_TEwg, TM_H_TEwg]= wg_1D_mode_profile (lambda, w, nTE_slab(1), nTE(1), nTE_slab(1), pts, M);
