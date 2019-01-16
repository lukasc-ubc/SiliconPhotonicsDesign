% wg_1D_analytic.m - Analytic solution of waveguide
% by Lumerical Solutions, http://www.lumerical.com/mode_online_help/slab_wg.m
% modified by Lukas Chrostowski, 2012
% See Yariv Photonics book, Chapter 3
% finds the TE and TM effective indices of a 3-layer waveguide

% usage:
%  - get effective indices for supported modes:
%  [nTE, nTM] = wg_1D_analytic2 (1.55e-6, 0.22e-6, 1.444, 3.47, 1.444)
%  - TEparam,TMparam: h, q, p parameters of the mode.

function [nTE,nTM,TEparam,TMparam]=wg_1D_analytic (lambda, t, n1, n2, n3)
k0 = 2*pi/lambda;
b0 = linspace( max([n1 n3])*k0, n2*k0, 1000);   %k0*n3 < b < k0*n2
b0 = b0(1:end-1);
te0=TE_eq(b0,k0,n1,n2,n3,t);
tm0=TM_eq(b0,k0,n1,n2,n3,t);

%TE
intervals=(te0>=0)-(te0<0);
izeros=find(diff(intervals)<0);
X0=[b0(izeros); b0(izeros+1)]';
[nzeros,scrap]=size(X0);
for i=1:nzeros
    nTE(i)=fzero(@(x) TE_eq(x,k0,n1,n2,n3,t),X0(i,:))/k0;
    [TEparam(i,1),TEparam(i,2),TEparam(i,3),TEparam(i,4)]= TE_eq(nTE(i)*k0,k0,n1,n2,n3,t);
end
nTE=nTE(end:-1:1);
TEparam=TEparam(end:-1:1,:);

%TM
intervals=(tm0>=0)-(tm0<0);
izeros=find(diff(intervals)<0);
X0=[b0(izeros); b0(izeros+1)]';
[nzeros,scrap]=size(X0);
for i=1:nzeros
    nTM(i)=fzero(@(x) TM_eq(x,k0,n1,n2,n3,t),X0(i,:))/k0;
    [TMparam(i,1),TMparam(i,2),TMparam(i,3),TMparam(i,4)]= TM_eq(nTM(i)*k0,k0,n1,n2,n3,t);
end
if nzeros>0
    nTM=nTM(end:-1:1);
    TMparam=TMparam(end:-1:1,:);
else
    nTM=[];
end

function [te0,h0,q0,p0]=TE_eq(b0,k0,n1,n2,n3,t)
h0 = sqrt( (n2*k0)^2 - b0.^2 );
q0 = sqrt( b0.^2 - (n1*k0)^2 );
p0 = sqrt( b0.^2 - (n3*k0)^2 );
%the objective is to find zeroes of te0 and tm0
te0 = tan( h0*t ) - (p0+q0)./h0./(1-p0.*q0./h0.^2);

function [tm0,h0,q0,p0]=TM_eq(b0,k0,n1,n2,n3,t)
h0 = sqrt( (n2*k0)^2 - b0.^2 );
q0 = sqrt( b0.^2 - (n1*k0)^2 );
p0 = sqrt( b0.^2 - (n3*k0)^2 );
pbar0 = (n2/n3)^2*p0;
qbar0 = (n2/n1)^2*q0;
tm0 = tan( h0*t ) - h0.*(pbar0+qbar0)./(h0.^2-pbar0.*qbar0);
