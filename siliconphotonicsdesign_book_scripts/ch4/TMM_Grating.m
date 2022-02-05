function Grating
%This file is used to plot the reflection/transmission spectrum.

% Grating Parameters
Period=310e-9;  % Bragg period
NG=200;    % Number of grating periods
L=NG*Period;    % Grating length
width0=0.5;     % mean waveguide width
dwidth=0.01;    % +/-  waveguide width
width1=width0 - dwidth;
width2=width0 + dwidth;
loss_dBcm=3;    % waveguide loss, dB/cm
loss=log(10)*loss_dBcm/10*100;

% Simulation Parameters:
span=30e-9;  % Set the wavelength span for the simultion
Npoints = 10000;

% from MODE calculations
switch 1
    case 1 % Strip waveguide; 500x220 nm, TE, 1550 nm
    	% neff_wavelength is the wavelength dependant effective index, found from MODE, in this case for a 500 nm waveguide at 1550 nm TE.
        neff_wavelength = @(w) 2.4379 - 1.1193 * (w*1e6-1.554) - 0.0350 * (w*1e6-1.554).^2; % 500x220 oxide strip waveguide
	% dneff_width is the grating strength parameter, expressed as delta n (rather than kappa).
	% it is defined as the grating strength relative to the 500 nm waveguide, e.g., for w = 490, the function returns the strength of the grating
	% for a grating that alternates between 490 and 500.
	% the TMM code below uses n2-n1, where n2 is the wider section (e.g., 510), and n1 is the narrower section (e.g., 490).
	% note that this function was calculated using eigenmodes, which doesn't give the correct reflection coefficient or grating strength.
	% a more accurate grating strength can be found from band structure simulations (e.g., FDTD Bloch mode approach), or experimental data.
        dneff_width = @(w) 10.4285*(w-0.5).^3 - 5.2487*(w-0.5).^2 + 1.6142*(w-0.5);
end

% Find Bragg wavelength using lambda_Bragg = Period * 2neff(lambda_bragg); 
% Assume neff is for the average waveguide width.
f = @(lambda) lambda - Period*2* (neff_wavelength(lambda)+(dneff_width(width2)+dneff_width(width1))/2);
wavelength0 = fzero(f,1550e-9);

wavelengths=wavelength0 + linspace(-span/2, span/2, Npoints);
n1=neff_wavelength(wavelengths)+dneff_width(width1);  % low index
n2=neff_wavelength(wavelengths)+dneff_width(width2);  % high index

[R,T]=TMM_Grating_RT(wavelengths, Period, NG, n1, n2, loss);
figure;
plot (wavelengths*1e6,[R, T],'LineWidth',3); hold all
plot ([wavelength0, wavelength0]*1e6, [0,1],'--');  % calculated bragg wavelength
xlabel('Wavelength [\mum]')
ylabel('Response');
axis tight;
%printfig ('PS-WBG')


function [R,T]=TMM_Grating_RT(wavelength, Period, NG, n1, n2, loss)
%Calculate the R and T versus wavelength
M=TMM_Grating_Matrix(wavelength, Period, NG, n1, n2, loss);
q=length(wavelength);
T=abs(ones(q,1)./squeeze(M(1,1,:))).^2;
R=abs(squeeze(M(2,1,:))./squeeze(M(1,1,:))).^2;


function T=TMM_Grating_Matrix(wavelength, Period, NG, n1, n2, loss)
% Calculate the total transfer matrix of the gratings
l=Period/2;
T_hw1=TMM_HomoWG_Matrix(wavelength,l,n1,loss);
T_is12=TMM_IndexStep_Matrix(n1,n2);
T_hw2=TMM_HomoWG_Matrix(wavelength,l,n2,loss);
T_is21=TMM_IndexStep_Matrix(n2,n1);
q=length(wavelength);
Tp=zeros(2,2,q); T=Tp; 
for i=1:length(wavelength)
	Tp(:,:,i)=T_hw2(:,:,i)*T_is21(:,:,i)*T_hw1(:,:,i)*T_is12(:,:,i);
	T(:,:,i)=Tp(:,:,i)^NG; % 1st order uniform Bragg grating

	% for an FP cavity, 1st order cavity, insert a high index region, n2.
	T(:,:,i)=Tp(:,:,i)^NG * (T_hw2(:,:,i))^1 * Tp(:,:,i)^NG * T_hw2(:,:,i);
end

function T_hw=TMM_HomoWG_Matrix(wavelength,l,neff,loss)
% Calculate the transfer matrix of a homogeneous waveguide.
beta=2*pi*neff./wavelength-1i*loss/2; %Complex propagation constant
T_hw=zeros(2,2,length(neff));
T_hw(1,1,:)=exp(1i*beta*l);	
T_hw(2,2,:)=exp(-1i*beta*l);


function T_is=TMM_IndexStep_Matrix(n1,n2)
% Calculate the transfer matrix for a index step from n1 to n2.
T_is=zeros(2,2,length(n1));
a=(n1+n2)./(2*sqrt(n1.*n2));
b=(n1-n2)./(2*sqrt(n1.*n2));
%T_is=[a b; b a];
T_is(1,1,:)=a;	T_is(1,2,:)=b;
T_is(2,1,:)=b; T_is(2,2,:)=a;


function printfig (pdf)
FONTSIZE=20;
set ( get(gca, 'XLabel'),'FontSize',FONTSIZE)
set ( get(gca, 'YLabel'),'FontSize',FONTSIZE)
set (gca, 'FontSize',FONTSIZE); box on;
print ('-dpdf', pdf); system([ 'pdfcrop ' pdf ' ' pdf '.pdf' ]);


