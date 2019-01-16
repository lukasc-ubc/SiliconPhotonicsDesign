% Plot the ring modulator spectrum
% Wei Shi, UBC, 2012
% weis@ece.ubc.ca
%
% RingMod_spectrum_plot;

c=299792458;
lambda=1e-9*(1530:0.1:1560);
Filter_type='add-drop';
pn_angle=2*pi; %
r=10e-6; Lc=0; L_pn=2*pi*r*pn_angle/(2*pi); w=500e-9; % WG parameters
pn_offset=0; ds_n_plus=1e-6; ds_p_plus=1e-6; % pn-junction design
T=25; V0=0; % temperature and voltage

[Ethru0 Edrop0 Qi0 Qc0 tau_rt0 Rj0 Cj0] = RingMod_spectrum (lambda, Filter_type, r, Lc, L_pn, w, pn_offset, ds_n_plus, ds_p_plus, T, V0);

figure;
plot(lambda*1e9, [10*log10(abs(Ethru0).^2); 10*log10(abs(Edrop0).^2)], 'linewidth', 2);
xlim([min(lambda) max(lambda)]*1e9);
set(gca, 'fontsize', 14);
xlabel({'\lambda (nm)'}, 'fontsize', 14);
ylabel({'Transmission (dB)'}, 'fontsize', 14);
legend('Through', 'Drop');

% zoom at one peak wavelength
lambda_zoom=1e-9*(1540.7:0.0025:1541);
V=-4:1:0;
lenV = length(V); lenLZ = length(lambda_zoom);
Ethru=zeros(lenV, lenLZ); Edrop=zeros(lenV, lenLZ);
A=zeros(lenV, lenLZ);
Qi=zeros(lenV, lenLZ);    Qc=zeros(lenV, lenLZ);
Cj=zeros(lenV,1);         Rj=zeros(lenV,1);
for i=1:lenV
    [Ethru(i,:) Edrop(i,:) Qi(i,:) Qc(i,:) tau_rt(i,:) Rj(i,:) Cj(i,:)]=RingMod_spectrum(lambda_zoom, Filter_type, r, Lc, L_pn, w, pn_offset, ds_n_plus, ds_p_plus, T, V(i));
end

Qt=1./(1./Qi+1./Qc);% total Q
tp=Qt./(c/1541e-9*2*pi); % photon lifetime
tp_av=sum(tp, 2)/(length(lambda_zoom)); % average photon lifetime across over the spectrum
fcq=1./(2*pi*tp_av);
fcj=1./(2*pi*Rj.*Cj);
fc=1./(1./fcq+1./fcj);

figure; plot(lambda_zoom*1e9, 10*log10(abs(Ethru).^2), 'linewidth', 2);
set(gca, 'fontsize', 14);
xlabel({'\lambda (nm)'}, 'fontsize', 14);
ylabel({'Transmission (dB)'}, 'fontsize', 14);
legend({cat(2, num2str(-V'), char(ones(length(V),1)*'V'))}, 'Location', 'best', 'fontsize', 14);

if strcmp(Filter_type,'add-drop')
    figure;
    plot(lambda_zoom*1e9, 10*log10(abs(Edrop).^2), 'linewidth', 2);
    set(gca, 'fontsize', 14);
    xlabel({'\lambda (nm)'}, 'fontsize', 14);
    ylabel({'Transmission (dB)'}, 'fontsize', 14);
    legend({cat(2, num2str(-V'), char(ones(length(V),1)*'V'))}, 'Location', 'best', 'fontsize', 14);
end

figure;
plot(-V, [fcq fcj fc]*1e-9, 'linewidth', 2);
set(gca, 'fontsize', 14);
xlabel({'Voltage (V)'}, 'fontsize', 14);
ylabel({'Cutoff frequency (GHz)'}, 'fontsize', 14);
legend({'\tau_p determined','p-n junction determined','f_c'}, 'Location', 'NorthWest', 'fontsize', 14);
