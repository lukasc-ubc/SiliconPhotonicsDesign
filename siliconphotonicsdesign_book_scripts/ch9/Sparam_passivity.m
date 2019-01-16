function VectorFitting_Sparam ()
clear; close all; FONTSIZE=20;

file='../DC_ringmod_type1_R=10,gap=180,Lc=0,wg=500,lambda=1550,mesh=2,angle=30.mat';
load (file)

c=3e8; wavelength=c./f*1e6;
fSCALING = 1e14; f=f/fSCALING; s=1i*2*pi*f;  Np=length(f);
hwait = waitbar(0,'Please wait...');

% average the S parameters for the symmetric parameters that were simulated twice,
% by considering amplitude and phase separately
S1221 = (abs(S12)+abs(S21))/2 .* exp ( 1i * (unwrap(angle(S12)) + unwrap(angle(S21)) ) /2 );
S12=S1221; S21 = S1221;
S4132 = (abs(S41)+abs(S32))/2 .* exp ( 1i * (unwrap(angle(S41)) + unwrap(angle(S32)) ) /2 );
S41=S4132; S32 = S4132;
S13=S31; S23=S41; S33=S11; S43=S21; S14=S32; S24=S42; S34=S12; S44=S22;


for i=1:Np
  Sparam = [ [ S11(i),S12(i),S13(i),S14(i)];
    [ S21(i),S22(i),S23(i),S24(i)];
    [ S31(i),S32(i),S33(i),S34(i)];
    [ S41(i),S42(i),S43(i),S44(i)] ] ;
  Test1(i) = norm(Sparam);
  Sparam_w (:,:,i)=Sparam;
end

% Check if S-Parameters are already passive; if not, perform Vector Fit,
% otherwise, export.
if ~isempty(find(Test1>1))
  
  % rational fit, sweep number of parameters, N:
  optsN=20:1:100; 
  opts.poletype='lincmplx'; opts.parametertype='S'; opts.stable=0; opts.Niter_out=10;
  Npassive=[]; rms3passive=[];
  for i=1:length(optsN)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opts.N=optsN(i);    % Vector Fit:
    [SER,rmserr,Hfit,opts2]=VFdriver(Sparam_w,s,[],opts);
    
    %%%%%%%%%%%%%%%%%%%%%%  Enforce Passivity:
    [SER,H_passive,opts3,wintervals]=RPdriver(SER,s,opts);
    % Note: added output parameter wintervals to RPdriver.
    
    % Find rms error:
    tell=0; Nc=length(SER.D);
    for col=1:Nc
      for row=col:Nc % makes assumption that S = S'
        tell=tell+1; % make a single vector:
        Sparam10(tell,:)= squeeze(Sparam_w(row,col,:)).';
        fit10(tell,:)= squeeze(Hfit(row,col,:)).';
        fitP(tell,:)= squeeze(H_passive(row,col,:)).';
      end
    end
    rms2(i) = sqrt(sum(sum(abs((Sparam10-fit10).^2))))/sqrt(4*Np);
    rms3(i) = sqrt(sum(sum(abs((Sparam10-fitP).^2))))/sqrt(4*Np);
    
    % Determine if the Passivity Enforcement was successful.
    if (rms3(i) < 1) && isempty(wintervals)
      Npassive=[Npassive optsN(i)]; rms3passive = [rms3passive rms3(i)];
    end
    
    waitbar (i/length(optsN), hwait); 
    if (rms3(i) < 1e-4) && isempty(wintervals) ; break; end
  end
  close (hwait);
  
  % Plot error versus number of fitting parameters:
  figure;
  semilogy(optsN(1:i),rms2,'ro-', 'LineWidth',2, 'MarkerSize',7); hold all;
  labels={}; labels{end+1} = 'Vector Fit, rms';
  plot (Npassive, rms3passive, 'kx', 'MarkerSize',14, 'LineWidth',3);
  labels{end+1} = 'Passivity-Enforced Fit, rms';
  xlabel('Fitting order, N');   ylabel('rms error');
  for ii=1:length(labels); labels{ii}=[labels{ii} ' ' char(31) ]; end;
  legend (labels,'Location','NorthEast');
  printfig(file,'convergence');
  
  % Plot residue values:
  figure;
  residues=sort(squeeze(prod(prod(abs(SER.R),1),2).^(1/4)), 1, 'descend');
  semilogy(residues,'-s', 'LineWidth',3, 'MarkerSize',8); 
  xlabel('Residue index'); ylabel('Residue magnitude');
  printfig(file,'residues');
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot S-Parameters + fit functions:
  figure;
  Nc=length(SER.D);
  for row=1:Nc
    for col= row:Nc
      dum1=squeeze(Sparam_w(row,col,:));
      dum2=squeeze(H_passive(row,col,:));
      h1=semilogy(wavelength,abs(dum1),'b','LineWidth',3); hold on
      h2=semilogy(wavelength,abs(dum2),'r-.','LineWidth',4);
      h3=semilogy(wavelength,abs(dum2-dum1),'g--','LineWidth',3);
    end
  end
  hold off
  xlabel('Wavelength'); ylabel('Amplitude [S]');
  axis tight; Yl =ylim; ylim ([Yl(1)*10 Yl(2)*2]);
  labels={}; labels{end+1} = 'FDTD S-Parameters';
  labels{end+1} = 'Passivity-enforced Fit'; labels{end+1}='Deviation';
  for i=1:length(labels); labels{i}=[labels{i} ' ' char(31)]; end;
  legend (labels,'Location','Best');
  printfig(file,'VF2a');
  
  figure;
  Nc=length(SER.D);
  for row= 1:Nc  %1
    for col= row:Nc  %3
      dum1=squeeze(Sparam_w(row,col,:));
      dum2=squeeze(H_passive(row,col,:));
      h1=semilogy(wavelength,unwrap(angle(dum1)),'b','LineWidth',3); hold on
      h2=plot(wavelength,unwrap(angle(dum2)),'r-.','LineWidth',4);
      h3=plot(wavelength,abs(unwrap(angle(dum2))-unwrap(angle(dum1))),'g--','LineWidth',3);
    end
  end
  hold off
  xlabel('Wavelength'); ylabel('Phase [S]');
  axis tight; Yl =ylim; ylim ([Yl(1)*10 Yl(2)*2]);
  labels={}; labels{end+1} = 'FDTD S-Parameters';
  labels{end+1} = 'Passivity-enforced Fit'; labels{end+1}='Deviation';
  for i=1:length(labels); labels{i}=[labels{i} ' ' char(31)]; end;
  legend (labels,'Location','Best');
  printfig(file,'VF2p');
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Passivity test results:
  for i=1:length(f)
    Test0(i) = norm(Sparam_w(:,:,i));
    Test1(i) = norm(Hfit(:,:,i));
    Test2(i) = norm(H_passive(:,:,i));
  end
  figure;
  plot (wavelength,Test0,'LineWidth',2); hold all;
  plot (wavelength,Test1,'--','LineWidth',2)
  plot (wavelength,Test2,'LineWidth',4)
  plot(wavelength,ones(length(f),1),'--','LineWidth',3);
  labels={}; labels{1} = 'FDTD S-Parameters'; labels{2} = 'Rational Fit'; labels{3} = 'Passivity-enforced Fit'; labels{4}='Passivity limit';
  for i=1:length(labels); labels{i}=[labels{i} ' ' char(31)]; end; legend (labels,'Location','Best');
  axis tight
  xlabel ('Wavelength [\mum]');
  ylabel ('S-Parameter Passivity Test');
  printfig(file,'passivitytest3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% export S parameters to INTERCONNECT
fid = fopen([file '.sparam'],'w');  Nc=4;
for row= 1:Nc
  for col= 1:Nc
    fprintf(fid,'%s\n',[ '(''port ' num2str(col) ''',''TE'',1,''port ' num2str(row) ''',1,''transmission'')' ]  );
    fprintf(fid,'%s\n',[ '(' num2str(Np) ',3)' ]  );
    dum2=squeeze(H_passive(row,col,:));
    mag = abs(dum2);
    phase = unwrap (angle(dum2)); % figure; plot (wavelength, phase);
    for i=1:Np
      fprintf(fid,'%g	%g	%g\n', f(i)*1e14, mag(i), phase(i));
    end
  end
end
fclose(fid);


function printfig (file, b)
global PRINT_titles;
PRINT_titles=0;
FONTSIZE=20;
set(get(gca,'xlabel'),'FontSize',FONTSIZE);
set(get(gca,'ylabel'),'FontSize',FONTSIZE);
set(get(gca,'title'),'FontSize',FONTSIZE-5);
set(gca,'FontSize',FONTSIZE-2);
if PRINT_titles==0
  delete(get(gca,'title'))
end
%a=strfind(file,'.'); file(a)=',';
pdf = [file(1:end-4) '_' b '.pdf'];
print ('-dpdf','-r300', pdf);
system([ 'pdfcrop ' pdf ' ' pdf ' &' ]);
% system(['acroread ' pdf '.pdf &']);


