%% SooMa ((SO)und insertion l(O)ss of (M)ultilayer pipeline j(A)cket)
%% Author: Zibo Liu
%% Email: zibo@kth.se
%% Date: 2021-05-18
%% License: An Open Source Code, please cite Zibo's relevant research papers after using the code.
%% Note that there may be assumptions when describing the acoustic response when developing this code, typically working in the low frequency range (<= 2000 Hz)


clear
clc
close
%% parameter
addpath(genpath('./data/'))
addpath(genpath('./src/'))


parameter_pressure_acoustics
parameter12


%% three layers (pipe + porpous layer + impervious jacket)
t_jacket = 0.0005;
t_damping = 0.00;
t34 = t_jacket + t_damping;

r1 = 0.15; %  see figure illustration
r2 = r1+0.0045; % ??? 
r3 = r2+0.05; 
r4 = r3+t34;

t12 = r2-r1;


l = 6; % length of the lagged pipe

%% impedance of the bare pipe
I12 = rho12.*t12.*(1+(0.025/r1)^2);
C12 = r1^2/E12/t12;


f_ring = sqrt(E12/rho12/(1-nu12^2))./(2*pi.*r1);
f_critical = c0^2/2/pi*(I12/(E12*t12^3/(12*(1-nu12^2))))^(1/2);
Ve_Z12 = 1i.*Ve_omega.*I12.*(1-(Ve_freq./f_critical).^2-(f_ring./Ve_freq).^2 );

Ve_Z12 = (1i.*Ve_omega*I12 + 1./(1i.*Ve_omega.*C12)); % simplified 

%% wave vector propagating along z axis, for bulk reaction
Ve_kz = Ve_k.*(1-2i.*rho0.*c0./Ve_Z12./Ve_k./r1).^(1/2);

%% porous material wave number and characteristic impedance
% parameter23_porous

rp = 2500*(1+1i*0.001); 

Ve_kp = Ve_k.*(1-1i.*rp./Ve_omega./rho0).^(1/2);
Ve_krp = (Ve_kp.^2 - Ve_kz.^2).^(1/2);
Ve_kr = (Ve_k.^2 - Ve_kz.^2).^(1/2);
% Ve_kr1 = (Ve_k.*2i.*rho0.*c0./Ve_Z12./r1).^(1/2);

Ve_Yp = rho0.*c0.*(1-1i.*rp./Ve_omega./rho0).^(1/2);

%% radiation impedance of the outer sheet
Ve_Z04 = 1i.*Ve_omega.*rho0./Ve_kr.*besselh(0, 2, Ve_kr*r4)./Ve_kr/r4; % 0 is the order, 2 is the second kind of hankle
Ve_Z02 = 1i.*Ve_omega.*rho0./Ve_kr.*besselh(0, 2, Ve_kr*r2)./Ve_kr/r2; % 0 is the order, 2 is the second kind of hankle

%% porous material transfer matrix
Ve_J02 = besselj(0, Ve_krp.*r2);
Ve_J03 = besselj(0, Ve_krp.*r3);
Ve_J12 = besselj(1, Ve_krp.*r2);
Ve_J13 = besselj(1, Ve_krp.*r3);

Ve_N02 = bessely(0, Ve_krp.*r2);
Ve_N03 = bessely(0, Ve_krp.*r3);
Ve_N12 = bessely(1, Ve_krp.*r2);
Ve_N13 = bessely(1, Ve_krp.*r3);

Ve_X = 1i.*Ve_krp./Ve_Yp./Ve_kp;
Ve_det = Ve_X.*(Ve_J13.*Ve_N03 - Ve_J03.*Ve_N13);

Ve_uT11 = Ve_J13.*Ve_N02 - Ve_J02.*Ve_N13;
Ve_uT12 = 1./Ve_X.*(Ve_J03.*Ve_N02 - Ve_J02.*Ve_N03);
Ve_uT21 = Ve_X.*(Ve_J12.*Ve_N13 - Ve_J13.*Ve_N12);
Ve_uT22 = (Ve_J12.*Ve_N03 - Ve_J03.*Ve_N12);

Ve_T11 = Ve_X./Ve_det.*Ve_uT11;
Ve_T12 = Ve_X./Ve_det.*Ve_uT12;
Ve_T21 = Ve_X./Ve_det.*Ve_uT21;
Ve_T22 = Ve_X./Ve_det.*Ve_uT22;

%% impedance of the outer sheet with/without damping
parameter34_jacket

Ve_Z34 = 1i.*Ve_omega.*(rho_jacket.*t_jacket+rho_damping*t_damping)*(1+1i*0.2);

%% transmission loss
% Ve_ratio = Ve_Z04.*(Ve_T11 + Ve_Z12.*Ve_T21) + Ve_T11.*Ve_Z34 + Ve_T12 + Ve_Z12.*(Ve_T21.*Ve_Z34 + Ve_T22);
Ve_ratio_lagged = Ve_Z12.*Ve_T21.*Ve_Z34;
Si = pi*r1^2;
S4 = 2*pi*r4*l;
Ve_TLlagged = 10*log10( Si/(2*rho0*c0) ./(1/2.*real(Ve_Z04).*S4).* abs(Ve_ratio_lagged).^2);

Ve_ratio_bare = Ve_Z12;
S2 = 2*pi*r2*l;
Ve_TLbare = 10*log10( Si/(2*rho0*c0) ./(1/2.*real(Ve_Z02).*S2).* abs(Ve_ratio_bare).^2);


%% insertion loss
Ve_IL = Ve_TLlagged - Ve_TLbare;

%% plot

load('Measured_insertion_loss') 

figure(2)
Ma = [Ve_freq,Ve_IL];
[Ve_freq_octave, Ve_IL_octave] = fun_octave(Ma);
pl_the_octave = semilogx(Ve_freq_octave, Ve_IL_octave,'k:','linewidth',2.5);
hold on 
pl_mea = semilogx(Measured_insertion_loss.freq, Measured_insertion_loss.caseA,'r--','linewidth',3.5);
hold off
plotxlabel = xlabel('Frequency (Hz)'); set(plotxlabel,'FontSize',16, 'interpreter', 'latex');
plotylabel = ylabel('Insertion loss (dB)'); set(plotylabel,'FontSize',16,'interpreter', 'latex');
plotlegend= legend([  pl_the_octave pl_mea],...
                    'Transfer Matrix Method','Measured'); set(plotlegend,'Location','Best','FontSize',12,'box','off','interpreter','latex');
axis([100 4000 -20 60]);  set(gca,'TickLabelInterpreter','latex'),% set(gcf, 'units','points','position',[350 350 420 210]);
filename = 'IL_case'; grid; set(gca,'TickLabelInterpreter','latex')
% savefigure(path_png, path_eps, path_fig, filename)




