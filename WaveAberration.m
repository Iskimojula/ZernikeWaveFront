%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class:     Psych 221/EE 362
% File:      WaveAberration
% Author:    Patrick Maeda
% Purpose:   Calculate and Plot Wave Aberration
% Date:      03.08.03	
%	
% Matlab 6.1:  03.09.03
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file calculates and plots the Wave Aberration specified by:
% jmax = highest mode in Wave Aberration Expansion
% Wrmsj = rms Wave Aberration coefficients of modes 0 to jmax
%
% The Zernike Polynomial definitions used are derived from:
% Thibos, L., Applegate, R.A., Schweigerling, J.T., Webb, R., VSIA Standards Taskforce Members,
% "Standards for Reporting the Optical Aberrations of Eyes"
% OSA Trends in Optics and Photonics Vol. 35, Vision Science and its Applications,
% Lakshminarayanan,V. (ed) (Optical Society of America, Washington, DC 2000), pp: 232-244. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wave Aberration definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Maximum Mode Number j')
jmax=14             %[INPUT] highest mode number from single indexing scheme

disp('RMS Wave Aberration Coefficients (micron), and Wavelength (nm)')
d=2;                %normalized pupil diameter

Wrmsj=[0 0 0 0.4164 0 0.135 0.074 -0.092 0.011...  %[INPUT] rms wavefront error coefficient in microns
       -0.12 -0.038 0.016 0.085 -0.06 0.047]'
 
lambda=570          %[INPUT] wavelength in nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Set-up normalized x,y grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xn=-1:0.02:1;   %normalized x-coordinates
yn=-1:0.02:1;   %normalized y-coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Wave Aberration function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W=zeros(length(xn),length(yn));
for j=0:jmax
   n=ceil((-3+sqrt(9+8*j))/2);   %highest power or order of the radial polynomial term
   m=2*j-n*(n+2);                %azimuthal frequency of the sinusoidal component
   W=W+Wrmsj(j+1)*zernike(n,m,xn,yn,d);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Wave Aberration function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
%subplot(2,1,1)
W=rot90(W);
W=flipud(W);
imagesc(xn,yn,W);
colormap gray
axis xy
axis square
%set(gca, 'TickDir', 'out')
title(['Wave Aberration Function'],'FontSize', 10);
xlabel('Normalized x pupil coordinate');
ylabel('Normalized y pupil coordinate');
%colormap gray
colormap jet

figure
%subplot(2,1,2)
meshc(xn,yn,W)
%view(-37.5,45) 
title(['Wave Aberration Function'],'FontSize', 10);
xlabel('Normalized x pupil coordinate');
ylabel('Normalized y pupil coordinate');
zlabel('RMS Wavefront Error ( \mum)');
%colormap('default')
colormap([0.8 0 0])

