%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class:     Psych 221/EE 362
% File:      WaveAberrationPSF
% Author:    Patrick Maeda
% Purpose:   Calculate and Plot PSF of Wave Aberration Function
% Date:      03.08.03	
%	
% Matlab 6.1:  03.09.03
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file calculates and plots the PSF of the Wave Aberration function specified by:
% jmax = highest mode in Wave Aberration Expansion
% Wrmsj = rms Wave Aberration coefficients of modes 0 to jmax
% d = pupil diameter in mm
% lambda = wavelength in nm
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
jmax=14                  %[INPUT] highest mode number from single indexing scheme

disp('Pupil Diameter (mm), RMS Wave Aberration Coefficients (micron), and Wavelength (nm)')
d=5.4;                   %[INPUT] pupil diameter in mm (3 to 8 mm)
PupilDiameter=d
Wrmsj=[0 0 0 0.4164 0 0.135 0.074 -0.092 0.011...  %[INPUT] rms wavefront error coefficient in microns
       -0.12 -0.038 0.016 0.085 -0.06 0.047]'
Wrmst=0;
for j=0:jmax
   Wrmst=Wrmst+Wrmsj(j+1)^2;
end
Wrmstotal=sqrt(Wrmst)    %total rms wavefront error in um

lambda=570               %[INPUT] wavelength in nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert units for calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wrmsj=Wrmsj*1e-3;          %rms wavefront error coefficients in mm
Wrmstotal=Wrmstotal*1e-3;  %total rms wavefront error in mm
lambda=lambda*1e-6;        %wavelength in mm
dw=d/lambda;               %pupil diameter in number of wavelengths
PRw=0.5*dw;                %pupil radius in number of wavelengths
apw=pi*PRw^2;              %pupil area in wavelength^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set-up x,y grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xwmin=-25000;  %minimum x-coordinate in number of wavelengths
xwmax=25000;   %maximum x-coordinate in number of wavelengths
ywmin=-25000;  %minimum y-coordinate in number of wavelengths
ywmax=25000;   %maximum y-coordinate in number of wavelengths
dxw=250;       %x-coordinate pixel width in number of wavelengths
dyw=250;       %y-coordinate pixel width in number of wavelengths

xw=xwmin:dxw:xwmax;   %x-coordinates in number of wavelengths
yw=ywmin:dyw:ywmax;   %y-coordinates in number of wavelengths
Imax=length(xw);
Jmax=length(yw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set-up circular pupil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for I=1:Imax    
    for J=1:Jmax
       P(I,J)=(sqrt(xw(I)^2+yw(J)^2) <= PRw);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Wave Aberration function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W=zeros(length(xw),length(yw));
for j=0:jmax
   n=ceil((-3+sqrt(9+8*j))/2);   %highest power or order of the radial polynomial term
   m=2*j-n*(n+2);                %azimuthal frequency of the sinusoidal component
   W=W+Wrmsj(j+1)*zernike(n,m,xw,yw,dw);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Compute PSF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PSF0=fft2(P)/apw;
PSF=fft2((P.*exp(-i*2*pi*W/lambda)))/apw;
PSF0=fftshift(PSF0);
PSF=fftshift(PSF);
PSF0=PSF0.*conj(PSF0);
PSF0=rot90(PSF0);
PSF=PSF.*conj(PSF);
PSF=rot90(PSF);
PSF0=flipud(PSF0);
PSF=flipud(PSF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot PSF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

umin=-0.002;    %minimum thetax-coordinate in radians
umax=0.002;     %maximum thetax-coordinate in radians
vmin=-0.002;    %minimum thetay-coordinate in radians
vmax=0.002;     %maximum thetay-coordinate in radians
du=0.00002;     %thetax-coordinate pixel width in radians
dv=0.00002;     %thetay-coordinate pixel width in radians
u=umin:du:umax;   %thetax-coordinates in radians
v=vmin:dv:vmax;   %thetay-coordinates in radians

figure
subplot(2,1,1)
scale=(2)^7/max(max(PSF0));
image(u*1000,v*1000,PSF0*scale)  %scaled for saturated display of image
%imagesc(u*1000,v*1000,PSF0)
axis image
xlabel('\theta_{x} (mrad)')
ylabel('\theta_{y} (mrad)')
axis square
axis xy
title(['PSF of Zero Aberration System, Pupil Diameter = ',...
        num2str(PupilDiameter), 'mm'], 'FontSize', 10);
colormap gray

%figure
subplot(2,1,2)
scale=(2)^7/max(max(PSF));
image(u*1000,v*1000,PSF*scale)   %scaled for saturated display of image
%imagesc(u*1000,v*1000,PSF)
axis image
axis xy
xlabel('\theta_{x} (mrad)')
ylabel('\theta_{y} (mrad)')
axis square
title(['PSF of Aberrated System, RMS Wavefront Error = ',num2str(Wrmstotal/lambda),'\lambda'],'FontSize', 10);
colormap gray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Plot PSF Cross-sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1)
plot(u*1000, PSF0((Imax+1)/2,:)/max(max(PSF0)))
xlabel('\theta_{x} (mrad)')
ylabel('Normalized Amplitude')
axis([umin*1000 umax*1000 0 1])
axis square
title(['PSF of Zero Aberration System, PD = ',...
        num2str(PupilDiameter), 'mm'], 'FontSize', 10);
subplot(2,2,2)
plot(v*1000, PSF0(:, (Jmax+1)/2)/max(max(PSF0)))
xlabel('\theta_{y} (mrad)')
ylabel('Normalized Amplitude')
axis([vmin*1000 vmax*1000 0 1])
axis square
title(['PSF of Zero Aberration System, PD = ',...
        num2str(PupilDiameter), 'mm'], 'FontSize', 10);
subplot(2,2,3)
plot(u*1000, PSF((Imax+1)/2,:)/max(max(PSF)))
xlabel('\theta_{x} (mrad)')
ylabel('Normalized Amplitude')
axis([umin*1000 umax*1000 0 1])
axis square
title(['PSF of Aberrated System, Wrms = ', num2str(Wrmstotal/lambda),'\lambda'],'FontSize', 10);
subplot(2,2,4)
plot(v*1000, PSF(:, (Jmax+1)/2)/max(max(PSF)))
xlabel('\theta_{y} (mrad)')
ylabel('Normalized Amplitude')
axis([vmin*1000 vmax*1000 0 1])
axis square
title(['PSF of Aberrated System, Wrms = ', num2str(Wrmstotal/lambda),'\lambda'],'FontSize', 10);
