%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class:     Psych 221/EE 362
% File:      ZernikePolynomialMTF
% Author:    Patrick Maeda
% Purpose:   Calculate and Plot MTF of the Wave Aberration Function
% Date:      03.04.03	
%	
% Matlab 6.1:  03.04.03
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file calculates and plots the MTF of the Wave Aberration function specified by:
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

xwmin=-15000;  %minimum x-coordinate in number of wavelengths
xwmax=15000;   %maximum x-coordinate in number of wavelengths
ywmin=-15000;  %minimum y-coordinate in number of wavelengths
ywmax=15000;   %maximum y-coordinate in number of wavelengths
dxw=150;       %x-coordinate pixel width in number of wavelengths
dyw=150;       %y-coordinate pixel width in number of wavelengths

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
% Compute PSF, OTF, and MTF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PSF0=fft2(P)/apw;
PSF=fft2((P.*exp(-i*2*pi*W/lambda)))/apw;
PSF0=PSF0.*conj(PSF0);
PSF=PSF.*conj(PSF);
OTF0=fft2(PSF0);
OTF0=OTF0/max(max(OTF0));
OTF=fft2(PSF);
OTF=OTF/max(max(OTF));
OTF0=fftshift(OTF0);
OTF=fftshift(OTF);
MTF0=abs(OTF0);
MTF0=rot90(MTF0);
MTF=abs(OTF);
MTF=rot90(MTF);
MTF0=flipud(MTF0);
MTF=flipud(MTF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot MTF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sxmin=-15000;    %minimum sx-coordinate in radians
sxmax=15000;     %maximum sx-coordinate in radians
symin=-15000;    %minimum sy-coordinate in radians
symax=15000;     %maximum sy-coordinate in radians
dsx=150;     %sx-coordinate pixel width in radians
dsy=150;     %sy-coordinate pixel width in radians
sx=sxmin:dsx:sxmax;   %sx-coordinates in radians
sy=symin:dsy:symax;   %sy-coordinates in radians

figure
subplot(2,1,1)
%mesh(sx*pi/180,sy*pi/180,MTF0)
contour(sx*pi/180,sy*pi/180,MTF0)
axis([-50 50 -50 50])
xlabel('s_{x} (cycle/deg)')
ylabel('s_{y} (cycle/deg)')
axis square
%zlabel('Contrast')
%set(gca,'XDir','rev','YDir','rev') 
title(['MTF of Zero Aberration System, ',...
        num2str(PupilDiameter), ' mm pupil'], 'FontSize', 10);
colormap([0.8 0 0])

%figure
subplot(2,1,2)
%mesh(sx*pi/180,sy*pi/180,MTF)
contour(sx*pi/180,sy*pi/180,MTF)
axis([-50 50 -50 50])
xlabel('s_{x} (cycle/deg)')
ylabel('s_{y} (cycle/deg)')
axis square
%zlabel('Contrast')
title(['MTF of Aberrated System, RMS Wavefront Error = ',num2str(Wrmstotal/lambda),'\lambda'],'FontSize', 10);
colormap([0.8 0 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Plot MTF Cross-sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sxaxismax=0.15*sxmax*pi/180;
%syaxismax=0.15*symax*pi/180;
sxaxismax=50;
syaxismax=50;

figure
subplot(2,2,1)
plot(sx*pi/180,MTF0((Imax+1)/2,:))
xlabel('s_{x} (cycle/deg)')
%ylabel('Contrast')
axis([0 sxaxismax 0 1])
axis square
title(['MTF of Zero Aberration System, ',...
        num2str(PupilDiameter), 'mm pupil'], 'FontSize', 10);
subplot(2,2,2)
plot(sy*pi/180,MTF0(:, (Jmax+1)/2))
xlabel('s_{y} (cycle/deg)')
%ylabel('Contrast')
axis([0 syaxismax 0 1])
axis square
title(['MTF of Zero Aberration System, ',...
        num2str(PupilDiameter), 'mm pupil'], 'FontSize', 10);
subplot(2,2,3)
plot(sx*pi/180,MTF((Imax+1)/2,:))
xlabel('s_{x} (cycle/deg)')
%ylabel('Contrast')
axis([0 sxaxismax 0 1])
axis square
title(['MTF of Aberrated System, Wrms = ', num2str(Wrmstotal/lambda),'\lambda'],'FontSize', 10);
subplot(2,2,4)
plot(sy*pi/180,MTF(:, (Jmax+1)/2))
xlabel('s_{y} (cycle/deg)')
%ylabel('Contrast')
axis([0 syaxismax 0 1])
axis square
title(['MTF of Aberrated System, Wrms = ', num2str(Wrmstotal/lambda),'\lambda'],'FontSize', 10);

