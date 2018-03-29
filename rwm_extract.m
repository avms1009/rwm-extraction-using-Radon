close all;
clear all;
clc;

c=2.9979*10^8;  % Propagation velocity
fo=4.5e9;       % Carrier frequency (4.5GHz)
landa=c/fo;     % Wavelength (60cm for fo = 4.5e9)
theta_inc=pi/4; % Loook angle
Vplat=200;         % Velocity of platform
Rcenter=2000;  % Range distance to center of target area
l=0.5;           % Antenna vertical length actual
L=1.5;             %Antenna azimuth length

Swath=(landa*Rcenter)/(l*cos(theta_inc));


Rnear=Rcenter-Swath/2;     % Range Near
Rfar=Rcenter+Swath/2;   % Rang Far
Tr=(2*Swath)/(c*9);     %Pulse duration=>window length is 10 Tr
Ts=(2*Rnear)/c-Tr/2;             % Start time of sampling
Tf=(2*Rfar)/c+Tr/2;           % End time of sampling

Bw=100e6;%therefore resolution is 1.5m
Kr=Bw/Tr;%Modulation Rate
Fs=2*Bw;
dt=1/Fs;
t=Ts:dt:Tf;
%%%%%target location%%%%%%%%%%%%%%%%
R0=Rcenter;
Timg=(landa*R0)/(L*Vplat);
eta0=Timg/4;
%%%%%%%%%%%%%%%%%%

max_PRF=(l*c)/(2*Rcenter*landa*tan(theta_inc));
min_PRF=2*Vplat/L;
%PRF=max_PRF+min_PRF/2;
PRF=600;
Dta=1/PRF;
eta=-Timg/2:Dta:Timg/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Va=10; %along track velocity(m/s)
Vc=50; %cross track velocity(m/s)
Aa=0.1; %along track acceleration(m/s^2)
Ac=0.1; %cross track acceleration(m/s^2)
X=R0*sin(theta_inc);
Vr=X*Vc/R0;
Ar=X*Ac/R0;  
% Target Estimation Parameter
Vr_estimated=Vr;
Va_estimated=Va;
Ar_estimated=Ar;



R=R0-Vr*(eta-eta0)+( ( ((Vplat-Va)^2-Ar*R0)*((eta-eta0).^2) )/(2*R0) );
Tdelay=2*R/c;
figure;
plot(eta,Tdelay);
hold on;

title('delay of received echo versus slow time');

counter=1;
Sig1=zeros(length(eta),length(t));
for Td=Tdelay
    Sig1(counter,:)=exp(-1j*2*pi*fo*Td)*exp( 1j*Kr*pi*((t-Td).^2) ) .* (-Tr/2+Td<=t & t<=Tr/2+Td);
    counter=counter+1;
end



figure;
imagesc(abs(Sig1));
figure;
mesh(abs(Sig1));

T=-Tr/2:dt:Tr/2;
St=exp(1j*Kr*pi*T.^2).*(-Tr/2<=T & T<=Tr/2);
len=length(t)+length(T)-1;
%Nfft=2^nextpow2(len);
Nfft=len;
St_f=fft(St,Nfft);
Ht=exp(-1j*Kr*pi*T.^2).*(-Tr/2<=T & T<=Tr/2);
H_f=fft(Ht,Nfft);

Sig2=zeros(length(eta),Nfft);

for counter=1:length(eta)
    Sig2(counter,:)=ifft( fft(Sig1(counter,:),Nfft) .* H_f , Nfft );
end

% s_fdc_elimination=exp(-1j*4*pi*Vr_estimated*eta/landa);
% s_fdc_elimination=conj(s_fdc_elimination');
% s_fdc_elimination=repmat(s_fdc_elimination,1,Nfft);
% Sig2=Sig2 .* s_fdc_elimination;

figure;
imagesc(abs(Sig2));
figure;
mesh(abs(Sig2));

RWM_max=Vr*Timg/2;
deltaT=2*RWM_max/c;
Nsample=deltaT/dt;
%Q=round(pi/(atan(PRF*Timg*c*dt/4*RWM_max)+pi/2));
Q=180;
Y=My_Hough(abs(Sig2),Q);
[a,b]=find(Y==max(max(Y)));
slope=-cot((pi/Q)*b);

RWM_max_estimated=-PRF*Timg*c*dt*slope/4;
Vr_estimated=RWM_max_estimated*2/Timg;