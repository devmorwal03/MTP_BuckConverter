%% clc
close all
clear all

%% converter parameters
%operating point
Rload=10;
Vdc0=350;
Uout0=200;
D0=Uout0/Vdc0;
Iload0=Uout0/Rload
iLf0=Iload0;
uRef0=Uout0;
Imax=2*Iload0;

P=bodeoptions;
P.FreqUnits='Hz';
P.PhaseWrapping='on';
P.Grid='on';
w=logspace(2,5,10000);
s=tf('s');

% filter paramters
fs=48e3;              % switching frequency
fsamp=48e3;          % sampling frequency
Tsamp=1/fsamp;
Lf=0.5*1/(fs)*1/(0.2*Iload0)*(1-0.5)*Vdc0     
Cf=1/(Lf*(0.05*2*pi*fs)^2)

Cdc=163e-6;
rCdc=0.2%0.2%200e-3;

Lg=50e-6;
rLg=0;
iLg0=0;

f0grid=1/(2*pi*sqrt(Cdc*Lg))

% total plant delay
%Td=0.5*1/fs+0.5*1/fsamp;   % delay of double update mode + fast sampling
Td=1.75/fs;
[num, den]=pade(Td,2);  % approximation for dead-time
Gt=tf(num, den);

% current controller design
kuff=1;
Gi=(Cf*Gt*s)/(1 - Gt*kuff + Cf*Lf*s^2);

wBWi=2*pi*fs*0.07;
Tni=Cf*Lf/Td;
kpi=wBWi*Lf;
Ri=kpi*(1+s*Tni)/(s*Tni);

Giol=Ri*Gi;
Gicl=feedback(Giol,1);

figure();
bode(Gi,w,P);
grid on;
hold on;
bode(Ri,w,P);
bode(Giol,w,P);
legend('G_{iL}','R_{i}','R_{i}*G_{iL}');
title('compensated open-loop tf of converter current');

% voltage controller design
Gu=1/(s*Cf)*Gicl;

wBWu=0.1*wBWi;
kpu=wBWu*Cf;
Tnu=10*1/(wBWu);
Ru=kpu*(1+s*Tnu)/(s*Tnu);

Guol=Ru*Gu;

figure();
bode(Gu,w,P);
grid on;
hold on;
bode(Ru,w,P);
bode(Guol,w,P);
legend('G_{uC}','R_{u}','R_{u}*G_{uC}');
title('compensated open-loop tf of converter output voltage');