%% buck_control_design.m
% This file is for the controller design of the buck converter stage
% calculation of input impedance


clear all;
close all;

P=bodeoptions;
P.FreqUnits='Hz';
P.PhaseWrapping='on';
P.Grid='on';
w=logspace(2,5,10000);
s=tf('s');

%operating point
Rload=20;
Vdc0=350;
Uout0=200;
Iload0=Uout0/Rload
iLf0=Iload0;
uRef0=Uout0;

% filter paramters
fs=24e3;              % switching frequency
fsamp=24e3;          % sampling frequency
Lf=1/fs*1/(0.4*Iload0)*(1-0.5)*Vdc0     
Cf=1/(Lf*(0.1*2*pi*fs)^2)

Cdc=163e-6;
rCdc=0.2%200e-3;
Ldc=0%0.05e-3;
rLdc=0%20e-3;

% total plant delay
%Td=0.5*1/fs+0.5*1/fsamp;   % delay of double update mode + fast sampling
Td=0.75/fs;
[num, den]=pade(Td,2);  % approximation for dead-time
Gt=tf(num, den);

% current controller design
kuff=1;
Gi=(Cf*Gt*s)/(1 - Gt*kuff + Cf*Lf*s^2);

wBWi=2*pi*fs*0.1;
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


kmff=0;kuff=0;
Zin=-(((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/ ...
   (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0)));
kmff=1;kuff=0;
Zinmff=-(((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/ ...
   (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0)));
kmff=0;kuff=1;
Zinuff=-(((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/ ...
   (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0)));
kmff=1;kuff=1;
Zinmffuff=-(((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/ ...
   (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0)));
% Zin=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*Vdc0^2)/(uRef0*((-Gt)*iL0*Ri*(1 + Rload*(Ru + Cf*s)) + uRef0 + Cf*Rload*s*uRef0));
% Zin=-((Rload*(Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*Vdc0^2)/(((-Rload)*(1 + Cf*Rload*s) + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*uRef0^2));
% Zinuff=-((Rload*(Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/(((-Rload)*(1 + Cf*Rload*s) + Gt*(Ri - Rload + Ri*Rload*(Ru + Cf*s)))*uRef0^2));
figure();
bode(Zin,w,P);
grid on;
hold on;
bode(Zinmff,w,P);
bode(Zinuff,w,P);
bode(Zinmffuff,w,P);
legend('Z_{in}','Z_{in,mff}','Z_{in,uff}','Z_{in,mffuff}');
title('Input impedance of converter');



%% effect of DC capacitance
kmff=0;kuff=0;
Zin=rLdc + Ldc*s + 1/((Cdc*s)/(1 + Cdc*rCdc*s) - (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + ...
       (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0))/((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2));
kmff=1;kuff=0;
Zinmff=rLdc + Ldc*s + 1/((Cdc*s)/(1 + Cdc*rCdc*s) - (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + ...
       (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0))/((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2));
kmff=0;kuff=1;
Zinuff=rLdc + Ldc*s + 1/((Cdc*s)/(1 + Cdc*rCdc*s) - (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + ...
       (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0))/((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2));
kmff=1;kuff=1;
Zinmffuff=rLdc + Ldc*s + 1/((Cdc*s)/(1 + Cdc*rCdc*s) - (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + ...
       (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0))/((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2));

figure();
bode(Zin,w,P);
grid on;
hold on;
bode(Zinmff,w,P);
bode(Zinuff,w,P);
bode(Zinmffuff,w,P);
legend('Z_{in}','Z_{in,mff}','Z_{in,uff}','Z_{in,mffuff}');
title('Input impedance incl. Cdc');


%% effect of current control bandwidth

kuff=1;
Gi=(Cf*Gt*s)/(1 - Gt*kuff + Cf*Lf*s^2);

wBWi=2*pi*fs*0.1;
Tni=Cf*Lf/Td;
kpi=wBWi*Lf;
Ri=kpi*(1+s*Tni)/(s*Tni);

Giol=Ri*Gi;
Gicl=feedback(Giol,1);


Gu=1/(s*Cf)*Gicl;

wBWu=0.1*wBWi;
kpu=wBWu*Cf;
Tnu=10*1/(wBWu);
Ru=kpu*(1+s*Tnu)/(s*Tnu);

Guol=Ru*Gu;


kmff=0;kuff=1;
Zinuff=-(((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/ ...
   (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0)));
% Zin=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*Vdc0^2)/(uRef0*((-Gt)*iL0*Ri*(1 + Rload*(Ru + Cf*s)) + uRef0 + Cf*Rload*s*uRef0));
% Zin=-((Rload*(Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*Vdc0^2)/(((-Rload)*(1 + Cf*Rload*s) + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*uRef0^2));
% Zinuff=-((Rload*(Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/(((-Rload)*(1 + Cf*Rload*s) + Gt*(Ri - Rload + Ri*Rload*(Ru + Cf*s)))*uRef0^2));
figure();
bode(Zinuff,w,P);
grid on;
hold on;


kuff=1;
Gi=(Cf*Gt*s)/(1 - Gt*kuff + Cf*Lf*s^2);

wBWi=2*pi*fs*0.033;
Tni=Cf*Lf/Td;
kpi=wBWi*Lf;
Ri=kpi*(1+s*Tni)/(s*Tni);

Giol=Ri*Gi;
Gicl=feedback(Giol,1);


Gu=1/(s*Cf)*Gicl;

wBWu=0.02*wBWi;
kpu=wBWu*Cf;
Tnu=10*1/(wBWu);
Ru=kpu*(1+s*Tnu)/(s*Tnu);

Guol=Ru*Gu;


kmff=0;kuff=1;
Zinuff=-(((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/ ...
   (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0)));
% Zin=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*Vdc0^2)/(uRef0*((-Gt)*iL0*Ri*(1 + Rload*(Ru + Cf*s)) + uRef0 + Cf*Rload*s*uRef0));
% Zin=-((Rload*(Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*Vdc0^2)/(((-Rload)*(1 + Cf*Rload*s) + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*uRef0^2));
% Zinuff=-((Rload*(Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/(((-Rload)*(1 + Cf*Rload*s) + Gt*(Ri - Rload + Ri*Rload*(Ru + Cf*s)))*uRef0^2));
bode(Zinuff,w,P);
grid on;
hold on;
legend('Z_{in,uff} wBWi=0.1*ws','Z_{in,uff} wBWi=0.02*ws');
title('Input impedance of converter, impact of current control bandwidth');


%% active damping
kuff=1;
Gi=(Cf*Gt*s)/(1 - Gt*kuff + Cf*Lf*s^2);

wBWi=2*pi*fs*0.1;
Tni=Cf*Lf/Td;
kpi=wBWi*Lf;
Ri=kpi*(1+s*Tni)/(s*Tni);

Giol=Ri*Gi;
Gicl=feedback(Giol,1);


Gu=1/(s*Cf)*Gicl;

wBWu=0.1*wBWi;
kpu=wBWu*Cf;
Tnu=10*1/(wBWu);
Ru=kpu*(1+s*Tnu)/(s*Tnu);

Guol=Ru*Gu;
Thp=1/(2*pi*200);
Tlp=1/(2*pi*2000);
Gad=1/100*s*Thp/(1+s*Thp)*1/(1+s*Tlp);

kmff=0;kuff=0;
Zin=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*uRef0*Vdc0^2)/ ...
  (uRef0^2*((-Gt)*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) - (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0) + ...
   Gad*Gt*Ri*(iLf0*(Rload + Lf*s + Cf*Lf*Rload*s^2) + uRef0 + Cf*Rload*s*uRef0)*Vdc0^2);
kmff=1;kuff=0;
Zinmff=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*uRef0*Vdc0^2)/ ...
  (uRef0^2*((-Gt)*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) - (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0) + ...
   Gad*Gt*Ri*(iLf0*(Rload + Lf*s + Cf*Lf*Rload*s^2) + uRef0 + Cf*Rload*s*uRef0)*Vdc0^2);
kmff=0;kuff=1;
Zinuff=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*uRef0*Vdc0^2)/ ...
  (uRef0^2*((-Gt)*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) - (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0) + ...
   Gad*Gt*Ri*(iLf0*(Rload + Lf*s + Cf*Lf*Rload*s^2) + uRef0 + Cf*Rload*s*uRef0)*Vdc0^2);
kmff=1;kuff=1;
Zinmffuff=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*uRef0*Vdc0^2)/ ...
  (uRef0^2*((-Gt)*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) - (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0) + ...
   Gad*Gt*Ri*(iLf0*(Rload + Lf*s + Cf*Lf*Rload*s^2) + uRef0 + Cf*Rload*s*uRef0)*Vdc0^2);

figure();
bode(Zin,w,P);
grid on;
hold on;
% bode(Zinmff,w,P);
bode(Zinuff,w,P);
% bode(Zinmffuff,w,P);
legend('Z_{in}','Z_{in,uff}');
title('Input impedance with active damping');


kmff=0;kuff=1;
Gad=0;
Zinuff0=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*uRef0*Vdc0^2)/ ...
  (uRef0^2*((-Gt)*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) - (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0) + ...
   Gad*Gt*Ri*(iLf0*(Rload + Lf*s + Cf*Lf*Rload*s^2) + uRef0 + Cf*Rload*s*uRef0)*Vdc0^2);

figure();
bode(Zinuff0,w,P);
grid on;
hold on;
% bode(Zinmff,w,P);
bode(Zinuff,w,P);
% bode(Zinmffuff,w,P);
legend('Z_{in,uff} without AD','Z_{in,uff} with AD');
title('Input impedance with active damping');






