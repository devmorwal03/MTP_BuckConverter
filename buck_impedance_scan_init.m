%% 
clc
close all
clear all
tic

P=bodeoptions;
P.FreqUnits='Hz';
P.PhaseWrapping='on';
P.Grid='on';
w=logspace(log10(1e1*2*pi),log10(1e4*2*pi),10000);
s=tf('s');

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

kmff=0;
kuff=1;

Zin=-(((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/ ...
   (uRef0*(Gt*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) + (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0)));

ZinCC=-(((1 - Gt*kuff + Cf*Lf*s^2 + Gt*Ri*(Ru + Cf*s))*Vdc0^2)/(uRef0*(Gt*iLf0*(kmff - kuff + Cf*kmff*Lf*s^2 + Ri*(Ru + Cf*s)) + Cf*(-1 + Gt*kmff)*s*uRef0)));

Zout1=s*Lg/(1+s^2*Lg*Cdc);

% Zin=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*Vdc0^2)/(uRef0*((-Gt)*iL0*Ri*(1 + Rload*(Ru + Cf*s)) + uRef0 + Cf*Rload*s*uRef0));
% Zin=-((Rload*(Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*Vdc0^2)/(((-Rload)*(1 + Cf*Rload*s) + Gt*Ri*(1 + Rload*(Ru + Cf*s)))*uRef0^2));
% Zinuff=-((Rload*(Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - Rload + Ri*Rload*(Ru + Cf*s)))*Vdc0^2)/(((-Rload)*(1 + Cf*Rload*s) + Gt*(Ri - Rload + Ri*Rload*(Ru + Cf*s)))*uRef0^2));
figure();
bode(Zin,w,P);
grid on;
hold on;
bode(ZinCC,w,P);
legend('Z_{in,uff}','Z_{in,CC}');
title('Input impedance of converter');

ZinEval=freqresp(Zin,w);
for i=1:length(ZinEval)
    ZinEval2(i)=ZinEval(1,1,i);
end

%%

 Thp=1/(2*pi*200);
 Tlp=1/(2*pi*2000);
 kad=0;
 Gad=kad*1/30*s*Thp/(1+s*Thp)*1/(1+s*Tlp);

 kmff=0;kuff=1;
 ZinAD=((Rload + Lf*s + Cf*Lf*Rload*s^2 + Gt*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s)))*uRef0*Vdc0^2)/ ...
   (uRef0^2*((-Gt)*iLf0*(Ri - kuff*Rload + Ri*Rload*(Ru + Cf*s) + kmff*(Rload + Lf*s + Cf*Lf*Rload*s^2)) - (-1 + Gt*kmff)*(1 + Cf*Rload*s)*uRef0) + ...
    Gad*Gt*Ri*(iLf0*(Rload + Lf*s + Cf*Lf*Rload*s^2) + uRef0 + Cf*Rload*s*uRef0)*Vdc0^2);

ZinADEval=freqresp(ZinAD,w);
for i=1:length(ZinADEval)
    ZinADEval2(i)=ZinADEval(1,1,i);
end

figure();
bode(Zin,w,P);
grid on;
hold on;
bode(ZinAD,w,P);
% bode(Zout1,w,P);

legend('Z_{in,uff}','Z_{in,uff,ad}');
title('Input impedance of converter with active damping');

Gadd=c2d(Gad,Tsamp,'tustin');
[NumAD,DenAD]=tfdata(Gadd,'v');

% 
% 
% 


%% perturbance
f_vec = logspace(1, 4, 40);

Ts = 500e-9;
Fs = 1/Ts;

Zin = zeros(1, length(f_vec));

%% method 1, frequency by frequency
%for droop = 0:1
    for k = 1:length(f_vec)

        f_test = f_vec(k);
% 
%         if (f_test < 100) % should be further improved to get better SNR
%             Vm = 1;
%         else
%             Vm = 0.1;
%         end

        Vm=20;

        w_test = 2*pi*f_test;

        t0 = 0.1; %startup
        N_periods = 10; % in order to reduce spectral leakage when calc. DFTs
        tend = t0 + N_periods/f_test;
        
        out = sim('buck_impedance_scan.slx', tend);
        t = out.buck.Time;
        ind0 = find(t>=t0, 1);

        i_in = out.buck.Data(ind0:end,1);
        v_in = out.buck.Data(ind0:end,2);

        N = length(i_in);
        freq = Fs*(0:((N - mod(N,2))/2))/N;

        fft_sig = fft(i_in);
        P2 = abs(fft_sig)/N;
        P1 = P2(1 : (N - mod(N,2))/2 + 1);
        P1(2:end-1) = 2*P1(2:end-1); 

        fft_sig_v = fft(v_in);
        P2_v = abs(fft_sig_v)/N;
        P1_v = P2_v(1 : (N - mod(N,2))/2 + 1);
        P1_v(2:end-1) = 2*P1_v(2:end-1); 

        delta_f = freq - f_test;
        index_f = find(abs(delta_f) == min(abs(delta_f))); %find index

        amp_i = P1(index_f); %amplitude of the current at f_vec(k)
        phase_i = angle(fft_sig(index_f)); %phase of the current at f_vec(k)

        amp_v = P1_v(index_f); %amplitude of the current at f_vec(k)
        phase_v = angle(fft_sig_v(index_f)); %phase of the current at f_vec(k)    
        Zin(k) = (amp_v*exp(1i*phase_v))/(amp_i*exp(1i*phase_i)); 
        k
    end
%end

figure();
subplot(2,1,1)
set(gcf, 'Position',  [100, 100, 1000, 700])
semilogx(f_vec, 20*log10(abs(Zin)), 'c*'); 
hold on;
semilogx(w/(2*pi), 20*log10(abs(ZinADEval2)));
ylim([0 60])
ylabel('|Z_{in}| [dB\Omega]')
%legend('Location','NorthEast');
title('Buck converter input impedance')
grid on
subplot(2,1,2)
semilogx(f_vec, (atan2(imag(Zin), real(Zin))*180/pi),  'c*'); 
hold on;
semilogx(w/(2*pi), (atan2(imag(ZinADEval2), real(ZinADEval2))*180/pi));
ylim([-190 190])
yticks([-180 -90 0 90 180])
xlabel('f [Hz]')
ylabel('arg(Z_{in}) [\circ]')
grid on

toc