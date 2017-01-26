%% preparing
clear all;
close all;
clc;

%% setting the parameters
a=5;b=6;c=9;
p1=-2+3i;
p2=-2-3i;
x0=[b;c;0;0;a];
time=0:0.010225:15;
Ku=3048.513;
Wu=5.4513;
UN=c*cos(b*time);
%% %% %% Part-1: Plant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% tf plant
 Np=[1 a];
 Dp=[1 28 298 1554 4401 6318];
 sysp=tf(Np,Dp);
%%% real plant
s = tf('s');  
Gr = (s+a)/((s+b)*((s+c)*(s+c))*(s^2+4*s+16.38));

%% %% %% Part-2: Controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% zigler-2

 Kpz2=3048.513;
 Wpz2=5.4513;
 pu=2*pi/Wpz2;
 syspz2=feedback(Kpz2*sysp,1);

 figure(21)
 step(syspz2)
 % maqadir ra baraye behtar shodan taqir midahim
 Kpidz2=Kpz2*0.6/1;
 tipidz2=pu*0.5*5;
 tdpidZ2=pu*0.125*2.5;
 
 sysZpid=tf([Kpidz2*tdpidZ2 Kpidz2 Kpidz2*tipidz2],[1 0]);
 sysZo=sysZpid*sysp;
 sysfZpid= feedback(sysZo,1);
 figure(22)
 step(sysfZpid)
  %%% info
infop2=stepinfo(sysfZpid)
yZpend2=step(sysfZpid);
endp2=yZpend2(end);
essp2=1-endp2
[Gmp2,Pmp2] = margin(sysZo)
% syspu2=feedback(sysZpid,sysp);
% [up2]=step(syspu2,10);
%  Lup2(1)=up2(1);
%   for i=1:1:151;
%       Lup2(i+1)=Lup2(i)+up2(i+1)*0.04032;
%   end
%   IEp2=Lup2(152);
figure(23);
pzplot(sysfZpid);
%%%realplant response
Grz=feedback(Gr*sysZpid,1);
figure(29);
step(Grz);
%%%noises
YR2=step(feedback(sysfZpid,1),time);
YD2=a*step(feedback(sysp,sysZpid),time);
YN2=lsim(feedback(-1*sysfZpid,-1),UN,time);
figure(28)
subplot(3,1,1)
plot(YR2+YD2)
subplot(3,1,2)
plot(YR2+YN2)
subplot(3,1,3)
plot(YR2+YD2+YN2)



%% Root Lucas
syms cc;
eqncc= 0.18==exp(-cc*pi/((1-cc^2)^0.5));
cc1= solve(eqncc,cc);
zeta=(((pi*1i - 2*log(3) + log(50))*(pi*1i + 2*log(3) - log(50)))^(1/2)*(log(3)*2i - log(50)*1i)*1i)/(- log(3)*log(50)*4i + pi^2*1i + log(3)^2*4i + log(50)^2*1i);
syms ww;
eqnww = 3==(pi-acos(zeta))/(ww*(1-zeta^2)^0.5);
Wn=0.7862531315;
alpha=acos(zeta);

kpid=650;ti=45/14;td=1/14 ;
syspid=tf([kpid*td kpid ti*kpid],[1 0]);
sysR=syspid*sysp;
SysfR=feedback(sysR,1);
figure(31)
rlocus(sysR);
% sgrid(zeta,Wn);

figure(32)
step(SysfR)

%%% info
infop3=stepinfo(SysfR)
yZpend3=step(SysfR);
endp3=yZpend3(end);
essp3=1-endp3
[Gmp3,Pmp3] = margin(sysR)
%  syspu3=sfeedback(syspid,sysp);
%  [up3]=step(syspu3,10);
%  Lup3(1)=up3(1);
%   for i=1:1:151;
%       Lup3(i+1)=Lup3(i)+up3(i+1)*0.04032;
%   end
%   IEp3=Lup3(152);
figure(33);
pzplot(SysfR);
%%%realplant response
Grr=feedback(Gr*syspid,1);
figure(39);
step(Grr);
%%%noises
YR3=step(feedback(SysfR,1),time);
YD3=a*step(feedback(sysp,syspid),time);
YN3=lsim(feedback(-1*SysfR,-1),UN,time);
figure(38)
subplot(3,1,1)
plot(YR3+YD3)
subplot(3,1,2)
plot(YR3+YN3)
subplot(3,1,3)
plot(YR3+YD3+YN3)


%% Frequency Response

 kFpid=130;
 tiFpid=13;
 tdFpid=4;
 sysFpid=tf([kFpid*tdFpid kFpid tiFpid*kFpid],[1 0]);
 sysFo=sysFpid*sysp;
 sysfFpid= feedback(sysFo,1);
 figure(41)
 margin(sysFo)
 
 figure(42)
 step(sysfFpid)
 %%% info
infop4=stepinfo(sysfFpid)
yZpend4=step(sysfFpid);
endp4=yZpend4(end);
essp4=1-endp4
[Gmp4,Pmp4] = margin(sysFo)
% syspiu4=sysFpi*(1-sysfFpi);
% [upi4]=step(syspiu4,10);
% Lupi4(1)=upi4(1);
%   for i=1:1:247;
%       Lupi4(i+1)=Lupi4(i)+upi4(i+1)*0.04032;
%   end
%   IEpi4=Lupi4(152);
  figure(43);
  pzplot(sysfFpid);
  %%%realplant response
  Grf=feedback(Gr*sysFpid,1);
figure(49);
step(Grf);
%%%noises
YR4=step(feedback(sysfFpid,1),time);
YD4=a*step(feedback(sysp,sysFpid),time);
YN4=lsim(feedback(-1*sysfFpid,-1),UN,time);
figure(48)
subplot(3,1,1)
plot(YR4+YD4)
subplot(3,1,2)
plot(YR4+YN4)
subplot(3,1,3)
plot(YR4+YD4+YN4)


%% all
figure(51);
step(sysfZpid);
hold on;
step(sysfFpid);
hold on;
step(SysfR);
