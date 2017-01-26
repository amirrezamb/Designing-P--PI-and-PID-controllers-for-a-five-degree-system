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
%ku ba sahih o khata
%%% calculations
 Kpz2=3048.513;
 Wpz2=5.4513;
 pu=2*pi/Wpz2;
 syspz2=feedback(Kpz2*sysp,1);
 figure(21)
 step(syspz2)
 %%%designing the controller
 % maqadir ra baraye behtar shodan taqir midahim
 Kpiz2=Kpz2*0.45/2.5;
 tipiz2=pu*0.83*3.5;
 sysZpi=tf([Kpiz2 Kpiz2*tipiz2],[1 0]);
 sysZo=sysZpi*sysp;
 sysfZpi= feedback(sysZo,1);
 figure(22)
 step(sysfZpi)
 %%% info
infop2=stepinfo(sysfZpi)
yZpend2=step(sysfZpi);
endp2=yZpend2(end);
essp2=1-endp2
[Gmp2,Pmp2] = margin(sysZo)
 syspu2=sysZpi*(1-sysfZpi);
 [up2]=step(syspu2,10);
 Lup2(1)=up2(1);
  for i=1:1:247;
      Lup2(i+1)=Lup2(i)+up2(i+1)*up2(i+1)*0.04032;
  end
IEp2=Lup2(248)
figure(23);
pzplot(sysfZpi);
%%%realplant response
Grz=feedback(Gr*sysZpi,1);
figure(29);
step(Grz);
%%%noises
YR2=step(feedback(sysfZpi,1),time);
YD2=a*step(feedback(sysp,sysZpi),time);
YN2=lsim(feedback(-1*sysfZpi,-1),UN,time);
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

kpi=150;ti=9; 
syspi=tf([kpi ti*kpi],[1 0]);
sysR=syspi*sysp;
SysfR=feedback(sysR,1);
figure(31)
rlocus(sysR);
sgrid(zeta,Wn);

figure(32)
step(SysfR)
%%% info
infop3=stepinfo(SysfR)
yZpend3=step(SysfR);
endp3=yZpend3(end);
essp3=1-endp3
[Gmp3,Pmp3] = margin(sysR)
 syspu3=syspi*(1-SysfR);
 [up3]=step(syspu3,10);
 Lup3(1)=up3(1);
  for i=1:1:189;
      Lup3(i+1)=Lup3(i)+up3(i+1)*up3(i+1)*0.04032;
  end
  IEp3=Lup3(190)
figure(33);
pzplot(SysfR);
%%%realplant response
Grr=feedback(Gr*syspi,1);
figure(39);
step(Grr);
%%%noises
YR3=step(feedback(SysfR,1),time);
YD3=a*step(feedback(sysp,syspi),time);
YN3=lsim(feedback(-1*SysfR,-1),UN,time);
figure(38)
subplot(3,1,1)
plot(YR3+YD3)
subplot(3,1,2)
plot(YR3+YN3)
subplot(3,1,3)
plot(YR3+YD3+YN3)


%% Frequency Response

 kFpi=72;
 tiFpi=20;
 sysFpi=tf([kFpi kFpi*tiFpi],[1 0]);
 sysFo=sysFpi*sysp;
 sysfFpi= feedback(sysFo,1);
 figure(41)
 margin(sysFo)
 
 figure(42)
 step(sysfFpi)
%%% info
infop4=stepinfo(sysfFpi)
yZpend4=step(sysfFpi);
endp4=yZpend4(end);
essp4=1-endp4
[Gmp4,Pmp4] = margin(sysFo)
syspiu4=sysFpi*(1-sysfFpi);
[upi4]=step(syspiu4,10);
Lupi4(1)=upi4(1);
  for i=1:1:247;
      Lupi4(i+1)=Lupi4(i)+upi4(i+1)*upi4(i+1)*0.04032;
  end
  IEpi4=Lupi4(248)
  figure(43);
  pzplot(sysfFpi);
  %%%realplant response
Grf=feedback(Gr*sysFpi,1);
figure(49);
step(Grf);
%%%noises
YR4=step(feedback(sysfFpi,1),time);
YD4=a*step(feedback(sysp,sysFpi),time);
YN4=lsim(feedback(-1*sysfFpi,-1),UN,time);
figure(48)
subplot(3,1,1)
plot(YR4+YD4)
subplot(3,1,2)
plot(YR4+YN4)
subplot(3,1,3)
plot(YR4+YD4+YN4)

%% all
figure(51);
step(sysfZpi);
hold on;
step(sysfFpi);
hold on;
step(SysfR);