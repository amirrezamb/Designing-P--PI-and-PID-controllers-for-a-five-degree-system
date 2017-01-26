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


%% zigler-1
% figure(11)
% step(sysp,15)
% hold on;
% ypz1=step(sysp,15);
% gz(1)=0;
% for i=1:1:1466
%     gz(i+1)=(ypz1(i+1)-ypz1(i))/0.010225;
% end
% Rpz1=max(gz);
% ylinez54=ypz1(54);
% Lpz1=0.5419-(ylinez54/Rpz1);
% 
%     
%% zigler-2
%%% calculations
Kpz2=3048.513;
Wpz2=5.4513;
pu=2*pi/Wpz2;
syspz2=feedback(Kpz2*sysp,1);
figure(21)
step(syspz2)
%%% desining the controller
% maqadir ra baraye behtar shodan taqir midahim
kZp=Kpz2*0.5/1;
sysfZp= feedback(kZp*sysp,1);
figure(22)
step(sysfZp);
%%% info
infop2=stepinfo(sysfZp)
yZpend2=step(sysfZp);
endp2=yZpend2(end);
essp2=1-endp2
[Gmp2,Pmp2] = margin(kZp*sysp)

 syspu2=kZp*(1-sysfZp);
 [up2]=step(syspu2,10);
 Lup2(1)=up2(1);
  for i=1:1:247;
      Lup2(i+1)=Lup2(i)+up2(i+1)*up2(i+1)*0.04032;
  end
  IEp2=Lup2(248)
figure(23);
pzplot(sysfZp);
%%%realplant response
Grz=feedback(Gr*kZp,1);
figure(29);
step(Grz);
%%%noises
YR2=step(feedback(sysfZp,1),time);
YD2=a*step(feedback(sysp,kZp),time);
YN2=lsim(feedback(-1*sysfZp,-1),UN,time);
figure(28)
subplot(3,1,1)
plot(YR2+YD2)
subplot(3,1,2)
plot(YR2+YN2)
subplot(3,1,3)
plot(YR2+YD2+YN2)

%% Root Lucas
%%% calculations
syms cc;
eqncc= 0.18==exp(-cc*pi/((1-cc^2)^0.5));
cc1= solve(eqncc,cc);
zeta=(((pi*1i - 2*log(3) + log(50))*(pi*1i + 2*log(3) - log(50)))^(1/2)*(log(3)*2i - log(50)*1i)*1i)/(- log(3)*log(50)*4i + pi^2*1i + log(3)^2*4i + log(50)^2*1i);
syms ww;
eqnww = 3==(pi-acos(zeta))/(ww*(1-zeta^2)^0.5);
Wn=0.7862531315;
alpha=acos(zeta);
%%% desining the controller
kRp=1500;
sysfRp=feedback(kRp*sysp,1);
figure(31)
rlocus(sysp);
sgrid(zeta,Wn);
figure(32)
step(sysfRp)
%%% info
infop3=stepinfo(sysfRp);
yZpend3=step(sysfRp);
endp3=yZpend3(end);
essp3=1-endp3;
[Gmp3,Pmp3] = margin(kRp*sysp);
 syspu3=kRp*(1-sysfRp);
 [up3]=step(syspu3,10);
 Lup3(1)=up3(1);
  for i=1:1:247;
      Lup3(i+1)=Lup3(i)+up3(i+1)*up3(i+1)*0.04032;
  end
  IEp3=Lup3(248)
figure(33);
pzplot(sysfRp);
%%%realplant response
Grr=feedback(Gr*kRp,1);
figure(39);
step(Grr);
%%%noises
YR3=step(feedback(sysfRp,1),time);
YD3=a*step(feedback(sysp,kRp),time);
YN3=lsim(feedback(-1*sysfRp,-1),UN,time);
figure(38)
subplot(3,1,1)
plot(YR3+YD3)
subplot(3,1,2)
plot(YR3+YN3)
subplot(3,1,3)
plot(YR3+YD3+YN3)

%% Frequency Response
%%% desining the controller
 kFp=1214;
 sysFo=kFp*sysp;
 sysfFp= feedback(sysFo,1);
 figure(41)
 margin(sysFo)
 
 figure(42)
 step(sysfFp)
%%% info 
infop4=stepinfo(sysfFp)
yZpend4=step(sysfFp);
endp4=yZpend4(end);
essp4=1-endp4
[Gmp4,Pmp4] = margin(sysFo)
 syspu4=kFp*(1-sysfFp);
 [up4]=step(syspu4,10);
 Lup4(1)=up4(1);
  for i=1:1:247;
      Lup4(i+1)=Lup4(i)+up4(i+1)*up4(i+1)*0.04032;
  end
IEp4=Lup4(248)
figure(43);
pzplot(sysfFp);
%%%realplant response
Grf=feedback(Gr*kFp,1);
figure(49);
step(Grf);
%%%noises
YR4=step(feedback(sysfFp,1),time);
YD4=a*step(feedback(sysp,kFp),time);
YN4=lsim(feedback(-1*sysfFp,-1),UN,time);
figure(48)
subplot(3,1,1)
plot(YR4+YD4)
subplot(3,1,2)
plot(YR4+YN4)
subplot(3,1,3)
plot(YR4+YD4+YN4)


%% all
figure(51)
step(sysfZp,15);
hold on
step(sysfFp,15);
hold on
step(sysfRp,15);