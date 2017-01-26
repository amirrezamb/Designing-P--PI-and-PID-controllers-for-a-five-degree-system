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
%% %% %% Part-1: Plant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1) two state space models
%%% tf plant
 Np=[1 a];
 Dp=[1 28 298 1554 4401 6318];
 sysp=tf(Np,Dp);

% it could also been done using the 2 following lines
% s = tf('s');  
% Gs = (s+a)/((s+b)*((s+c)*(s+c))*(s-p1)*(s-p2))

%% first system
 A1=[-28 1 0 0 0;
     -298 0 1 0 0;
     -1554 0 0 1 0;
     -4401 0 0 0 1;
     -6318 0 0 0 0];
 B1=[0;0;0;1;a];
 C1=[1 0 0 0 0];
 D1=[0];
 sysSS1=ss(A1,B1,C1,D1);
 %%% transfer function for first system
 [bb1,aa1] = ss2tf(A1,B1,C1,D1);
 sysL1=tf(bb1,aa1);
 
 %%% second system
  A2=[0 1 0 0 0;
      0 0 1 0 0;
      0 0 0 1 0;
      0 0 0 0 1;
      -6318 -4401 -1554 -298 -28];
 B2=[0;0;0;0;1];
 C2=[a 1 0 0 0];
 D2=[0];
 sysSS2=ss(A2,B2,C2,D2);
 %%% transfer function for second system
 [bb2,aa2] = ss2tf(A2,B2,C2,D2);
 sysL2=tf(bb2,aa2);
 
 
%% 2)eigenvalues & vectores
%%% first system
% e1 = eig(A1)
[V1,D1] = eig(A1);
for i=1:5
    mode1 (i)=exp( D1(i,i));
end
%%% second system
% e2 = eig(A2)
[V2,D2] = eig(A2);
for i=1:5
    mode2 (i)=exp( D2(i,i));
end
%% 3)reponses
%%% first system
%%%% step
[Yinitial1] = initial(sysSS1,x0,15);
[YUstep1] = step(sysSS1,15);
y11=Yinitial1+YUstep1;
figure(11)
plot(y11)
SteadyStateError1=y11(end);
%%%% impulse
[YUimpulse1]=impulse(sysSS1,15);
y12=Yinitial1+YUimpulse1;
figure(12)
plot(y12)
%%%% frequency
%time----->>setting the parameters
UF=a*sin(b*time);
[yFLs1]=lsim(sysSS1,UF,time,x0);
y13=yFLs1;
figure(13)
lsim(sysSS1,UF,time,x0)
%%% second system
%%%% step
[Yinitial2] = initial(sysSS2,x0,15);
[YUstep2] = step(sysSS2,15);
y21=Yinitial2+YUstep2;
figure(21)
plot(y21)
SteadyStateError2=y21(end);
%%%% impulse
[YUimpulse2]=impulse(sysSS2,15);
y22=Yinitial2+YUimpulse2;
figure(22)
plot(y22)
%%%% frequency
%time----->>setting the parameters
UF2=a*sin(b*time);
[yFLs2]=lsim(sysSS2,UF2,time,x0);
y23=yFLs2;
figure(23)
lsim(sysSS2,UF2,time,x0)
%% 4)Routh
syms krh
eqn = ((1071.90-0.095*krh)* (0.82*krh+4175.36)-242.5*(5*krh+6318))/(1071.90-0.095*krh) == 0;
solKrh = solve(eqn,krh)
% z=+ (15*194821601759227699797813005273572441^(1/2))/856519558037504 - 4014323550079942485/856519558037504
%% 5) root lucas
[rloc,kloc]=rlocus(sysp);
[xxx,yyy]=find(real(rloc)<0.1 & real(rloc)>-0.00001);
Wuloc=rloc(4,26);
Kuloc=kloc(26);
figure(55)
rlocus(sysp);
%% 6) bode & nyquist
Kbode=5000;
figure(1)
bode(sysp)
[Gm,Pm,Wgm,Wpm] = margin(sysp);
figure(2)
nyquist(sysp)
figure(3)
syspprime=tf(Kbode*Np,Dp);
bode(syspprime)
[Gm,Pm,Wgm,Wpm] = margin(syspprime);
figure(4)
nyquist(syspprime)
%% %% %% End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% help
% syms x
% eqn = (x+b)*((x+C)^2)*(X^2+) == 1;
% solx = solve(eqn,x)
% Np=[
% stepinfo(T);
