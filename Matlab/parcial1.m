%% Ejercicio 1
close all; clear all; clc

syms g1 g2 g3 g4 g5 g6 h2 h6 h4

p1=g1*g2*g3*g4*g5;
p2=g1*g6;

l1=-g2*h2;
l2=-g6*h6;
l3=-g1*g2*g3*g4*h4;


delta=1-(l1+l2+l3)+((l1*l2)+(l1*l3));

delta1=1-l2;
delta2=1;

M=((p1*delta1)+(p2*delta2))/delta;
M= simplify(M);

pretty(M)

%% Ejercicio 2
close all; clear all; clc

syms s
G=6/((s^2)+350*s);
Ct=50;
C=30*s + 1440;

F_LA=G/(1+(G*Ct));

FT=simplify(C*F_LA);
FT=collect(FT, s);
pretty(FT)

% Sistema tipo 0 para una entrada  escalon ess=1/1+kp

kp=FT;

s=0;

kp=eval(kp)
ess=1/(1+kp)



%% Ejercicio 3
close all; clear all; clc
syms s D theta V


J = 0.0000020; %[Kgm2], Momento de inercia de la esfera

m = 250; %[g], Masa de la esfera

tau = 0.002; %[s], Constante de tiempo del servo

R = 3; %[cm], Radio de la esfera

L = 0.450; %[m], Longitud de la barra

K = 36; %[rad/v], Ganancia del servo

g = 9.81; %[m/s^2], Aceleración de la gravedad


eq1= (J/R^2)*D*(s^2) == -(m*g/L)*theta;
eq2= theta == (K/tau)*V;


[D_sol,theta_sol] = solve([eq1, eq2],[D, theta]);

D_V = simplify(D_sol/V)

pretty(D_V)  


%% Ejercicio 4
close all; clear all; clc

k=1.5;
ymax=2.15;
tp=0.05;

mp=(ymax-k)/k;

psita=sqrt(log(mp)^2 / (pi^2 + log(mp)^2))

wn=pi / (tp*sqrt(1-psita^2))

%% Ejercicio 5
close all; clear all; clc

s=tf('s');

G1=1.5/((0.03*s+1)^2);
G2=1.5/((0.06*s+1)^2);
G3=1.5/((0.015*s+1)^2);
G4=1.5/((0.015*s-1)^2);
G5=1.5/(0.015*s+1);
step(G3, 0.15)

ts=0.0875;

tau=(ts/4)