% 1. Dado el siguiente sistema se pide:
% 1.1. Determinar las constantes de error de posición, velocidad y aceleración ( K p ,Kv , Ka )
% 1.2. Determinar el error de estado estable para los distintos tipos de entradas.

%Sistema 1
close all; clear all; clc

syms s

G = 50 / ((1+0.5*s)*(1+2*s));

kp=G;
kv=simplify(s*G);
ka=simplify((s^2)*G);

s=0;
kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)

%%
%Sistema 2
close all; clear all; clc

syms s

G = 2 / s*((1+0.1*s)*(1+0.5*s));

kp=G;
kv=simplify(s*G);
ka=simplify((s^2)*G);

s=0;
kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)

%% Sistema 3

close all; clear all; clc

syms s

G = 1 / s*((s^2)+(4*s)+200);

kp=G;
kv=simplify(s*G);
ka=simplify((s^2)*G);

s=0;
kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)

%% Sistema 4

close all; clear all; clc

syms s

G = (30*(1+2*s)*(1+4*s)) / s*((s^2)+(2*s)+10);

kp=G;
kv=simplify(s*G);
ka=simplify((s^2)*G);

s=0;
kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)

%% Sistema 5

close all; clear all; clc

syms s

G = (10*(1+s)) / s*(s+4)*(4*(s^2)+(6*s)+1);

kp=G;
kv=simplify(s*G);
ka=simplify((s^2)*G);

s=0;
kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)

%% Sistema 6

kp=K;

%% Sistema 7

close all; clear all; clc

syms s

G = (10*(1+s)) / (s^2)*(s+5)*(s+6);

kp=G;
kv=simplify(s*G);
ka=simplify((s^2)*G);

s=0;
kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)

%% Sistema 8

close all; clear all; clc

syms s

G = (10*(1+s)) / (s^3)*((s^2)+5*s+5);

kp=G;
kv=simplify(s*G);
ka=simplify((s^2)*G);

s=0;
kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)

%% Sistema 12

close all; clear all; clc

syms s

G = (3*s) /( (s^2)*(s+6));

kp=G;
kv=simplify(s*G);
ka=simplify((s^2)*G);

s=0;
kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)


%% Ejercicio 3

close all; clear all; clc

syms s k kt

g1=100/((0.2*s)+1);
g2=1/(20*s);

g1lc=g1/(1+(g1*kt));
g=collect(k*g1lc*g2, s);

kp=g;
kv=simplify(s*g);
ka=simplify((s^2)*g);


s=0;
kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)

k=50;
kt=5;


kp=eval(kp)
kv=eval(kv)
ka=eval(ka)
%Escalon
ess_e=1/(1+kp)
%Rampa
ess_r=(1/kv)
%Parabola
ess_p=(1/ka)

%% Ejercicio 5

close all; clear all; clc

s=tf('s');

g=3000/((s^2) + 156.25*s + 1837.5);
step(g,2)