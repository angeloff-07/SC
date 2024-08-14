clear all; close all; clc;

pkg load symbolic

%6. Encontrar las Transformadas de Laplace de las siguientes funciones. Verificar utilizando Octave.

syms t a w

g1=laplace(dirac(t));

g2=laplace(heaviside(t));

g3=laplace(exp(-2*t));

g4=laplace(5*exp(-5*t));

g5=laplace(1-exp(-2*t));

g6=laplace((t*sin(2*t))+(3*exp(-10*t)));

g7=laplace(heaviside(t-2)*exp(-5*(t-2)));

g8=laplace(exp(-a*t)*cos(w*t));


% 7. Encontrar las Transformadas Inversas de Laplace de las siguientes funciones. Verificar utilizando Octave.

syms s

g1=ilaplace(2/(s+3));

%...

g11=ilaplace(((100*(s+2))/(s*((s^2)+4)*(s+1)))*exp(-s));

%8. Sean sistemas modelados por las siguientes funciones de transferencia. Determinar el
%valor final de las salidas para entradas escalón unitario usando la propiedad del
%Teorema del Valor Final. Simular la respuesta de los sistemas.

pkg load control
s=tf('s');

g1=tf(5,[1 2]); %TVF=2.5
step(g1);

g2=tf(1,[1 5 6]); %TVF=0.166
step(g2);

#g3 no pude usar ni inputdelay ni pade

g4=tf([5 5],[1 1 2]); %TVF=2.5
step(g4);

g5=(5*(s+2))/((s+3)*(s+4)); %TVF=0.83
step(g5, 10);

g6=5/s; %TVF=infinito

step(g6, 100);

g7=(12*(s+2))/(s*(s+4)); %TVF=infinito
step(g7);

g8=s/(s+40);  %TVF=infinito
step(g8);
