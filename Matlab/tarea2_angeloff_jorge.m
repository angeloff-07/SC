%% Sistema 1: Sistema a lazo cerrado con realimentación tacométrica
clear all; close all; clc;

s = tf ('s');

kp = 31; kd = 16;

g1 = (53 * (s+1)) / ((s + 10) * (s + 100));

g2 = 1 / s;

ft_lc_int = feedback( g1 , kd );

ft_la = kp * ft_lc_int * g2;

ft_lc_ext = feedback ( ft_la , 1 );

step ( ft_lc_ext );


%%  Sistema 2: Sistema de control en cascada
clear all; close all; clc;


pi=(s+0.1)/s;

p=12;

g1=27/(s+200);

g2=5/(s+0.1);

ft_la_int=p*g1;

ft_lc_int=feedback(ft_la_int,1);

ft_la_ext=pi*ft_lc_int*g2;

ft_lc_ext=feedback(ft_la_ext,1);

step(ft_lc_ext);

%% Sistema 3: Sistema a lazo cerrado con controlador PID
clear all; close all; clc;

syms kp ti td wn psita s

g1=(((ti*td*(s^2)))+(ti*s)+1)/(ti*s);

g2=(wn^2)/((s^2)+(2*psita*wn*s)+(wn^2));

ft_la=kp*g1*g2;

ft_lc=ft_la/(1+ft_la)

pretty(ft_lc)

ft_lc=factor(ft_lc)

pretty(ft_lc)

%% Sistema 4: Sistema a lazo cerrado con controlador PID

syms k  t s kp ti td 

g1=k/((t*s)+1);

control= (1+(1/(ti))+(td*s))*kp;

ft_la=control*g1;

ft_lc=ft_la/(1+ft_la)

factor(ft_lc)

simplify(ft_lc)

pretty(ft_lc)

%% Ejercicio 3
% 3.1. Determinar las cuatro funciones de transferencia que modelan el sistema.

s=tf('s');
L=1e-6; R=2; Ka=0.042; J=10e-6; B=0.3e-5; Kb=0.042;

g1=1/((L*s)+R);

g2=1/((J*s)+B);

% % a) W_E con T=0
ft_la1=g1*Ka*g2;

ft_lc1=feedback(ft_la1, Kb)


% b) W_T
ft_la2=g2;
ft_lc2=feedback(ft_la2,(Ka*Kb*g1));
ft_lc2=ft_lc2*-1

% c) I_E
ft_la3=g1;
ft_lc3=feedback(ft_la3, (g2*Ka*Kb))

% d) I_T

ft_la4=g2*(-Kb)*g1;
ft_lc4=feedback(ft_la4,-Ka)

% 3.2. Simular la respuesta del sistema para e(t ) = 24u(t ) y t L (t ) = 0.01u (t − 2) .

% Simular la respuesta para e(t) = 24u(t) y tL(t) = 0.01u(t-2)
t = 0:0.01:2.5;
e_t = 24 * ones(size(t));  % Entrada escalón de 24V
tl_t = 0.01 * (t >= 2);    % Perturbación escalón a t=2s

% Respuesta de la velocidad W(s) y la corriente I(s)
% Puede verse que la respuesta real para cada salida corresponde a la suma
% de las entradas.

W_response = lsim(ft_lc1, e_t, t) + lsim(ft_lc2, tl_t, t);
I_response = lsim(ft_lc3, e_t, t) + lsim(ft_lc4, tl_t, t);

% Graficar los resultados
figure;
subplot(2,1,1);
plot(t, W_response,'LineWidth', 1.5);
title('Respuesta de la Velocidad W(t)');
xlabel('Tiempo [s]');
ylabel('Velocidad [rad/s]');

subplot(2,1,2);
plot(t, I_response,'LineWidth', 1.5);
title('Respuesta de la Corriente I(t)');
xlabel('Tiempo [s]');
ylabel('Corriente [A]');

%Cada plot representa correctamente los fenomenos, al arrancar el motor se
%produce un pico de corriente debido a los bobinados y a los dos segundos
%se produce una perturbacion que podria ser una carga, disminuyendo la
%velocidad del rotor.


%% Ejercicio 4
% 4.1. Determinar las funciones de transferencia W ( s ) / Wr ( s ) y W ( s ) / TL ( s ) .

Kp=20;
W_Wr=feedback((Kp*ft_lc1),1);

W_TL=-feedback(g2,(Ka*(Kp+Kb)*g1));

% 4.2. Simular la respuesta del sistema para Wr (t ) = 300u (t ) y t L (t ) = 0.01u (t − 2) .

% Simular la respuesta para wref(t) = 300u(t) y tL(t) = 0.01u(t-2)
t = 0:0.01:2.5;
e_t = 300 * ones(size(t));  % Entrada escalón de 300V
tl_t = 0.01 * (t >= 2);    % Perturbación escalón a t=2s

% Respuesta de la velocidad W(s) y la corriente I(s)
% Puede verse que la respuesta real para cada salida corresponde a la suma
% de las entradas.

W_response = lsim(W_Wr, e_t, t) + lsim(W_TL, tl_t, t);


% Graficar los resultados
figure;

plot(t, W_response,'LineWidth', 1.5);
title('Respuesta de la Velocidad W(t)');
xlabel('Tiempo [s]');
ylabel('Velocidad [rad/s]');

%Puede verse que ajustar la variable Kp=20 nos acercamos a los 300r/seg
%que tenemos de referencia 

%% Diagramas de Flujo de Señal y Álgebra de Mason
% 6. Demostrar que los siguientes sistemas son equivalentes.

close all; clear all; clc
% Definición de las funciones.
syms G H 

P=G; %caminos directos: 1
L1=-G*H; %lazos: 1
delta=1-L1; %Δ=1−(suma de todas las ganancias de bucles individuales)+(suma de las ganancias de combinaciones de dos bucles disjuntos)−(suma de las ganancias de combinaciones de tres bucles disjuntos)+…

delta1=1;

M=(P*delta1)/delta

%% Sistema 2.
P=G;
L1=-G*H;
delta=1-L1;
delt1=1;
M=(P*delta1)/delta

%% 7. Demostrar que los siguientes sistemas no son equivalentes.
clear all; close all; %clc;

% Sistema 1
syms G1 G2 G3 H1 H2 H3

P=G1*G2*G3;
L1=-G1*H1;
L2=-G2*H2;
L3=-G3*H3;

delta=1-(L1+L2+L3)+((L1*L2)+(L1*L3)+(L2*L3))-(L1*L2*L3);
delta1= 1;

M=(P*delta1)/delta;

M=simplify(M)

%% Sistema 2

P=G1*G2*G3;
L1=-G1*H1;
L2=-G2*H2;
L3=-G3*H3;

delta=1-(L1+L2+L3)+(L1*L3);
delta1=1;


M=(P*delta1)/delta;

M=simplify(M)
%Los sistemas no son equivalentes por el calculo de delta, estos son
%distintos

%% 10. Aplicar la Regla de Mason para encontrar las siguientes Funciones de Transferencia.
%Sistema 1
%Y5/Y1
syms G1 G2 G3 G4 H1 H2 H3

P1=G1*G2*G3;
P2=G4*G3;
L1=-G1*H1;
L2=-G3*H2;
L3=-G1*G2*G3*H3;
L4=-G4*G3*H3;

delta=1-(L1+L2+L3+L4)+(L1*L2);
delta1=1;
delta2=1;

M=((P1*delta1)+(P2*delta2))/delta;
M=simplify(M);
pretty(M)

%% Y4/Y1
syms G1 G2 G3 G4 H1 H2 H3

P1=G1*G2;
P2=G4;

L1=-G1*H1;
L2=-G3*H2;
L3=-G1*G2*G3*H3;
L4=-G4*G3*H3;

delta=1-(L1+L2)+(L1*L2);
delta1=1;
delta2=1;

M=((P1*delta1)+(P2*delta2))/delta;
M=simplify(M);
pretty(M)


%% Y2/Y1
syms G1 G2 G3 G4 H1 H2 H3

P1=1;

L1=-G1*H1;
L2=-G3*H2;
L3=-G1*G2*G3*H3;
L4=-G4*G3*H3;

delta=1-(L1+L2)+(L1*L2);
delta1=1-(L2);

M=(P1*delta1)/delta;
M=simplify(M);
pretty(M)

%% Sistema 2  

%Y5/Y1
syms G1 G2 G3 G4 H1 H2 H3 H4

P1=G1*G2*G3;
P2=G4*G3;
L1=-G1*H1;
L2=-G3*H2;
L3=-G1*G2*G3*H3;
L4=-G4*G3*H3;
L5=-H4;

delta=1-(L1+L2+L3+L4+L5)+(L1*L2+L1*L5);
delta1=1;
delta2=1;

M=((P1*delta1)+(P2*delta2))/delta;
M=simplify(M);
pretty(M)

%% Y4/Y1
syms G1 G2 G3 G4 H1 H2 H3 H4

P1=G1*G2;
P2=G4;

L1=-G1*H1;
L2=-G3*H2;
L3=-G1*G2*G3*H3;
L4=-G4*G3*H3;
L5=-H4;

delta=1-(L1+L2+L3+L4+L5)+(L1*L2+L1*L5);
delta1=1-L5;
delta2=1-L5;

M=((P1*delta1)+(P2*delta2))/delta;
M=simplify(M);
pretty(M)

%% Y2/Y1
syms G1 G2 G3 G4 H1 H2 H3 H4

P1=1;

L1=-G1*H1;
L2=-G3*H2;
L3=-G1*G2*G3*H3;
L4=-G4*G3*H3;
L5=-H4;

delta=1-(L1+L2+L3+L4+L5)+(L1*L2+L1*L5);

delta1=1-(L2+L5);


M=((P1*delta1)+(P2*delta2))/delta;
M=simplify(M);
pretty(M)

%% Sistema 3
% Y5/Y1
syms s

P1=10*(1/s)*(1/s);
P2=10*5;

L1=-5/s;
L2=-10;
L3=-(1/s)*(1/s);
L4=-5;
L5=-5*(10*s)*(-5);

delta=1-(L1+L2+L3+L4+L5);

delta1=1;
delta2=1;

M=((P1*delta1)+(P2*delta2))/delta;
M=simplify(M);
pretty(M)

%% Y4/Y1

syms s

P1=10*(1/s)*(1/s);
P2=10*5;

L1=-5/s;
L2=-10;
L3=-(1/s)*(1/s);
L4=-5;
L5=-5*(10*s)*(-5);

delta=1-(L1+L2+L3+L4+L5);

delta1=1;
delta2=1;

M=((P1*delta1)+(P2*delta2))/delta;
M=simplify(M);
pretty(M)

%% Y2/Y1

syms s

P1=10;

L1=-5/s;
L2=-10;
L3=-(1/s)*(1/s);
L4=-5;
L5=-5*(10*s)*(-5);

delta=1-(L1+L2+L3+L4+L5);

delta1=1-L2;

M=(P1*delta1)/delta;
M=simplify(M);
pretty(M)

%% Sistema 5
% Y5/Y1
syms s;
P1=(1/s)*(1/s)*30;
P2=5*(1/s);

L1=-(1/s);
L2=-30;
L3=-(1/s)*(1/s)*30;
L4=-5*(1/s);
L5=-10;

delta=1-(L1+L2+L3+L4+L5)+((L1*L5)+(L1*L4)*(L2*L5)+(L3*L5));

delta1=1-L5;
delta2=1-L1;

M=((P1*delta1)+(P2*delta2))/delta;
M= simplify(M);

pretty(M)

%% Y4/Y1
syms s;
P1=-(1/s)*(1/s);
P2=5*(1/s);

L1=-(1/s);
L2=-30;
L3=-(1/s)*(1/s)*30;
L4=-5*(1/s);
L5=-10;

delta=1-(L1+L2+L3+L4+L5)+((L1*L5)+(L1*L4)*(L2*L5)+(L3*L5));

delta1=1-L5;
delta2=1            ;

M=((P1*delta1)+(P2*delta2))/delta;
M= simplify(M);

pretty(M)

%% Y2/Y1

syms s;
P1=1;


L1=-(1/s);
L2=-30;
L3=-(1/s)*(1/s)*30;
L4=-5*(1/s);
L5=-10;

delta=1-(L1+L2+L3+L4+L5)+((L1*L5)+(L1*L4)*(L2*L5)+(L3*L5));

delta1=1-(L1+L2+L5);


M=((P1*delta1))/delta;
M= simplify(M);

pretty(M)


