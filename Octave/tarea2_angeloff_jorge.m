clear all; close all; clc;

% pkg load symbolic
% pkg load control

s = tf ('s');

%Sistema 1: Sistema a lazo cerrado con realimentación tacométrica

% kp = 31; kd = 16;
% 
% g1 = (53 * (s+1)) / ((s + 10) * (s + 100));
% 
% g2 = 1 / s;
% 
% ft_lc_int = feedback( g1 , kd );
% 
% ft_la = kp * ft_lc_int * g2;
% 
% ft_lc_ext = feedback ( ft_la , 1 );
% 
% step ( ft_lc_ext );


% Sistema 2: Sistema de control en cascada

% pi=(s+0.1)/s;
% 
% p=12;
% 
% g1=27/(s+200);
% 
% g2=5/(s+0.1);
% 
% ft_la_int=p*g1;
% 
% ft_lc_int=feedback(ft_la_int,1);
% 
% ft_la_ext=pi*ft_lc_int*g2;
% 
% ft_lc_ext=feedback(ft_la_ext,1);
% 
% step(ft_lc_ext);

% Sistema 3: Sistema a lazo cerrado con controlador PID

% syms kp ti td wn psita s
% 
% g1=(((ti*td*(s^2)))+(ti*s)+1)/(ti*s);
% 
% g2=(wn^2)/((s^2)+(2*psita*wn*s)+(wn^2));
% 
% ft_la=kp*g1*g2;
% 
% ft_lc=ft_la/(1+ft_la)
% 
% pretty(ft_lc)
% 
% ft_lc=factor(ft_lc)
% 
% pretty(ft_lc)

% Sistema 4: Sistema a lazo cerrado con controlador PID

% syms k  t s kp ti td 
% 
% g1=k/((t*s)+1);
% 
% control= (1+(1/(ti))+(td*s))*kp;
% 
% ft_la=control*g1;
% 
% ft_lc=ft_la/(1+ft_la)
% 
% factor(ft_lc)
% 
% simplify(ft_lc)
% 
% pretty(ft_lc)

% % Ejercicio 3
% % 3.1. Determinar las cuatro funciones de transferencia que modelan el sistema.
% 
% s=tf('s');
% L=1e-6; R=2; Ka=0.042; J=10e-6; B=0.3e-5; Kb=0.042;
% 
% g1=1/((L*s)+R);
% 
% g2=1/((J*s)+B);
% 
% % % a) W_E con T=0
% ft_la1=g1*Ka*g2;
% 
% ft_lc1=feedback(ft_la1, Kb)
% 
% 
% % b) W_T
% ft_la2=g2;
% ft_lc2=feedback(ft_la2,(Ka*Kb*g1));
% ft_lc2=ft_lc2*-1
% 
% % c) I_E
% ft_la3=g1;
% ft_lc3=feedback(ft_la3, (g2*Ka*Kb))
% 
% % d) I_T
% 
% ft_la4=g2*(-Kb)*g1;
% ft_lc4=feedback(ft_la4,-Ka)
% 
% % 3.2. Simular la respuesta del sistema para e(t ) = 24u(t ) y t L (t ) = 0.01u (t − 2) .
% 
% % % Simular la respuesta para e(t) = 24u(t) y tL(t) = 0.01u(t-2)
% % t = 0:0.01:2.5;
% % e_t = 24 * ones(size(t));  % Entrada escalón de 24V
% % tl_t = 0.01 * (t >= 2);    % Perturbación escalón a t=2s
% % 
% % % Respuesta de la velocidad W(s) y la corriente I(s)
% % % Puede verse que la respuesta real para cada salida corresponde a la suma
% % % de las entradas.
% % 
% % W_response = lsim(ft_lc1, e_t, t) + lsim(ft_lc2, tl_t, t);
% % I_response = lsim(ft_lc3, e_t, t) + lsim(ft_lc4, tl_t, t);
% % 
% % % Graficar los resultados
% % figure;
% % subplot(2,1,1);
% % plot(t, W_response,'LineWidth', 1.5);
% % title('Respuesta de la Velocidad W(t)');
% % xlabel('Tiempo [s]');
% % ylabel('Velocidad [rad/s]');
% % 
% % subplot(2,1,2);
% % plot(t, I_response,'LineWidth', 1.5);
% % title('Respuesta de la Corriente I(t)');
% % xlabel('Tiempo [s]');
% % ylabel('Corriente [A]');
% 
% %Cada plot representa correctamente los fenomenos, al arrancar el motor se
% %produce un pico de corriente debido a los bobinados y a los dos segundos
% %se produce una perturbacion que podria ser una carga, disminuyendo la
% %velocidad del rotor.
% 
% 
% %Ejercicio 4
% % 4.1. Determinar las funciones de transferencia W ( s ) / Wr ( s ) y W ( s ) / TL ( s ) .
% 
% Kp=20;
% W_Wr=feedback((Kp*ft_lc1),1);
% 
% W_TL=-feedback(g2,(Ka*(Kp+Kb)*g1));
% 
% % 4.2. Simular la respuesta del sistema para Wr (t ) = 300u (t ) y t L (t ) = 0.01u (t − 2) .
% 
% % Simular la respuesta para wref(t) = 300u(t) y tL(t) = 0.01u(t-2)
% t = 0:0.01:2.5;
% e_t = 300 * ones(size(t));  % Entrada escalón de 300V
% tl_t = 0.01 * (t >= 2);    % Perturbación escalón a t=2s
% 
% % Respuesta de la velocidad W(s) y la corriente I(s)
% % Puede verse que la respuesta real para cada salida corresponde a la suma
% % de las entradas.
% 
% W_response = lsim(W_Wr, e_t, t) + lsim(W_TL, tl_t, t);
% 
% 
% % Graficar los resultados
% figure;
% 
% plot(t, W_response,'LineWidth', 1.5);
% title('Respuesta de la Velocidad W(t)');
% xlabel('Tiempo [s]');
% ylabel('Velocidad [rad/s]');

%%Puede verse que ajustar la variable Kp=20 nos acercamos a los 300r/seg
%%que tenemos de referencia 

