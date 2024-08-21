% Unidad N° 3 
%% Ejercicio 1
close all; clear all; clc
% 
% s=tf('s');
% % 
% Ys_Xs=(5*(s+1))/((5*s^2)+(2*s)+3);
% 
% polos=pole(Ys_Xs);

% if real(polos)<0
%     disp('Sistema estable')
% 
% else
%     disp('Sistema inestable')
% end

% Ys_Xs=roots([5 2 3]);
% step(Ys_Xs,50)

%% Ejercicio 2
% Aplicando TVF siendo E(s)=10V/s

% w_s=(3000*10)/1837.5;
%%
syms zeta omega_n

% Ecuaciones
eq1 = 2*zeta*omega_n == 156.25;
eq2 = omega_n^2 == 1837.5;

% Resolver las ecuaciones
sol = solve([eq1, eq2], [zeta, omega_n]);


zeta = 1.8225

omega_n = 42.8661

K = 3000 / omega_n ^ 2 

ts = 4 / (zeta * omega_n)
%%
s=tf('s');
W_V=3000/(s^2+156.25*s+1837.5);
step(W_V,1)

% pzmap(W_V);
