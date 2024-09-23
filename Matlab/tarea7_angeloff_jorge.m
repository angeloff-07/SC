%% Estabilidad Absoluta
% 11. Estudiar la estabilidad de los sistemas cuyas ecuaciones características son:
close all; clear all; clc


N=[1 2.996 3 10.998];
p=roots(N)

%Inestable tiene polos con parte real positiva
%% 

close all; clear all; clc


% syms s
% N= collect((2*( s + 1)^2),s);

N=[2 4 2]
p=roots(N)

%Estable, tiene un polo doble en -1

%% 

close all; clear all; clc
syms s k

% N=3*s^2 + (2 + k )*s + 1; 
N=[3 (2+k) 1];
p=roots(N);

%Su estabilidad depende del valor que tome k

%% 12. Para los sistemas que se muestran a continuación:
close all; clear all; clc

s=tf('s');

G1=(50*s +8)/(s^3 + 11*s^2 + 23*s -8);
rlocus(G1)

%Inestable tiene un polo en 0.303

step(G1)

%% 

close all; clear all; clc

s=tf('s');

G2=(50*s +8)/(s^3 + 11*s^2 + 23*s +8);
rlocus(G2)

%Estable, no tiene polos con parte real positiva
step(G2)

%%
 close all; clear all; clc

s=tf('s');

G3=(50*s - 8)/(s^3 + 11*s^2 + 23*s +8);
rlocus(G3)

%Estable, no tiene polos con parte real positiva
 step(G3)

 %%
 close all; clear all; clc

s=tf('s');

G4=(50*s + 8)/(s^3 + 23*s +8);
rlocus(G4)

%Inestable tiene polos complejos conjugados en 0.173+-4.81i
 step(G4)

 %%
  close all; clear all; clc

s=tf('s');

G8=(s)/(2*s^5 + 4*s^4 + 2*s^3 + 4*s^2 + s +6);
rlocus(G8)

%Inestable tiene polos complejos conjugados en 0.579 +- 0.847i
 step(G8)

 %%
  close all; clear all; clc

s=tf('s');

G9=(54)/(s^3 + 13*s^2 + 55*s +75);
rlocus(G9)

%Inestable para ganacias mayores a 11.9

  step(G9)

  %%  Estabilidad Relativa
% 14. Determinar si los sistemas cuyas ecuaciones características se listan a continuación
% son estables. Verificar usando Octave

close all; clear all; clc


% syms s
% N= collect((2*( s + 1)^2),s);

N=[2 4 2]
p=roots(N)

%Estable, tiene un polo doble en -1

%%
close all; clear all; clc
syms k
N=[1 11 10 110]
p=roots(N)
%Estable

%%
close all; clear all; clc
syms k
N=[2 2 3 1 3 2 1]
p=roots(N)

%Inestable

%%

close all; clear all; clc
syms k
N=[1 3 -2 7 12]
p=roots(N)
%Inestable

%% 15. Determinar los límites de estabilidad del siguiente sistema:
close all; clear all; clc

syms k s 
G=1/(s+1)^3;
FdTLC=k*G/(1+k*G);
FdTLC=collect(simplify(FdTLC),'s')

G=zpk([],[-1 -1 -1],1)
step(feedback(7*G,1),feedback(8*G,1),feedback(9*G,1),100);grid
legend('k=7','k=8','k=9','location','best')

%Puede verse que con k=9 inestable y diverge, k=8 inestable y oscila, k=7
%estable y converge
%% 


