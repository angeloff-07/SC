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

% G=zpk([],[-1 -1 -1],1)
% step(feedback(7*G,1),feedback(8*G,1),feedback(9*G,1),100);grid
% legend('k=7','k=8','k=9','location','best')

%Puede verse que con k=9 inestable y diverge, k=8 inestable y oscila, k=7
%estable y converge
%% 16. Determinar la estabilidad relativa de los siguientes sistemas. Validar los resultados obtenidos.
%Sistema 1

 close all; clear all; clc

 s=tf('s');

 g1=10/(s^2 + 10)

 rlocus(g1)

%Sistema estable para cualquier ganancia k

%% Sistema 2

 close all; clear all; clc

 s=tf('s');

 g2=(10*s + 20)/(s^2 + 120*s + 10)

 rlocus(g2)
%Sistema estable para cualquier ganancia k

%% Sistema 3

 close all; clear all; clc

 s=tf('s');

 g3=(45)/(s^3 + 12*s^2 + 10*s + 45)

 rlocus(g3)
%Sistema inestable para cualquier ganancia k>1.6666, calculado con criterio
%R-H y chequeado con LDR
step(feedback(1.60*g3,1),1000)

%% Sistema 4

 close all; clear all; clc

 s=tf('s');

 g4=((s + 10)*(s+ 20))/((s +1)*(s + 5))

 rlocus(g4)
%Sistema estable para cualquier ganancia k

%% Sistema 5

 close all; clear all; clc

 s=tf('s');

 g5=((s + 10)*(s+ 20))/((s - 1)*(s + 5))

 rlocus(g5)
%Sistema inestable para cualquier ganancia k<0.025, calculado con criterio
%R-H y chequeado con LDR
step(feedback(0.02*g5,1),10)

%% Sistema 6

 close all; clear all; clc

 s=tf('s');

 g6_1=(1)/(s + 10);
 g6_2=(s)/(s^2 + s + 1);
 g6=g6_1*g6_2
 rlocus(g6)
%Sistema estable para cualquier ganancia k

%% Ejemplo LDR, Técnica del lugar de raíces
close all; clear all; clc
G=tf([1 10],conv([1 5],[1 4 8]));
figure
% Lugar de raíces.
rlocus(G); sgrid(0.707,[4 8 12])
axis([-25 10 -25 25]);
% ASINTOTAS
h=line([0.5 0.5],[-25 25]) %interseccion de las asintotas y de donde hasta donde
set(h,'LineStyle','-.')
