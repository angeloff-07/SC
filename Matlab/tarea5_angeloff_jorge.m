%% 6. Esbozar la respuesta temporal al escalón unitario de los sistemas caracterizados por
% las siguientes funciones de transferencia.
% ej1
close all; clear all; clc
s=tf('s');

g1=15/(5*s + 1);
step(g1,100)

%% ej4
close all; clear all; clc
s=tf('s');
g1=15/(s + 5);
step(g1,20)

%% ej5
close all; clear all; clc
s=tf('s');

g5=15/(5*s+1);
g5.inputdelay=2;

step(g5,50)

%% ej6
close all; clear all; clc
s=tf('s');
g6=625/(s^2+60*s+625);
step(g6, 10);

%% ej10
close all; clear all; clc
s=tf('s');
g10=187500/(s^2+600*s+62500); 
step(g10, 0.1);
 %% ej11
close all; clear all; clc
s=tf('s');
g11=0.1875/(s^2+0.6*s+0.0625);
step(g11,100);

%% ej13
close all; clear all; clc
s=tf('s');
g13=0.0625/(s^2+0.3*s+0.0625);
step(g13, 100);

%Sobrepasamiento (OverShoot=9.46%)
os=9.46/100;

k=1;
ymax=1.09;
tp=16;

mp=(ymax-k)/k;

psita=sqrt(log(mp)^2 / (pi^2 + log(mp)^2));

wn=pi / (tp*sqrt(1-psita^2));

%% ej16
close all; clear all; clc
s=tf('s');

g16=0.0625/(s^2+0.0625);
step(g16,500);

%% ej17
close all; clear all; clc
s=tf('s');
g17=0.25/(s^2+0.25);
step(g17,500);

%% ej19
close all; clear all; clc
s=tf('s');
g19=2/(3*s+1)^2;step(g19,50);