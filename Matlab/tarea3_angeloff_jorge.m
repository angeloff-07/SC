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

%% Ejercicio 4
close all; clear all; clc

syms  XA XS W s

MA = 500; % Masa del automóvil
MS = 50;  % Masa de la suspensión
BA = 80;  % Amortiguador
KR = 10;  % Resorte entre automóvil y suspensión
KC = 50;  % Resorte de la cubierta

eq1 = (MA * s^2) * XA == BA * s * (XS - XA) + KR * (XS - XA);
eq2 = (MS * s^2) * XS == BA * s * (XA - XS) + KR * (XA - XS) + KC * (W - XS);

% Resolver el sistema para XA y XS
[XA_sol, XS_sol] = solve([eq1, eq2], [XA, XS]);

% Calcular la función de transferencia Xa/W
Xa_W = simplify(XA_sol / W)
pretty(Xa_W)
%% 
close all; clear all; clc
s=tf('s');

num = [ 8 , 1 ];
den = [ 50, 88, 61, 8, 1];

Xa_W=tf(num, den)

step(Xa_W, 200);grid
title('Respuesta al escalón del sistema de amortiguación')
xlabel('Tiempo [s]'); ylabel('Amplitud [m]')

%% Ejercicio de clase teorica
close all; clear all; clc

syms s M1 M2 B1 B2 B3 K F Y1 Y2

eq1 = ( M2 * (s^2) * Y2 ) + ( B3 * s *Y2) + ( B2 * s * ( Y2 - Y1 ) ) + ( K * ( Y2 - Y1)) == 0;
eq2 = ( M1 * (s^2) * Y1 ) + ( B1 * s *Y1) + ( B2 * s * ( Y2 - Y1 ) ) + ( K * ( Y2 - Y1)) == F;

[Y1_SOL, Y2_SOL] =solve ([eq1, eq2], [Y1, Y2]);

Y1_F = simplify(Y1_SOL/F);
pretty(Y1_F)

Y2_F = simplify(Y2_SOL/F);
pretty(Y2_F)


%% ejercicio 5
close all; clear all; clc

syms v ra i la kb wm jm wm ki bm s

%suponiendo tl=0 condicion incial nula

eq1 = v == ra * i + la * s * i +kb * wm;
eq2 = jm * s *wm == ki * i - bm * wm;

[w_sol, i_sol] = solve([eq1, eq2],[wm, (i)]);

w_v = simplify(w_sol/v);

pretty(w_v)

theta_v = w_v * s;

%% 
close all; clear all; clc

syms v ra i la kb wm jm wm ki bm s tl

%suponiendo tl=0 condicion incial nula

eq1 = 0 == ra * i + la * s * i +kb * wm;
eq2 = jm * s *wm == ki * i - bm * wm - tl;

[w_sol, i_sol] = solve([eq1, eq2],[wm, (i)]);

w_tl = simplify(w_sol/tl);

disp('w_tl = ')

pretty(w_tl)

%% 

close all; clear all; clc

s = tf('s');

ra = 5.8; 
la = 135.e-6;
kb = 14.48e-3; 
jm = 1.6e-2;
ki = 14.48e-3;
bm = 10e-5;


G = ki / (( bm * ra ) + ( kb * ki ) + ( jm * la * (s^2) ) + (bm * la * s) + (jm * ra * s))

step(12*G, 1000);

%% Ejercicio 7

close all; clear all; clc

s = tf('s');

ra = 5.8; 
la = 135.e-6;
kb = 14.48e-3; 
jm = 1.6e-2;
ki = 14.48e-3;
bm = 10e-5;


tm = ( ra * jm ) / ((ra * bm) + ( ki * kb ));

G1 = 1 / ((tm * s) + 1 )
G2 = 1 / (((tm/2)*s )+1) %% llega mas rapido al valor de regimen

step(G1, G2)

%% Ejercicio 6

close all; clear all; clc

syms s m1 x1 b1 u k1 k2 x2 m2 

eq1 = m1 * s^2 *x1 == b1 * ( s * u - s * x1 ) + k1 * ( u - x1 ) +k2 * ( x2 - x1);
eq2 = m2 * s^2 * x2 == k2 * ( x1 - x2 );

[x1_sol, x2_sol] = solve ( [ eq1, eq2 ],[ x1 , x2 ]);

x1_u = simplify(x1_sol/u);
x2_u = simplify(x2_sol/u);

x1_u = collect(x1_u, s); %%para ordenar las potencias
x2_u = collect(x2_u, s);

disp('x1_u = ')
pretty(x1_u)

disp('x2_u = ')
pretty(x2_u)

%% Ejercicio 8
close all; clear all; clc

syms m M i s theta l g x f b

eq1 = i * s^2 * theta + m * l^2 * s^2 * theta == m * g * l * theta + m * l * s^2 * x;
eq2 = (M + m) * s^2 * x == f - b * s *x + m * l * s^2 * theta;

[theta_sol, x_sol] = solve ([eq1, eq2],[theta, x]);

theta_f = simplify(theta_sol / f);

theta_f = collect(theta_f, s);

disp('theta_f = ')
pretty(theta_f)

%% Ejercicio 9
close all; clear all; clc

syms m s theta c g tm

syms r i v l kb theta_m j ka b 

eq1 = m * s^2 * theta + c * s * theta  + g * theta == tm;

[theta_sol] = solve([eq1], [theta]);
theta_sol = simplify(theta_sol)

theta_tm = simplify(theta_sol / tm);
theta_tm = collect(theta_tm, s);
disp('theta_tm = ')
pretty(theta_tm)

tm = ka * i;

eq2 = v == r * i +l * s * i + kb * s * theta_m;

eq3 = j * s^2 * theta_m == ka * i - b * s * theta_m;

[i_sol, theta_m_sol] = solve([eq2, eq3],[i, theta_m]);

tm_v = simplify((ka * i_sol) / v);
tm_v = collect(tm_v, s);
disp('tm_v = ')
pretty(tm_v)

theta_v = theta_tm * tm_v;
theta_v = collect(theta_v, s);
disp('theta_v = ')
pretty(theta_v)

%%
close all; clear all; clc

syms e s i1 i2 r1 r2 l1 l2 c

eq1 = e == i1 * ( r1 + s * l1 + (1/(s*c))) + i2 * (-1/(s*c));
eq2 = 0 == i1 * ( -1/(s*c) ) + i2 * ( r2 + s * l2 + 1/(s*c));

[i1_sol, i2_sol] = solve([eq1, eq2],[i1, i2]);

% i1_sol=simplify(i1_sol);
% i1_sol=collect(i1_sol,s);
% disp('i1_sol =')
% pretty(i1_sol)
% 
% i2_sol=simplify(i2_sol);
% i2_sol=collect(i2_sol,s);
% disp('i2_sol =')
% pretty(i2_sol)

er_e=(i2_sol*r2)/e;

er_e=simplify(er_e);
er_e=collect(er_e,s);
disp('er_e =')
pretty(er_e)
%%
close all; clear all; clc

S = tf('s');
R1 = 100;
R2 = 250;
C = 1e-6;
L1 = 100e-3;
L2 = 100e-3;

G= R2/((S^3)*C*L1*L2 + (S^2)*((C*L1*R2)+(C*L2*R1)) + S*(L1 + L2 +C*R1*R2) + (R1+R2));
G

step(G, 10)
%%