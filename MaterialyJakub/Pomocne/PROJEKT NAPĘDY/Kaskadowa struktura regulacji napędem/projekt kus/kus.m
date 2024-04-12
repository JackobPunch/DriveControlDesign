%% KUS
%%

clear all; close all; clc

%% PARAMETRY ZNAMIONOWE
Tp = 0.1;
Pn = 54000;
Un = 440;
In = 130.7;
nn = 1450;
wn = nn*pi/30;
Rt = 0.254;
Lt = 0.00163;
Js = 0.97;
lambdan = 2;
p = 50;



%% PARAMETRY OBLICZONE

T = Lt / Rt;
psien = (Un-Rt*In)/wn;
J = 8 * Js;
B = J * Rt/(psien.^2);

% Wyznaczenie nastaw regulatorów ci¹g³ych
Y = 10/(2.5*In);
Kt = 10/(1.5*wn);
Kp = (1.5*Un)/10;

% Regulator pr¹du
beta = lambdan / p;
m = 0.5 * B * (1-sqrt(1-4*T/B));
kz = (B-m-beta)/(Y*(B-m));
V = (beta*Y*Kp*B)/(((B-m)-beta)*Rt);
uz0 = lambdan*In*Y*(B-m)/(B-m-beta);

% Regulator prêdkoœci
TR = 4 * beta;
dwm = 0.02 * wn;
kw = Js/(psien*kz*Kt*dwm);
Kw = J/(2*Kt*kz*beta*psien);

%%%
Id = lambdan * In;
dIt = p * In;

% Transmitancje uk³adu silnika pr¹du sta³ego

opt = stepDataOptions('InputOffset',0,'StepAmplitude',Un);

Gwu = [tf([1/psien],[B*T B 1])]
Gwm = [tf([(Rt/psien^2)*T Rt/psien^2],[B*T B 1])]
Giu = [tf([B/Rt 0],[B*T B 1])]
Gim = [tf([1/psien],[B*T B 1])]
Gdiu = [tf([B/Rt 0 0], [B*T B 1])]

[y1,t1]=step(Gwu,opt);
WN = wn*ones(length(t1),1);
figure(1)
plot(t1,y1)
hold on
plot(t1,WN,'r--')

[y2,t2]=step(Giu,opt);
id = Id*ones(length(t2),1);
title('OdpowiedŸ skokowa prêdkoœci k¹towej dla napiêcia U_N');
xlabel('t [s]'); ylabel('\omega [rad/s]');
grid on;

figure(2)
plot(t2,y2); hold on
plot(t2,id,'r--')

[y3,t3]=step(Gdiu,opt);
dit = dIt*ones(length(t3),1);
title('OdpowiedŸ skokowa pr¹du twornika dla napiêcia U_N');
xlabel('t [s]'); ylabel('I [A]');
grid on;
figure(3)
plot(t3,y3); hold on
plot(t3,dit,'--r')
title('Pochodna pr¹du twornika dla napiêcia U_N')
xlabel('t [s]'); ylabel('dI/dt [A/s]')
grid on;