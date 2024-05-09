clear; close all; clc;
%dane z tabelki
PN=54e3;
UN=440;
IN=130.7;
nN=1450; wN=nN*pi/30;
Rt=0.254;
Lt=1.63;
Js=0.97;

%obliczenia tych ograniczeń i innych
lambdaN=1.8;
p=50;
Y=10/(2.5*IN);
KT=10/(1.2*wN);
Kp=1.5*UN/10;

psie=(UN-IN*Rt)/(wN);

T=0.001*Lt/Rt;

%sprwadzamy czy B>4T - przyjęliśmy J=Js+0*Js bo już jest od razu większe 4
%razy

T4= 4*T %4*T

B=(Js+0*Js)*Rt/(psie^2)

%transmitancje - tak żeby nie było, że nie mamy

%definicujemy jednostkę Laplace'a

s=tf('s');

%GwU(s)

GwU=tf([1/psie],[B*T, B, 1])

%GIU

GIU=tf([B/Rt, 0], [B*T, B, 1])

%GdIU - transmitancja z pochodną prądu

GdIU=GIU*s

sim("Ne_proj_zad1.slx")


%============wykresy z Simulinka=============

%prędkość
figure
plot(GwU.time, GwU.signals.values)
title('Odpowiedz skokowa transmitancji G_{wU}')
xlabel('t[s]')
ylabel('w(t) [rad/s]')

%prąd
figure
plot(GIU.time, GIU.signals.values)
hold on
plot(GIU.time,ones(1,length(GIU.time)).*IN*lambdaN,'g-')
hold on
title('Odpowiedz skokowa transmitancji G_{IU}')
hold on
xlabel('t[s]')
ylabel('I_t(t) [A]')


%pochodna prądu
figure
plot(GdIU.time, GdIU.signals.values)
hold on
plot(GdIU.time,ones(1,length(GdIU.time)).*IN*p,'g-')
hold on
title('Odpowiedz skokowa transmitancji G_{dIU}')
hold on
xlabel('t[s]')
ylabel('dI_t/dt [A]')




