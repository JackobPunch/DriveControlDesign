clc
clear all
close all

T=20;
c = 0;
b = 1;
yref = 0.5;
Kr = 2;
Ti = 2;
Td = 0.2;
Tdf = 0.1*Td;
Tt = 0.5;
ogr = 1;

Tmin = 0.25;
N = 5;
Tp = (Tmin)/N;

%sim('pid.mdl')
%n=1496;
%dt=100/1496;
%t=0:dt:100-dt;
%figure;
%hold on, grid on
%plot(t,ciagly.signals.values)
%legend('anti-windup','windup');xlabel('[s]');
%hold off
sim('gotowe.mdl')

figure;
subplot(211)
hold on, grid on
plot(tout,out.signals.values)
title('Przed zmiana');
legend('Ciagly anti-windup','Ciagly windup','Dyskretny ZOH','Dyskretny Tustin');
xlabel('[s]');

c = 1;
b = 1;
yref = 0.5;
Kr = 2;
Ti = 2;
Td = 0.2;
Tdf = 0.1*Td;
Tt = 0.5;
ogr = 1;

Tmin = 0.25;
N =10;
Tp = (Tmin)/N;


sim('gotowe.mdl')

subplot(212)
hold on, grid on
plot(tout,out.signals.values)
title('Po zmianie');
legend('Ciagly anti-windup','Ciagly windup','Dyskretny ZOH','Dyskretny Tustin');
xlabel('[s]');