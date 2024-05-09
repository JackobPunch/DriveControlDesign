close all, clear all, clc
%% Dane projektu 
%% Dane silnika
Pn=17; % [kW] !!
nn=700; % [obr/min]
Un=230; % [V]
In=85; % [A]
Js=0.75; % [kgm2]
Rt=0.253; % [Ohm]
Lt=1.9*10^(-3); % [H]

%% obliczone parametry 
w_n=nn*2*pi/60;
psi_en=(Un-Rt*In)/w_n; % [Wb]
Mn=psi_en*In; % [Nm]
T=Lt/Rt;
J=22*Js;
B=J*Rt/psi_en^2;

% przyjete dane
lambda=1.8;
p=50;
% Ograniczenia
Imax = lambda*In;
dIdtmax = p*In;
omega_max = w_n;


% założone dane do regulatora
Y=10/(2.5*In);
Kt=10/(1.2*w_n);
Kp=0.15*Un;
% _________________
J=11*Js; %[kgm2]
% _________________
tau0=3.3*10^-3;


%% Model silnika (transmitancje)
% model_silnika

%% Model przekształtnika tyrystorowego
% ważne: allmargin()

% %% Regulatory 
% % obliczenia potrzebne do simulinka
% statyzm=0.02*omega_n
% regulatory % uruchomienie skryptu regulatory.m

%% Inne
% u_z0 = lambda*In/kz*1.015
%% transmitancje modeli
GwU=tf(1/psi_en,[B*T B 1])
GwM=tf((Rt/psi_en^2)*[T 1],[B*T B 1])
GIU=tf([B/Rt 0],[B*T B 1])
GIUdt=tf([B/Rt 0 0],[B*T B 1])
GIM=tf(1/psi_en,[B*T B 1])
%% odpowiedzi skokowe 
t=linspace(0,2,1000);
odpI=step(GIU,t);odpI=Un*odpI;
odpdI=step(GIUdt,t); odpdI=Un*odpdI;
odpw=step(GwU,t); odpw=Un*odpw;
figure
plot(t,odpI,t,Imax.*ones(1,length(t)));
title('odpowiedź skokowa prądu twornika');
xlabel('t [s]');ylabel('I[A]');
legend('prąd twornika','Ograniczenie prądu twornika')
 grid on
figure
plot(t,odpdI,t,dIdtmax.*ones(1,length(t)));
title('odpowiedź skokowa pochodnej prądu twornika');
xlim([0 0.2])
xlabel('t [s]');ylabel('I[A]');
legend('pochodna prądu twornika','Ograniczenie pochodnej prądu twornika')
grid on
figure
plot(t,odpw,t,w_n.*ones(1,length(t)));
title('odpowiedź skokowa prędkości obrotowej');
xlabel('t [s]');ylabel('\omega [rad/s]');
legend('prędkość obrotowa','znamionowa prędkość obrotowa')
grid on



% Do obczajenia czy sie przyda

%PUNKT 6
% 
% %nastawy regulatora prądu
% 
% T1=0.5*B*(1-sqrt(1-(4*T/B)));
% beta=lambdan/p;
% B1=B-T1;
% Kz=(B1-beta)/(Y*B1);
% m=T1;
% V=((Kp*Y)/R1)*((B*beta)/(B1-beta));
% uz0=lambdan*In*((Y*B1)/(B1-beta))
% 
% %regulator predkosci P:
% 
% KwP=Mn/(psie*Kz*KT*6.078)
% 
% %regulator predkosci PI:
% 
% KwPI=J/(2*beta*psie*KT*Kz)
% TrPI=4*beta
