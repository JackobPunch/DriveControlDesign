%% Dane projektu 
%% Dane silnika
Pn=26; % [kW] !!
nn=625; % [obr/min]
Un=230; % [V]
In=132; % [A]
Js=1.25; % [kgm2]
Rt=0.202; % [Ohm]
Lt=1.9; % [H]
%% obliczone parametry 
omega_n=nn*2*pi/60;
psi_en=(Un-Rt*In)/omega_n; % [Wb]
Mn=psi_en*In; % [Nm]
T=Lt/Rt;
J=22*Js;
B=J*Rt/psi_en^2

% przyjete dane
lambda=1.8;
p=50;
% Ograniczenia
Imax = lambda*In;
dIdtmax = p*In;
omega_max = omega_n;


% założone dane do regulatora
Y=10/(2.5*In);
Kt=10/(1.2*omega_n);
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
xlabel('t [s]');ylabel('I[A]');
legend('pochodna prądu twornika','Ograniczenie pochodnej prądu twornika')
grid on
figure
plot(t,odpw,t,omega_n.*ones(1,length(t)));
title('odpowiedź skokowa prędkości obrotowej');
xlabel('t [s]');ylabel('\omega [rad/s]');
legend('prędkość obrotowa','znamionowa prędkość obrotowa')
grid on