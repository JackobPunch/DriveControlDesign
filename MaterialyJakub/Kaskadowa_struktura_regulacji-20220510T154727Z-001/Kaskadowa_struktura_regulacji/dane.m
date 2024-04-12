%% Dane projektu 
% Plik, skrypt ten w założeniu ma się atomatycznie wczytywać po uruchomieniu
% projektu

% Ważne: ze wzgledu na kolizje oznaczeń zakomentować to co nie "idzie" do
% simulinka

%% Dane silnika
Pn=17; % [kW] !!
nn=1500; % [obr/min]
Un=220; % [V]
In=88; % [A]
Js=0.275; % [kgm2]
Rt=0.15; % [Ohm]
Lt=0.01875; % [H]
omega_n=nn*2*pi/60;
psi_en=(Un-Rt*In)/omega_n; % [Wb]
Mn=psi_en*In; % [Nm]
% przyjete dane
lambda=2;
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

%% Regulatory 
% obliczenia potrzebne do simulinka
statyzm=0.02*omega_n
regulatory % uruchomienie skryptu regulatory.m

%% Inne
u_z0 = lambda*In/kz*1.015
