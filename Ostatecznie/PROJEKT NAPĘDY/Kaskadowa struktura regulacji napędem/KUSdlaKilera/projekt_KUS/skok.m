%% Skrypt do schematu labolatoriÛw z KUS
% Autor: Tomasz Kielar, Tomasz Telesz

clear all; close all; clc;
%% dane programu

Pn = 54000;  %moc znamionowa [W]
Un = 440;   %napiÍcie znamionowe [V]
In = 130.7;    %prπd znamionowy [A]
n = 1450;  %prÍdkoúÊ znamionowa [obr/min]
Rt = 0.254;  %rezystancja twornika [Ohm]
Lt = 1.63*10^-3;   %indukcujnosc twornika [H]
Js = 0.97;  %moment bezwladnosci [kg*m^2]
Lambda = 2;
p = 50;
tau0 = 3.3e-3;
Tend = 6;


%Fp = 
%% Wzory

%Jmr
Jmr = 5*Js;

%J
J = Js+Jmr;

%prÍdkoúÊ kπtowa
wn = 2*pi*n/60;

%psi
psi = ((Un-Rt*In)*30)/(pi*n);

%B
B= J*Rt/(psi^2);

%T
T = Lt/Rt;

%Mn moment obciπøenia
Mn = psi * In; 

%T1
T1 = (B/2)*(1-sqrt(1-4*T/B));

%B1
B1 = B - T1;

%Beta
Beta = Lambda/p;

%Y
Y = 10/(2.5*In);

%kz
kz = (B1-Beta)/(Y*B1);

%Kt
Kt = 10/(1.2*wn);

%deltawn
deltawn = wn*0.04;

%Kp
Kp = 0.15*Un;

%Uz0
Uz0 = Lambda*In*(Y*B1)/(B1-Beta);

TR = 4*Beta;

%% Badanie odpowiedzi uk≥adu na skok jednostkowy

sim Simulinkasdasd.slx
%% Obliczanie ograniczeÒ

% ograniczenie pradu
time = ones(1, length(Przebiegi.time));
Id = Lambda * In * time;

%ograniczenie pochodnej pradu
derI = p * In* time;

%ograniczenie predkosci obrotowej
wLimit = (n*2*pi/60)*time;

%% ODPOWIEDè NA SKOK JEDNOSTKOWY
f = figure; set(f,'name','Odpowiedz skokowa predkosci katowej','numbertitle','off', 'Position', [500, 200, 700, 500]);

%Predkosc obrotowa
% subplot(3,1,1); 
plot(Przebiegi.time, Przebiegi.signals.values(:,1));
title('Odpowiedü skokowa prÍdkoúci kπtowej'); xlabel('Czas [s]'); ylabel('prÍdkoúÊ kπtowa [rad/s]'); hold on; grid on
plot(Przebiegi.time, wLimit,'-r')
legend('prÍdkoúÊ kπtowa','ograniczenie prÍdkoúci kπtowej')

% % ograniczenie predkosci obrotowej
% subplot(3,2,1); plot(Przebiegi.time,wLimit, '-r');
% hleg = legend('PredkoúÊ obrotowa','Ograniczenie','Location','NorthEast');
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

% subplot(3,1,2); 
f = figure; set(f,'name','Odpowiedz, pochodna pradu','numbertitle','off');
plot(Przebiegi.time, Przebiegi.signals.values(:,4));
title('Pochodna prπdu twornika'); xlabel('Czas [s]'); ylabel('WartoúÊ [A]'); grid on
hold on;
plot(Przebiegi.time, derI,'r')
legend('pochodna prπdu ','ograniczenie wartoúci pochodnej')

% subplot(3,2,3); plot(Przebiegi.time, Przebiegi.signals.values(:,3));
% title('GIU'); xlabel('Czas [s]'); ylabel('Odpowiedü obiektu');
% 
% %Pochodna prπdu
% subplot(3,2,4); plot(Przebiegi.time, Przebiegi.signals.values(:,4));
% title('dGIU/dt'); xlabel('Czas [s]'); ylabel('Odpowiedü obiektu'); hold on;
% 
% %Ograniczenie pochodnej prπdu
% subplot(3,2,4); plot(Przebiegi.time, derI, '-r');
% hleg = legend('Pochodna prπdu','Ograniczenie','Location','NorthEast');
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%Prπd
% subplot(3,1,3); 
f = figure; set(f,'name','Odpowiedz skokowa pradu twornika','numbertitle','off');
plot(Przebiegi.time, Przebiegi.signals.values(:,3));
title('Odpowiedz skokowa prπdu twornika'); xlabel('Czas [s]'); ylabel('WartoúÊ znamionowa [A]'); hold on; grid on
plot(Przebiegi.time, Id,'-r');
legend('prπd twornika','ograniczenie prπdu twornika')

%Ograniczenie prπdu
% subplot(3,1,3); 
% plot(Przebiegi.time, Id, '-r');
% hleg = legend('Prπd','Ograniczenie','Location','NorthWest');
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]);
% % axis([0 1.5 0 450]);
% hold off;