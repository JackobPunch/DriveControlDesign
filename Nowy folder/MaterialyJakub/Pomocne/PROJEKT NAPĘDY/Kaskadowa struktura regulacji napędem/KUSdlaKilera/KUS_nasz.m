%% Skrypt do schematu labolatori�w z KUS
% Autor: Adam Miarka, Damian Kad�uczka
% data 10.03.2016 grupa 29

clear all; close all; clc;
%% dane programu

Js = 2.1;
Rt = 0.112;
Lt = 0.00159;
n = 585;
In = 215;
Un = 230;
Pn = 43000;
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

%pr�dko�� k�towa
wn = 2*pi*n/60;

%psi
psi = ((Un-Rt*In)*30)/(pi*n);

%B
B= J*Rt/(psi^2);

%T
T = Lt/Rt;

%Mn moment obci��enia
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

%% Badanie odpowiedzi uk�adu na skok jednostkowy

sim Simulinkasdasd.slx
%% Obliczanie ogranicze�

% ograniczenie pradu
time = ones(1, length(Przebiegi.time));
Id = Lambda * In * time;

%ograniczenie pochodnej pradu
derI = p * In* time;

%ograniczenie predkosci obrotowej
wLimit = n*time;

%% ODPOWIED� NA SKOK JEDNOSTKOWY
f = figure; set(f,'name','Odpowiedz skokowa predkosci katowej','numbertitle','off');

%Predkosc obrotowa
% subplot(3,1,1); 
plot(Przebiegi.time, Przebiegi.signals.values(:,1));
title('Pr�dko�� k�towa'); xlabel('Czas [s]'); ylabel('Warto�� [rad/s]'); hold on; grid on

% % ograniczenie predkosci obrotowej
% subplot(3,2,1); plot(Przebiegi.time,wLimit, '-r');
% hleg = legend('Predko�� obrotowa','Ograniczenie','Location','NorthEast');
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

% subplot(3,1,2); 
f = figure; set(f,'name','Odpowiedz, pochodna pradu','numbertitle','off');
plot(Przebiegi.time, Przebiegi.signals.values(:,4));
title('Pochodna pr�du'); xlabel('Czas [s]'); ylabel('Warto�� [A]'); grid on

% subplot(3,2,3); plot(Przebiegi.time, Przebiegi.signals.values(:,3));
% title('GIU'); xlabel('Czas [s]'); ylabel('Odpowied� obiektu');
% 
% %Pochodna pr�du
% subplot(3,2,4); plot(Przebiegi.time, Przebiegi.signals.values(:,4));
% title('dGIU/dt'); xlabel('Czas [s]'); ylabel('Odpowied� obiektu'); hold on;
% 
% %Ograniczenie pochodnej pr�du
% subplot(3,2,4); plot(Przebiegi.time, derI, '-r');
% hleg = legend('Pochodna pr�du','Ograniczenie','Location','NorthEast');
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%Pr�d
% subplot(3,1,3); 
f = figure; set(f,'name','Odpowiedz skokowa pradu twornika','numbertitle','off');
plot(Przebiegi.time, Przebiegi.signals.values(:,5));
title('Pr�d'); xlabel('Czas [s]'); ylabel('Warto�� znamionowa [A]'); hold on; grid on

%Ograniczenie pr�du
% subplot(3,1,3); 
plot(Przebiegi.time, Id, '-r');
hleg = legend('Pr�d','Ograniczenie','Location','NorthWest');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]);
% axis([0 1.5 0 450]);
hold off;


%% Obliczanie nastaw ci�g�ego regulatora pr�dko�ci

%Regulator typu P Kwp = Kri 
%Kw = Mn/psi*kz*Kt*deltawn;

%Regulator typu PI
Kw = J/(2*Kt*kz*Beta*psi);

%% Obliczanie nastaw ci�g�ego regulatora pr�du

m = T1;
V = (Beta*(Y*Kp*B))/((B1-Beta)*Rt);

%% Rozpocz�cie symulacji uk�ad�w z regulatorem ci�g�ym dla r�znych rodzaj�w obci��enia

%% ROZRUCH UK�ADU BEZ MOMENTU OBCIA�ENIA ***REGULACJA CI�G�A***
sim Model_simulink_no_Torque.slx;

f = figure; set(f,'name','No Torque','numbertitle','off');

subplot(2,2,1); plot(no_Torque.time,no_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(no_Torque.time,no_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(no_Torque.time,no_Torque.signals.values(:,3));
title('pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(no_Torque.time,no_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(no_Torque.time,no_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

% drukuj� powi�kszenie w celu zweryfikowania sta�ej czasowej uk�adu (beta)
f = figure; set(f,'name','Bez obci��enia - powi�kszenie','numbertitle','off');
plot(no_Torque.time,no_Torque.signals.values(:,4)); hold on; grid on
plot(no_Torque.time,2*In*0.632*ones(length(no_Torque.time),1),'--r');
plot(0.04*ones(100,1),linspace(0,400,100),'--r'); title('Pr�d silnika - weryfikacja sta�ej czasowej(Beta)');
xlabel('time[s]'); ylabel('I[A]'); grid on;
xlim([0 0.1]);

%% ROZRUCH ZE ZNAMIONOWYM MOMENTEM CZYNNYM ***REGULACJA CI�G�A***
sim Model_simulink_static_Torque.slx;

f = figure; set(f,'name','static_Torque','numbertitle','off');

subplot(2,2,1); plot(static_Torque.time,static_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(static_Torque.time,static_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(static_Torque.time,static_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(static_Torque.time,static_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(static_Torque.time,static_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off; 

% obserwujemy nawr�t silnika
f = figure; set(f,'name','Obserwujemy nawr�t silnika - moment czynny','numbertitle','off');
plot(static_Torque.time,static_Torque.signals.values(:,3)); grid on
title('Pr�dko�� obrotowa - nawr�t silnika dla obci��enia znamionowym momentem czynnym');
xlabel('Czas[s]'); ylabel('n[obr/min]'); grid on;
% axis([0 0.1 -0.1 0.3]);

%%  PRACA UK�ADU Z UDAROWYM MOMENTEM OBCI��ENIA ***REGULACJA CI�G�A***
sim Model_simulink_impact_Torque.slx;

f = figure; set(f,'name','impact_Torque','numbertitle','off'); title('TYYYTU�');

subplot(2,2,1); plot(impact_Torque.time,impact_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(impact_Torque.time,impact_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(impact_Torque.time,impact_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(impact_Torque.time,impact_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(impact_Torque.time,impact_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%%  PRACA UK�ADU Z BIERNYM MOMENTEM OBCI��ENIA ***REGULACJA CI�G�A***
sim Model_simulink_passive_Torque.slx;

f = figure; set(f,'name','passive_Torque','numbertitle','off'); title('TYYYTU�');

subplot(2,2,1); plot(passive_Torque.time,passive_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(passive_Torque.time,passive_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(passive_Torque.time,passive_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(passive_Torque.time,passive_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(passive_Torque.time,passive_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

% obserwujemy nawr�t silnika
f = figure; set(f,'name','Obserwujemy nawr�t silnika - moment bierny','numbertitle','off');
plot(passive_Torque.time,passive_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa - brak nawrotu silnika dla znamionowego momentu biernego');
xlabel('Czas[s]'); ylabel('n[obr/min]'); grid on;
% axis([0 0.1 -0.1 0.3]);

%%  PRACA UK�ADU Z PULSACYJNYM MOMENTEM OBCI��ENIA ***REGULACJA CI�G�A***
sim Model_simulink_pulse_Torque.slx;

f = figure; set(f,'name','pulse_Torque','numbertitle','off'); title('TYYYTU�');

subplot(2,2,1); plot(pulse_Torque.time,pulse_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(pulse_Torque.time,pulse_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(pulse_Torque.time,pulse_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(pulse_Torque.time,pulse_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(pulse_Torque.time,pulse_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%% Punkt 3 oblcizanie transmitancji zast�pczej uk�adu ci�g�ego

% reg_momentu*przekszta�tnik*silnik*bezw�adno��
RegPredkosc = tf([Kw*TR Kw],[TR,0]);        % transmitacja regulatora pr�dko�ci
RegMoment = tf([m 1],[V 0]);                % transmitancja regulatora momentu(pr�du)
PrzeksztTyryst = tf([Kp],[tau0 1]);         % transmitancja przekszta�tnika tyrystorowego
Silnik = tf([B 0],[Rt*B1*T1 Rt*(B1+T1) Rt]);% transmitancja silnik
Bezwlad = tf([1],[J 0]);                    % transmitancja momentu bezwladnosci
sys = series(RegMoment,PrzeksztTyryst);     % liczymy szeregowo transmitancje regulatora predkosci i przekszta�tnika tyrystorowego
sys = series(sys,Silnik);                   % liczymy szeregowo transmitancj� poprzedniego uk�adu i silnika
sys = feedback(sys, Y);                     % liczymy sprz�enie zwrotne uzyskanej transmitancji oraz Y
sys = series(RegPredkosc,sys);              % liczymy szeregowo transmitancj� poprzedniego uk�adu i regulatora pr�dko�ci
disp('Transmitancja otwartego uk�adu regulacji = ');
sys = series(sys,Bezwlad)                  % liczymy szeregowo transmitancj� poprzedniego uk�adu i momentu bezwladnosci
%% Charakterystyki bodego otwartego uk�adu
f = figure; set(f,'name','Bode plot dla uk�adu otwartego','numbertitle','off');
bode(sys);grid on;                          % charakterystyki Bodego

%% Charakterystyka nyquista otwartego uk�adu
f = figure; set(f,'name','Nyquist plot dla uk�adu otwartego','numbertitle','off');
subplot(2,1,1);
nyquist(sys);grid on;hold on                     % charakterystyki nyquista

tt = linspace(0,2*pi,200);
aa=0; bb=0; rr=1;
xx = aa+rr*cos(tt); yy = bb+rr*sin(tt);
plot(xx,yy,'-','LineWidth',2);
axis equal;

subplot(2,1,2);
nyquist(sys);

% axis([-0.025 0 -1e-3 1e-3]); 
%% liczymy transmitancj� zamkni�tego uk�adu w celu wyznaczenia wg
disp('Transmitancja zamkni�tego uk�adu regulacji = ');
sys = feedback(sys, Kt)

%% Charakterystyka Bodego uk�adu zamkni�tego
f = figure; set(f,'name','Bode plot dla uk�adu zamkni�tego','numbertitle','off');
bode(sys);grid on;                          % charakterystyki Bodego
[Gm,Pm,Wgm,Wpm] = margin(sys);

%% Charakterystyka nyquista uk�adu zamkni�tego oraz obliczanie zapasu fazy i wzmocnienia
f = figure; set(f,'name','Nyquist plot dla uk�adu zamkni�tego','numbertitle','off');
subplot(2,1,1);
nyquist(sys);grid on; title(sprintf('Zapas wzmocnienia %f[razy]\nZapas fazy %f[�]\nOp�nienie grupowe %f[s]',Gm,Pm,Pm*pi/180/Wpm));
subplot(2,1,2);
nyquist(sys); axis([-0.025 0 -1e-3 1e-3]); 
wg = 8.32;                                  %[rad/s] odczytane z wykresu
wgp = wg*100;
fgp = 2*pi*wgp;
Tp = 1/fgp;

N=20;
Tp = Beta/N;                                
%Ti = TR;                                     % czas zdwojenie do dyskretyzacji (4*Beta)
%% Liczymy nastawy dyskretnego regulatora obrot�w
K1 = Kw;
K2 = Kw*(Tp/TR-1);

%% Liczymy nastawy dyskretnego regulatora pr�du
K3 = m/V;
K4 = (Tp-m)/V;
%% Wykonujemy symulacj� dla �le dobranego(zbyt ma�ego) czasu pr�bkowania Tp i obserwujemy zniekszta�cenie przebieg�w
Tp = 0.01;
Q = N*2/2e15;
% ROZRUCH UK�ADU BEZ MOMENTU OBCI��ENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_no_torque.slx;

f = figure; set(f,'name','Discrete No Torque1','numbertitle','off');

subplot(2,2,1); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%% Liczymy b��d kwantyzacji dla dobrze dobranego Tp
%% Wykonujemy symulacj� dla zbyt du�ego kwantyzatora i obserwujemy rezultaty
Tp = Beta/N;
Q = 0.05;                   

% ROZRUCH UK�ADU BEZ MOMENTU OBCI��ENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_no_torque.slx;

f = figure; set(f,'name','Discrete No Torque2','numbertitle','off');

subplot(2,2,1); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%% Nast�pnie wykonujemy symulacj� dla ma�ego kwantyzatora i obserwujemy brak zak��ce�
Q = N*2/2e15;

% ROZRUCH UK�ADU BEZ MOMENTU OBCI��ENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_no_torque.slx;

f = figure; set(f,'name','Discrete No Torque3','numbertitle','off');

subplot(2,2,1); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;
%% Rozpocz�cie symulacji uk�ad�w z regulatorem dyskretnym dla r�nych rodzaj�w obcia�enia

%% ROZRUCH UK�ADU BEZ MOMENTU OBCI��ENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_no_torque.slx;

f = figure; set(f,'name','Discrete No Torque4','numbertitle','off');

subplot(2,2,1); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%% DZIA�ANIE UK�ADU Z UDAROWYM MOMENTEM OBCI��ENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_impact_torque.slx;

f = figure; set(f,'name','Discrete Impact Torque','numbertitle','off');

subplot(2,2,1); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on
%title('\frac{1}{2}','Interpreter','latex') 

%naniesienie ograniczenia
subplot(2,2,4); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%% ROZRUCH UK�ADU BEZ MOMENTU OBCI��ENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_no_torque.slx;

f = figure; set(f,'name','Discrete No Torque5','numbertitle','off');

subplot(2,2,1); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%% DZIA�ANIE UK�ADU Z PULSACYJNYM MOMENTEM OBCI��ENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_pulse_torque.slx;

f = figure; set(f,'name','Discrete Pulse Torque','numbertitle','off');

subplot(2,2,1); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on
%title('\frac{1}{2}','Interpreter','latex') 

%naniesienie ograniczenia
subplot(2,2,4); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;

%% ROZRUCH UK�ADU ZE ZNAMIONOWYM MOMENTEM OBCI��ENIA ***REGULACHA DYSKRETNA***
deltaIn=(psi*V*Mn)/(psi^2*V+J*Kp*Y);
Uz0 = (Lambda*In+deltaIn)*(Y*B1)/(B1-Beta);

sim Model_simulink_dyskretny_static_torque.slx;

f = figure; set(f,'name','Discrete Static Torque','numbertitle','off');

subplot(2,2,1); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,1));
title('Napi�cie silnika'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,2); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,2));
title('Napi�cie sterownika pr�du'); xlabel('Czas [s]'); ylabel('U[V]'); grid on

subplot(2,2,3); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,3));
title('Pr�dko�� obrotowa'); xlabel('Czas [s]'); ylabel('n[obr/min]'); grid on

subplot(2,2,4); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,4));
title('Pr�d silnika'); xlabel('Czas [s]'); ylabel('I[A]'); hold on; grid on

%naniesienie ograniczenia
subplot(2,2,4); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,6),'-r');
hleg = legend('Warto�� pr�du','Ograniczenie','Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1]); hold off;



