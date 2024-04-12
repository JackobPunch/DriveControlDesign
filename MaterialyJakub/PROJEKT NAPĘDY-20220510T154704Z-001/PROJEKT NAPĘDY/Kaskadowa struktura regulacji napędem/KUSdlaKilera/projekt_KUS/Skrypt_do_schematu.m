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
tau0 = 3.3e-4;
Tend = 20;
%Fp = 
%% Wzory

%Jmr
Jmr = 5*Js;

%J
J = Js+Jmr;

%prÍdkoúÊ kπtowa
wn = (2*pi*n)/60;

%psi
psi = (Un-(Rt*In))/wn;

%B
B= (J*Rt)/(psi^2);

%T
T = Lt/Rt;

%Mn moment obciπøenia
Mn = psi * In; 

%T1
T1 = (B/2)*(1-sqrt(1-(4*T/B)));

%B1
B1 = B - T1;

%Beta
Beta = Lambda/p;

%Y
Y = 7.8/(2.5*In);

%kz
kz = (B1-Beta)/(Y*B1);

%Kt
Kt = 10/(1.2*wn);

%deltawn
deltawn = wn*0.05;

%Kp
Kp = 0.15*Un;

%Uz0
Uz0 = (Lambda*In)*((Y*B1)/(B1-Beta));

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
wLimit = n*time;

%% ODPOWIEDè NA SKOK JEDNOSTKOWY
f = figure; set(f,'name','Odpowiedzi obiektu na skok jendostkowy','numbertitle','off', 'Position', [500, 100, 700, 600],'Color','white');

%Predkosc obrotowa
subplot(311); plot(Przebiegi.time, Przebiegi.signals.values(:,1));
title('PrÍdkoúÊ kπtowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]'); hold on;  grid on;

%ograniczenie predkosci obrotowej
% subplot(3,2,1); plot(Przebiegi.time,wLimit, '-r');
% hleg = legend('PredkoúÊ obrotowa','Prπd maksymalny','Location','NorthEast');
% set(hleg,'FontAngle','normal'); hold off;

subplot(312); plot(Przebiegi.time, Przebiegi.signals.values(:,2));
title('GwM'); xlabel('Czas [s]'); ylabel('Odpowiedü obiektu');  grid on;

% subplot(3,2,3); plot(Przebiegi.time, Przebiegi.signals.values(:,3));
% title('GIU'); xlabel('time [s]'); ylabel('Object response');
% 
% %Pochodna prπdu
% subplot(3,2,4); plot(Przebiegi.time, Przebiegi.signals.values(:,4));
% title('dGIU/dt'); xlabel('time [s]'); ylabel('Object response'); hold on;
% 
% %Ograniczenie pochodnej prπdu
% subplot(3,2,4); plot(Przebiegi.time, derI, '-r');
% hleg = legend('Pochodna prπdu','Prπd maksymalny','Location','NorthEast');
% set(hleg,'FontAngle','normal'); hold off;

%Prπd
subplot(313); plot(Przebiegi.time, Przebiegi.signals.values(:,5));
title('GIM'); xlabel('Czas [s]'); ylabel('Odpowiedü obiektu'); hold on;  grid on;

%Ograniczenie prπdu
subplot(313); plot(Przebiegi.time, Id, '-r');
hleg = legend('Prπd','Prπd maksymalny','Location','NorthWest');
set(hleg,'FontAngle','normal');
%axis([0 1.5 0 260]);
hold off;  grid on;

%% Obliczanie nastaw ciπg≥ego regulatora prÍdkoúci

%Regulator typu P 
%Kwp = Kri 
%Kw = Mn/(psi*kz*Kt*deltawn);

%Regulator typu PI
Kw = J/(2*Kt*kz*Beta*psi);

%% Obliczanie nastaw ciπg≥ego regulatora prπdu

m = T1;
V = (Beta*(Y*Kp*B))/((B1-Beta)*Rt);

%% RozpoczÍcie symulacji uk≥adÛw z regulatorem ciπg≥ym dla rÛznych rodzajÛw obciπøenia

%% ROZRUCH UK£ADU BEZ MOMENTU OBCIAØENIA ***REGULACJA CI•G£A***
sim Model_simulink_no_Torque.slx;

f = figure; set(f,'name','2a.Brak moment obciπøenia','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(no_Torque.time,no_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(no_Torque.time,no_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(no_Torque.time,no_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(no_Torque.time,no_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(no_Torque.time,no_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

% drukujÍ powiÍkszenie w celu zweryfikowania sta≥ej czasowej uk≥adu (beta)
f = figure; set(f,'name','2a.Bez obciπøenia - powiÍkszenie','numbertitle','off', 'Position', [500, 100, 700, 300],'Color','white');
plot(no_Torque.time,no_Torque.signals.values(:,4)); hold on;  grid on;
plot(no_Torque.time,2*In*0.632*ones(length(no_Torque.time),1),'--r');
plot(0.04*ones(100,1),linspace(0,250,100),'--r'); title('Prπd silnika - weryfikacja sta≥ej czasowej(Beta)');
xlabel('Czas [s]'); ylabel('I [A]'); grid on;
xlim([0 0.1]);

%% ROZRUCH ZE ZNAMIONOWYM MOMENTEM CZYNNYM ***REGULACJA CI•G£A***
sim Model_simulink_static_Torque.slx;

f = figure; set(f,'name','2b.Znamionowy moment czynny','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(static_Torque.time,static_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(static_Torque.time,static_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(static_Torque.time,static_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(static_Torque.time,static_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(static_Torque.time,static_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

% obserwujemy nawrÛt silnika
f = figure; set(f,'name','Obserwujemy nawrÛt silnika - moment czynny','numbertitle','off', 'Position', [500, 100, 700, 400],'Color','white');
subplot(211); plot(static_Torque.time,static_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa - nawrÛt silnika dla obciπøenia znamionowym momentem czynnym');
xlabel('Czas [s]'); ylabel('\omega [rad/s]'); grid on;
axis([0 10 -1 200]);

subplot(212); plot(static_Torque.time,static_Torque.signals.values(:,3));
xlabel('Czas [s]'); ylabel('\omega [rad/s]'); grid on;
axis([0 0.08 -1 1]);

%%  PRACA UK£ADU Z UDAROWYM MOMENTEM OBCI•ØENIA ***REGULACJA CI•G£A***
sim Model_simulink_impact_Torque.slx;

f = figure; set(f,'name','2a.Udarowy moment obciπøenia','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white'); 
title('TYYYTU£');

subplot(411); plot(impact_Torque.time,impact_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(impact_Torque.time,impact_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(impact_Torque.time,impact_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(impact_Torque.time,impact_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(impact_Torque.time,impact_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

%%  PRACA UK£ADU Z BIERNYM MOMENTEM OBCI•ØENIA ***REGULACJA CI•G£A***
sim Model_simulink_passive_Torque.slx;

f = figure; set(f,'name','2c.Bierny moment obciπøenia','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white'); 
title('TYYYTU£');

subplot(411); plot(passive_Torque.time,passive_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(passive_Torque.time,passive_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(passive_Torque.time,passive_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(passive_Torque.time,passive_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(passive_Torque.time,passive_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

% obserwujemy nawrÛt silnika
f = figure; set(f,'name','Obserwujemy nawrÛt silnika - moment bierny','numbertitle','off', 'Position', [500, 100, 700, 400],'Color','white');
subplot(211); plot(passive_Torque.time,passive_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa - brak nawrotu silnika dla znamionowego momentu biernego');
xlabel('Czas [s]'); ylabel('\omega [rad/s]'); grid on;
xlim([0 10])

subplot(212); plot(passive_Torque.time,passive_Torque.signals.values(:,3));
xlabel('Czas [s]'); ylabel('\omega [rad/s]'); grid on;
axis([0 0.04 -0.05 0.1]);


%%  PRACA UK£ADU Z PULSACYJNYM MOMENTEM OBCI•ØENIA ***REGULACJA CI•G£A***
sim Model_simulink_pulse_Torque.slx;

f = figure; set(f,'name','Pulsacyjny moment obciπøenia','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white'); 
title('TYYYTU£');

subplot(411); plot(pulse_Torque.time,pulse_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(pulse_Torque.time,pulse_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(pulse_Torque.time,pulse_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(pulse_Torque.time,pulse_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(pulse_Torque.time,pulse_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

%% Punkt 3 oblcizanie transmitancji zastÍpczej uk≥adu ciπg≥ego

% reg_momentu*przekszta≥tnik*silnik*bezw≥adnoúÊ
RegPredkosc = tf([Kw*TR Kw],[TR,0]);        % transmitacja regulatora prÍdkoúci
RegMoment = tf([m 1],[V 0]);                % transmitancja regulatora momentu(prπdu)
PrzeksztTyryst = tf([Kp],[tau0 1]);         % transmitancja przekszta≥tnika tyrystorowego
Silnik = tf([B 0],[Rt*B1*T1 Rt*(B1+T1) Rt]);% transmitancja silnik
Bezwlad = tf([1],[J 0]);                    % transmitancja momentu bezwladnosci
sys = series(RegMoment,PrzeksztTyryst);     % liczymy szeregowo transmitancje regulatora predkosci i przekszta≥tnika tyrystorowego
sys = series(sys,Silnik);                   % liczymy szeregowo transmitancjÍ poprzedniego uk≥adu i silnika
sys = feedback(sys, Y);                     % liczymy sprzÍøenie zwrotne uzyskanej transmitancji oraz Y
% sys = series(psi,sys);
% sys = series(RegPredkosc,sys);
% sys = series(sys,Bezwlad);
sys = series(RegPredkosc,sys);              % liczymy szeregowo transmitancjÍ poprzedniego uk≥adu i regulatora prÍdkoúci
disp('Transmitancja otwartego uk≥adu regulacji = ');
sys = series(sys,Bezwlad)                  % liczymy szeregowo transmitancjÍ poprzedniego uk≥adu i momentu bezwladnosci
%% Charakterystyki bodego otwartego uk≥adu
f = figure; set(f,'name','3.Bode plot dla uk≥adu otwartego','numbertitle','off', 'Position', [500, 100, 700, 400],'Color','white');

bode(sys);grid on;                          % charakterystyki Bodego

%% Charakterystyka nyquista otwartego uk≥adu
P = nyquistoptions;
P.ShowFullContour = 'off'; 

f = figure; set(f,'name','3.Nyquist plot dla uk≥adu otwartego','numbertitle','off', 'Position', [500, 100, 700, 400],'Color','white');
subplot(2,1,1);
nyquist(sys,P);                       % charakterystyki nyquista
subplot(2,1,2);
nyquist(sys,P); grid on;
axis([-0.1 0 -5e-3 5e-3]); 
%% liczymy transmitancjÍ zamkniÍtego uk≥adu w celu wyznaczenia wg
disp('Transmitancja zamkniÍtego uk≥adu regulacji = ');
sys = feedback(sys, Kt)

%% Charakterystyka Bodego uk≥adu zamkniÍtego
f = figure; set(f,'name','3.Bode plot dla uk≥adu zamkniÍtego','numbertitle','off', 'Position', [500, 100, 700, 400],'Color','white');

bode(sys); grid on;                          % charakterystyki Bodego
[Gm,Pm,Wgm,Wpm] = margin(sys);

%% Charakterystyka nyquista uk≥adu zamkniÍtego oraz obliczanie zapasu fazy i wzmocnienia

f = figure; set(f,'name','3.Nyquist plot dla uk≥adu zamkniÍtego','numbertitle','off', 'Position', [500, 100, 700, 400],'Color','white');
subplot(2,1,1);
nyquist(sys,P); grid on; 
title(['\fontsize{12}Gm = ', num2str(round(20*log10(Gm),3)), ' [dB]' ,'  ,  ', 'Pm = ', num2str(round(Pm,3)), ' [∞]' ,'  ,  ', '\theta_m_a_x = ', num2str(round(Pm*pi/180/Wpm,6)) , ' s']);
subplot(2,1,2);
nyquist(sys,P); grid on;
axis([-0.1 0 -5e-3 5e-3]); 
wg = 50;                                  %[rad/s] odczytane z wykresu
wgp = wg*100;
fgp = 2*pi*wgp;
Tp1 = 1/fgp;

N=20;
Tp = Beta/N;                                
%Ti = TR;                                     % czas zdwojenie do dyskretyzacji (4*Beta)
%% Liczymy nastawy dyskretnego regulatora obrotÛw
K1 = Kw;
K2 = Kw*(Tp/TR-1);

%% Liczymy nastawy dyskretnego regulatora prπdu
K3 = m/V;
K4 = (Tp-m)/V;
%% Wykonujemy symulacjÍ dla üle dobranego(zbyt ma≥ego) czasu prÛbkowania Tp i obserwujemy zniekszta≥cenie przebiegÛw
Tp = 0.01;
Q = N*2/2e15;
% ROZRUCH UK£ADU BEZ MOMENTU OBCI•ØENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_no_torque.slx;

f = figure; set(f,'name','4.Dysktretny brak moment obciπøenia1','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]'); grid on; hold on;
xlim([0 10])

subplot(412); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on; hold on;
xlim([0 10])

subplot(413); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on; hold on;
xlim([0 10])

subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on; hold on;
xlim([0 10])


% PorÛwnanie 
Tp = Beta/N;
Q = N*2/2e15;
sim Model_simulink_dyskretny_no_torque.slx;

subplot(411); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1),'Color','magenta');
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]'); grid on;
hleg = legend('Nieprawid≥owo dobrany (T_p = 0.01)','Prawid≥owo dobrany (T_p = 0.002)','Location','NorthEast');
set(hleg,'FontAngle','normal'); 
xlim([0 10])

subplot(412); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2),'Color','magenta');
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
hleg = legend('Nieprawid≥owo dobrany (T_p = 0.01)','Prawid≥owo dobrany (T_p = 0.002)','Location','NorthEast');
set(hleg,'FontAngle','normal');
xlim([0 10])

subplot(413); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3),'Color','magenta');
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
hleg = legend('Nieprawid≥owo dobrany (T_p = 0.01)','Prawid≥owo dobrany (T_p = 0.002)','Location','NorthEast');
set(hleg,'FontAngle','normal');
xlim([0 10])

subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4),'Color','magenta');
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
hleg = legend('Nieprawid≥owo dobrany (T_p = 0.01)','Prawid≥owo dobrany (T_p = 0.002)','Location','NorthEast');
set(hleg,'FontAngle','normal');
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('Nieprawid≥owo dobrany (T_p = 0.01)','Prawid≥owo dobrany (T_p = 0.002)','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

%% Liczymy b≥πd kwantyzacji dla dobrze dobranego Tp
%% Wykonujemy symulacjÍ dla zbyt duøego kwantyzatora i obserwujemy rezultaty
Tp = Beta/N;
Q = 0.05;

% ROZRUCH UK£ADU BEZ MOMENTU OBCI•ØENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_no_torque.slx;

f = figure; set(f,'name','Dyskretny, brak moment obciπøenia2','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on; hold on;
xlim([0 10])

subplot(412); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on; hold on;
xlim([0 10])

subplot(413); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on; hold on;
xlim([0 10])

subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on; hold on;
xlim([0 10])

% PorÛwnanie
Tp = Beta/N;
Q = N*2/2e15;
sim Model_simulink_dyskretny_no_torque.slx;

subplot(411); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1),'Color','magenta');
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]'); grid on;
hleg = legend('Nieprawid≥owo dobrany (Q = 0.05)','Prawid≥owo dobrany (Q = 2e-14)','Location','NorthEast');
set(hleg,'FontAngle','normal'); 
xlim([0 10])

subplot(412); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2),'Color','magenta');
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
hleg = legend('Nieprawid≥owo dobrany (Q = 0.05)','Prawid≥owo dobrany (Q = 2e-14)','Location','NorthEast');
set(hleg,'FontAngle','normal');
xlim([0 10])

subplot(413); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3),'Color','magenta');
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
hleg = legend('Nieprawid≥owo dobrany (Q = 0.05)','Prawid≥owo dobrany (Q = 2e-14)','Location','NorthEast');
set(hleg,'FontAngle','normal');
xlim([0 10])

subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4),'Color','magenta');
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
hleg = legend('Nieprawid≥owo dobrany (Q = 0.05)','Prawid≥owo dobrany (Q = 2e-14)','Location','NorthEast');
set(hleg,'FontAngle','normal');
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('Nieprawid≥owo dobrany (Q = 0.05)','Prawid≥owo dobrany (Q = 2e-14)','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

% %% NastÍpnie wykonujemy symulacjÍ dla ma≥ego kwantyzatora i obserwujemy brak zak≥ÛceÒ
% Q = N*2/2e15;
% 
% % ROZRUCH UK£ADU BEZ MOMENTU OBCI•ØENIA ***REGULACHA DYSKRETNA***
% sim Model_simulink_dyskretny_no_torque.slx;
% 
% f = figure; set(f,'name','Dyskretny, brak moment obciπøenia3','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');
% 
% subplot(411); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
% title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
% xlim([0 10])
% 
% subplot(412); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
% title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
% xlim([0 10])
% 
% subplot(413); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
% title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
% xlim([0 10])
% 
% subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
% title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
% xlim([0 10])
% 
% %naniesienie ograniczenia
% subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
% hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');  grid on;
% set(hleg,'FontAngle','normal'); hold off;
% xlim([0 10])

%% RozpoczÍcie symulacji uk≥adÛw z regulatorem dyskretnym dla rÛønych rodzajÛw obciaøenia

%% ROZRUCH UK£ADU BEZ MOMENTU OBCI•ØENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_no_torque.slx;

f = figure; set(f,'name','Dyskretny, brak moment obciπøenia4','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');  grid on;
set(hleg,'FontAngle','normal'); hold off;
xlim([0 10])


%% ROZRUCH UK£ADU BEZ MOMENTU OBCI•ØENIA ***REGULACHA DYSKRETNA***
sim Model_simulink_dyskretny_no_torque.slx;

f = figure; set(f,'name','6. brak moment obciπøenia5','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');  grid on;
set(hleg,'FontAngle','normal'); hold off;
xlim([0 10])

%% ROZRUCH UK£ADU ZE ZNAMIONOWYM MOMENTEM OBCI•ØENIA ***REGULACHA DYSKRETNA***
Q = N*2/2e15;
deltaIn=(psi*V*Mn)/(psi^2*V+J*Kp*Y);
Uz0 = (Lambda*In+deltaIn)*(Y*B1)/(B1-Beta);
sim Model_simulink_dyskretny_static_torque.slx;

f = figure; set(f,'name','6.Znamionowy moment obciπøenia','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');  grid on;
set(hleg,'FontAngle','normal'); hold off;
xlim([0 10])

f = figure; set(f,'name','6.Przybliøenie','numbertitle','off', 'Position', [500, 000, 700, 300],'Color','white');

plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;

%naniesienie ograniczenia
plot(discrete_static_Torque.time,discrete_static_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');  grid on;
set(hleg,'FontAngle','normal'); hold off;
axis([0.7 1 280 320])

% %% DZIA£ANIE UK£ADU Z PULSACYJNYM MOMENTEM OBCI•ØENIA ***REGULACHA DYSKRETNA***
% sim Model_simulink_dyskretny_pulse_torque.slx;
% 
% f = figure; set(f,'name','Dyskretny, pulsacyjny moment obciπøenia','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');
% 
% subplot(411); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,1));
% title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
% xlim([0 10])
% 
% subplot(412); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,2));
% title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
% xlim([0 10])
% 
% subplot(413); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,3));
% title('prÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
% xlim([0 10])
% 
% subplot(414); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,4));
% title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
% %title('\frac{1}{2}','Interpreter','latex') 
% xlim([0 10])
% 
% %naniesienie ograniczenia
% subplot(414); plot(discrete_pulse_Torque.time,discrete_pulse_Torque.signals.values(:,6),'-r');
% hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
% set(hleg,'FontAngle','normal'); hold off;  grid on;
% xlim([0 10])

%% ROZRUCH UK£ADU BEZ MOMENTU OBCI•ØENIA ***REGULACHA DYSKRETNA***
Q = 1/2^2
sim Model_simulink_dyskretny_no_torque_CYKL.slx;

f = figure; set(f,'name','5.CYKL,  brak moment obciπøenia5','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(discrete_no_Torque.time,discrete_no_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');  grid on;
set(hleg,'FontAngle','normal'); hold off;
xlim([0 10])

%% CYKL GRANICZNY: ROZRUCH UK£ADU Z BIERNYM MOMENTEM OBCIAØENIA
Q = 1/2^8;
sim Model_simulink_dyskretny_passive_Torque_CYKL.slx;

f = figure; set(f,'name','5.CYKL,  bierny moment obciπøenia','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white'); 
title('TYYYTU£');

subplot(411); plot(passive_Torque.time,passive_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(passive_Torque.time,passive_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(passive_Torque.time,passive_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(passive_Torque.time,passive_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(passive_Torque.time,passive_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])


%% CYKL GRANICZNY: ROZRUCH UK£ADU Z UDAROWYM MOMENTEM OBCI•ØENIA
Q = 1/2^8;
sim Model_simulink_dyskretny_impact_torque_CYKL.slx;

f = figure; set(f,'name','5.CYKL,  udarowy moment obciπøenia','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white'); 
title('TYYYTU£');

subplot(411); plot(passive_Torque.time,passive_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(passive_Torque.time,passive_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(passive_Torque.time,passive_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(passive_Torque.time,passive_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(passive_Torque.time,passive_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

%% CYKL GRANICZNY: ROZRUCH UK£ADU ZE ZNAMIONOWYM MOMENTEM OBCIAØENIA
deltaIn=(psi*V*Mn)/(psi^2*V+J*Kp*Y);
Uz0 = (Lambda*In+deltaIn)*(Y*B1)/(B1-Beta);
Q = 1/2^8;
sim Model_simulink_dyskretny_static_torque_CYKL.slx;

f = figure; set(f,'name','5.CYKL,  znamionowy moment obciaøenia','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white'); 
title('TYYYTU£');

subplot(411); plot(passive_Torque.time,passive_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(passive_Torque.time,passive_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(passive_Torque.time,passive_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(passive_Torque.time,passive_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(passive_Torque.time,passive_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

%% DZIA£ANIE UK£ADU Z UDAROWYM MOMENTEM OBCI•ØENIA ***REGULACHA DYSKRETNA***
Tp = Beta/N;
Q = N*2/2e15;
Uz0 = (Lambda*In)*((Y*B1)/(B1-Beta));
sim Model_simulink_dyskretny_impact_torque.slx;

f = figure; set(f,'name','7.Udarowy moment obciπøenia, z kwantyzerem, DYSKRETNY','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
%title('\frac{1}{2}','Interpreter','latex') 
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');  grid on;
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])

%%  PRACA UK£ADU Z UDAROWYM MOMENTEM OBCI•ØENIA ***REGULACJA CI•G£A***
sim Model_simulink_impact_Torque.slx;

f = figure; set(f,'name','7.Udarowy moment obciπøenia CI•G£Y','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white'); 
title('TYYYTU£');

subplot(411); plot(impact_Torque.time,impact_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(impact_Torque.time,impact_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(impact_Torque.time,impact_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(impact_Torque.time,impact_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(impact_Torque.time,impact_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])


%% DZIA£ANIE UK£ADU Z UDAROWYM MOMENTEM OBCI•ØENIA ***REGULACHA DYSKRETNA***
Tp = Beta/N;
Q = N*2/2e15;
Uz0 = (Lambda*In)*((Y*B1)/(B1-Beta));
sim Model_simulink_dyskretny_impact_torque_NO_QUANTIZER.slx;

f = figure; set(f,'name','7.Udarowy moment obciπøenia, bez kwantyzera, DYSKRETNY','numbertitle','off', 'Position', [500, 000, 700, 800],'Color','white');

subplot(411); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,1));
title('NapiÍcie silnika'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(412); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,2));
title('NapiÍcie sterownika prπdu'); xlabel('Czas [s]'); ylabel('U [V]');  grid on;
xlim([0 10])

subplot(413); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,3));
title('PrÍdkoúÊ obrotowa'); xlabel('Czas [s]'); ylabel('\omega [rad/s]');  grid on;
xlim([0 10])

subplot(414); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,4));
title('Prπd silnika'); xlabel('Czas [s]'); ylabel('I [A]'); hold on;  grid on;
%title('\frac{1}{2}','Interpreter','latex') 
xlim([0 10])

%naniesienie ograniczenia
subplot(414); plot(discrete_impact_Torque.time,discrete_impact_Torque.signals.values(:,6),'-r');
hleg = legend('WartoúÊ prπdu','Prπd maksymalny','Location','NorthEast');  grid on;
set(hleg,'FontAngle','normal'); hold off;  grid on;
xlim([0 10])