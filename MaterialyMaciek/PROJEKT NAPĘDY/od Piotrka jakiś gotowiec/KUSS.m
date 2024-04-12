clear all; close all
%% Dane wejœ›ciowe
Pn = 26e3; % Moc znamionowa
Un = 230;  % Znamionowe napiecie zasilania
In = 132;   % Prad znamionowy
nn = 625; % Znamionowa predkosc obrotowa
Rt = 0.202; % Rezystancja uogolniona
Lt = 1.9e-3; % Indukcyjnosc
Js = 1.25; % Moment bezwladnosci wirnika

L = Lt; % Indukcyjnosc calkowita
R = Rt; % Rezystancja uogolniona
J = 29*Js; % Moment bezwladnosci silnika i agregatu

wn = 2*pi*nn/60;

lambda_n = 1.8; % Stala ograniczajaca prad
p = 50; % Dopuszczalna krotnosc pradu
Beta = lambda_n/p; % Stala czasowa przebiegu pradu twornika

psi_e = (Un-In*R)/wn; % znamionowy strumieï¿½ skojarzony rotacyjnie z uzwojeniem tworika

B = J*R/(psi_e^2); % sta³a elektromechaniczna napêdu
T = L/R; % elektromagnetyczna staï¿½a czasowa
Tm = 0; % sta³a czasowa obwodu wzbudzenia
Id= lambda_n*In; % dopuszczalny pr¹d twornika
Mn=psi_e*In;  % Moment znamionowy
Mm=Mn; % Moment mechaniczny


%Sprawdzanie warunku B > 4T
if B > 4*T
    disp('Ok')
else
    disp('Zle')
end

% Obliczenie transmitancji
G_wU = tf([0 0 1/psi_e],[B*T B 1]);
G_wM = tf([0 Rt/(psi_e^2)*T Rt/(psi_e^2)],[B*T B 1]);
G_IU = tf([0 B/Rt 0],[B*T B 1]);
G_IM = tf([0 0 1/psi_e],[B*T B 1]);
% 
%sim('odpimp.slx')
  sim('pkt1.slx')
% 
 w_wektor=ones(1,length(dane.time))*wn; % ograniczenie prêdkoœci
 I_wektor=ones(1, length(dane.time))*In; % ograniczenie pr¹du
 I_pochodna_wektor=ones(1, length(dane.time))*p*In; % ograniczenie pochodnej pr¹du

  figure;
  subplot(311)
  plot(dane.time, dane.signals.values(:,1),dane.time, w_wektor,'--')
  title('Przebieg prêdkoœci k¹towej'); xlabel('Czas [s]'); ylabel('\omega[rad/s]'); grid on;
  legend('Prêdkoœæ k¹towa', 'Ograniczenie')
%  
  subplot(312)
  plot(dane.time, dane.signals.values(:,2), dane.time, I_pochodna_wektor,'--')
 title('Przebieg pochodnej pr¹du twornika'); xlabel('Czas [s]'); ylabel('I[A/s]'); grid on;
  legend('Pochodna pr¹du twornika', 'Ograniczenie')
  subplot(313)
  plot(dane.time, dane.signals.values(:,3),dane.time, I_wektor,'--');
  title('Przebieg pr¹du twornika'); xlabel('Czas [s]'); ylabel('I[A]'); grid on;
  legend('Pr¹d twornika', 'Ograniczenie')


% Nastawy regulatora pr¹du

Beta = lambda_n/p;

Y=10/(2.5*In);
Kt=10/(1.2*wn);
wz=wn*Kt;
Kp=1.5*Un/10;
T1=0.5*B*(1-sqrt(1-4*T/B));
B1=B-T1;
kz=(B1-Beta)/(Y*B1);
m=T1;
V=Beta*Y*Kp*B/((B1-Beta)*Rt);
Uz0=lambda_n*In*Y*B1/(B1-Beta);

% Nastawy regulatora prêdkoœci PI

TR=4*Beta;
Kw=J/(2*Kt*kz*Beta*psi_e);

% Rozruch bez momentu obci¹¿enia i obci¹¿enie udarowe(stabilizacja)
sim('symulacja.slx')

w_wektor=ones(1,length(omega.time))*wn; % ograniczenie prêdkoœci
I_wektor=ones(1, length(Is.time))*Id; % ograniczenie pr¹du
I_pochodna_wektor=ones(1, length(Is.time))*p*In; % ograniczenie pr¹du
Un_wektor=ones(1,length(Us.time))*Un; % ograniczenie napiêcia

% Przebieg napiêæia
f=figure;
set(f,'name','Rozruch bez momentu obci¹¿eniia i obci¹¿enie udarowe(stabilizacja)')
subplot(311)
plot(Us.time(:,1), Us.signals.values(:,1), Us.time, Un_wektor,'--');
title('Przebieg napiêcia')
xlabel('Czas [s]')
ylabel('Napiêcie [V]')
legend('Napiêcie', 'Ograniczenie')
grid on;
% 
%Przebieg pr¹du twornika
subplot(312)
plot(Is.time(:,1), Is.signals.values(:,1), Is.time, I_wektor,'--');
title('Przebieg pr¹du twornika')
xlabel('Czas [s]')
ylabel('Pr¹d twornika [A]')
legend('Pr¹d', 'Ograniczenie')
grid on;
% 
% Przebieg prêdkoœci k¹towej
subplot(313)
plot(omega.time(:,1), omega.signals.values(:,1), omega.time, w_wektor,'--');
title('Przebieg prêdkoœci k¹towej')
xlabel('Czas [s]')
ylabel('Prêdkoœæ k¹towa [rad/s]')
legend('Prêdkoœæ k¹towa', 'Ograniczenie')
grid on;
% Rozruch ze znamionowym momentem czynnym
 sim('symulacja2.slx')
 
 w_wektor=ones(1,length(omega.time))*wn; % ograniczenie prêdkoœci
 I_wektor=ones(1, length(Is.time))*In; % ograniczenie pr¹du
 Un_wektor=ones(1,length(Us.time))*Un; % ograniczenie napiêcia
 
% Przebieg napiêæia
 f=figure;
 set(f,'name', 'Rozruch ze znamionowym momentem czynnym')
 subplot(311)
 plot(Us.time(:,1), Us.signals.values(:,1), Us.time, Un_wektor,'--');
 title('Przebieg napiêcia')
 xlabel('Czas [s]')
 ylabel('Napiêcie [V]')
 legend('Napiêcie', 'Ograniczenie')
 grid on;
% 
% Przebieg pr¹du twornika
 subplot(312)
 plot(Is.time(:,1), Is.signals.values(:,1), Is.time, I_wektor,'--');
 title('Przebieg pr¹du twornika')
 xlabel('Czas [s]')
 ylabel('Pr¹d twornika [A]')
 legend('Pr¹d', 'Ograniczenie')
 grid on;
% 
% Przebieg prêdkoœci k¹towej
 subplot(313)
 plot(omega.time(:,1), omega.signals.values(:,1), omega.time, w_wektor,'--');
 title('Przebieg prêdkoœci k¹towej')
 xlabel('Czas [s]')
 ylabel('Prêdkoœæ k¹towa [rad/s]')
 legend('Prêdkoœæ k¹towa', 'Ograniczenie')
 grid on;
% 
% % Rozruch ze znamionowym momentem biernym
 sim('symulacja1.slx')
% 
 w_wektor=ones(1,length(omega.time))*wn; % ograniczenie prêdkoœci
 I_wektor=ones(1, length(Is.time))*In; % ograniczenie pr¹du
 Un_wektor=ones(1,length(Us.time))*Un; % ograniczenie napiêcia
 
% Przebieg napiêæia
 f=figure;
 set(f,'name', 'Rozruch ze znamionowym momentem biernym')
 subplot(311)
 plot(Us.time(:,1), Us.signals.values(:,1), Us.time, Un_wektor,'--');
 title('Przebieg napiêcia')
 xlabel('Czas [s]')
 ylabel('Napiêcie [V]')
 legend('Napiêcie', 'Ograniczenie')
 grid on;
% 
% Przebieg pr¹du twornika
 subplot(312)
 plot(Is.time(:,1), Is.signals.values(:,1), Is.time, I_wektor,'--');
 title('Przebieg pr¹du twornika')
 xlabel('Czas [s]')
 ylabel('Pr¹d twornika [A]')
 legend('Pr¹d', 'Ograniczenie')
 grid on;
% 
% Przebieg prêdkoœci k¹towej
 subplot(313)
 plot(omega.time(:,1), omega.signals.values(:,1), omega.time, w_wektor,'--');
 title('Przebieg prêdkoœci k¹towej')
 xlabel('Czas [s]')
 ylabel('Prêdkoœæ k¹towa [rad/s]')
 legend('Prêdkoœæ k¹towa', 'Ograniczenie')
 grid on;
% 
% % Zapas fazy i modu³u
% 
%reg_momentu*przeksztaltnik*silnik*bezwladnosci
tau0=3.3e-3
RegPredkosc = tf([Kw*TR Kw],[TR,0]);        % transmitacja regulatora prêdkoœci
RegMoment = tf([m 1],[V 0]);                % transmitancja regulatora momentu(pr¹du)
PrzeksztTyryst = tf([Kp],[tau0 1]);         % transmitancja przekszta³tnika tyrystorowego
Silnik = tf([B 0],[Rt*B1*T1 Rt*(B1+T1) Rt]);% transmitancja silnik
Bezwlad = tf([1],[J 0]);                    % transmitancja momentu bezwladnosci
uklad = series(RegMoment,PrzeksztTyryst);     % liczymy szeregowo transmitancje regulatora predkosci i przekszta³tnika tyrystorowego
uklad = series(uklad,Silnik);                   % liczymy szeregowo transmitancjê poprzedniego uk³adu i silnika
uklad = feedback(uklad, Y);                     % liczymy sprzê¿enie zwrotne uzyskanej transmitancji oraz Y
uklad = series(RegPredkosc,uklad);              % liczymy szeregowo transmitancjê poprzedniego uk³adu i regulatora prêdkoœci
disp('Transmitancja otwartego uk³adu regulacji: ');
uklad = series(uklad,Bezwlad)                  % liczymy szeregowo transmitancjê poprzedniego uk³adu i momentu bezwladnosci
% 
% Charakterystyki bodego otwartego uk³adu
 f = figure; set(f,'name','Bode plot dla uk³adu otwartego','numbertitle','off');
 bode(uklad);grid on;                          % charakterystyki Bodego
 
 %Charakterystyka nyquista otwartego uk³adu
 f = figure; set(f,'name','Nyquist plot dla uk³adu otwartego','numbertitle','off');
 subplot(2,1,1);
 nyquist(uklad);grid on;                       % charakterystyki nyquista
 subplot(2,1,2);
 nyquist(uklad); axis([-0.85 0 -0.1 0.1]);  
% 
% liczymy transmitancjê zamkniêtego uk³adu w celu wyznaczenia wg
 disp('Transmitancja zamkniêtego uk³adu regulacji: ');
 uklad = feedback(uklad, Kt)
 
 %Charakterystyka Bodego uk³adu zamkniêtego
 f = figure; set(f,'name','Bode plot dla uk³adu zamkniêtego','numbertitle','off');
 bode(uklad);grid on;                          % charakterystyki Bodego
 [Gm,Pm,Wgm,Wpm] = margin(uklad);
 
 %Charakterystyka nyquista uk³adu zamkniêtego oraz obliczanie zapasu fazy i wzmocnienia
 f = figure; set(f,'name','Nyquist plot dla uk³adu zamkniêtego','numbertitle','off');
 subplot(2,1,1);
 nyquist(uklad);grid on; title(sprintf('Zapas wzmocnienia %f[razy]\nZapas fazy %f[stopni]\nMaksymalne opóŸnienie %f[s]',Gm,Pm,Pm*pi/180/Wpm));
 subplot(2,1,2);
 nyquist(uklad); axis([-0.85 0 -0.1 0.1]); 
%            
% % Dokonac dyskretyzacji regulatorów dzia³ania ci ²ag³ego i wykona ´ c symulacje dla róÿznych czasów prób- ´
% wz = 1.25;                                  %[rad/s] odczytane z wykresu
% zwp = wz*10;
% fgp = 2*pi*zwp;
% Tp = 1/fgp;
% 
% N=20;
% Tp = Beta/N; 
 
% Nastawy dyskretnego regulatora prêdkoœci obrotowej
% K1 = Kw;
% K2 = Kw*(Tp/TR-1);
% Nastawy dyskretnego regulatora pr¹du
% K3 = m/V;
% K4 = (Tp-m)/V;
% Symulacja dla Ÿle dobranego (zbyt du¿ego) czasu próbkowania
% Tp = 0.5;
% sim('symulacja_dyskretny.slx')
% 
