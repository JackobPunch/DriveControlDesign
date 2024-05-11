clear all; close all
%% Dane wej��ciowe
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

psi_e = (Un-In*R)/wn; % znamionowy strumie� skojarzony rotacyjnie z uzwojeniem tworika

B = J*R/(psi_e^2); % sta�a elektromechaniczna nap�du
T = L/R; % elektromagnetyczna sta�a czasowa
Tm = 0; % sta�a czasowa obwodu wzbudzenia
Id= lambda_n*In; % dopuszczalny pr�d twornika
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
 w_wektor=ones(1,length(dane.time))*wn; % ograniczenie pr�dko�ci
 I_wektor=ones(1, length(dane.time))*In; % ograniczenie pr�du
 I_pochodna_wektor=ones(1, length(dane.time))*p*In; % ograniczenie pochodnej pr�du

  figure;
  subplot(311)
  plot(dane.time, dane.signals.values(:,1),dane.time, w_wektor,'--')
  title('Przebieg pr�dko�ci k�towej'); xlabel('Czas [s]'); ylabel('\omega[rad/s]'); grid on;
  legend('Pr�dko�� k�towa', 'Ograniczenie')
%  
  subplot(312)
  plot(dane.time, dane.signals.values(:,2), dane.time, I_pochodna_wektor,'--')
 title('Przebieg pochodnej pr�du twornika'); xlabel('Czas [s]'); ylabel('I[A/s]'); grid on;
  legend('Pochodna pr�du twornika', 'Ograniczenie')
  subplot(313)
  plot(dane.time, dane.signals.values(:,3),dane.time, I_wektor,'--');
  title('Przebieg pr�du twornika'); xlabel('Czas [s]'); ylabel('I[A]'); grid on;
  legend('Pr�d twornika', 'Ograniczenie')


% Nastawy regulatora pr�du

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

% Nastawy regulatora pr�dko�ci PI

TR=4*Beta;
Kw=J/(2*Kt*kz*Beta*psi_e);

% Rozruch bez momentu obci��enia i obci��enie udarowe(stabilizacja)
sim('symulacja.slx')

w_wektor=ones(1,length(omega.time))*wn; % ograniczenie pr�dko�ci
I_wektor=ones(1, length(Is.time))*Id; % ograniczenie pr�du
I_pochodna_wektor=ones(1, length(Is.time))*p*In; % ograniczenie pr�du
Un_wektor=ones(1,length(Us.time))*Un; % ograniczenie napi�cia

% Przebieg napi��ia
f=figure;
set(f,'name','Rozruch bez momentu obci��eniia i obci��enie udarowe(stabilizacja)')
subplot(311)
plot(Us.time(:,1), Us.signals.values(:,1), Us.time, Un_wektor,'--');
title('Przebieg napi�cia')
xlabel('Czas [s]')
ylabel('Napi�cie [V]')
legend('Napi�cie', 'Ograniczenie')
grid on;
% 
%Przebieg pr�du twornika
subplot(312)
plot(Is.time(:,1), Is.signals.values(:,1), Is.time, I_wektor,'--');
title('Przebieg pr�du twornika')
xlabel('Czas [s]')
ylabel('Pr�d twornika [A]')
legend('Pr�d', 'Ograniczenie')
grid on;
% 
% Przebieg pr�dko�ci k�towej
subplot(313)
plot(omega.time(:,1), omega.signals.values(:,1), omega.time, w_wektor,'--');
title('Przebieg pr�dko�ci k�towej')
xlabel('Czas [s]')
ylabel('Pr�dko�� k�towa [rad/s]')
legend('Pr�dko�� k�towa', 'Ograniczenie')
grid on;
% Rozruch ze znamionowym momentem czynnym
 sim('symulacja2.slx')
 
 w_wektor=ones(1,length(omega.time))*wn; % ograniczenie pr�dko�ci
 I_wektor=ones(1, length(Is.time))*In; % ograniczenie pr�du
 Un_wektor=ones(1,length(Us.time))*Un; % ograniczenie napi�cia
 
% Przebieg napi��ia
 f=figure;
 set(f,'name', 'Rozruch ze znamionowym momentem czynnym')
 subplot(311)
 plot(Us.time(:,1), Us.signals.values(:,1), Us.time, Un_wektor,'--');
 title('Przebieg napi�cia')
 xlabel('Czas [s]')
 ylabel('Napi�cie [V]')
 legend('Napi�cie', 'Ograniczenie')
 grid on;
% 
% Przebieg pr�du twornika
 subplot(312)
 plot(Is.time(:,1), Is.signals.values(:,1), Is.time, I_wektor,'--');
 title('Przebieg pr�du twornika')
 xlabel('Czas [s]')
 ylabel('Pr�d twornika [A]')
 legend('Pr�d', 'Ograniczenie')
 grid on;
% 
% Przebieg pr�dko�ci k�towej
 subplot(313)
 plot(omega.time(:,1), omega.signals.values(:,1), omega.time, w_wektor,'--');
 title('Przebieg pr�dko�ci k�towej')
 xlabel('Czas [s]')
 ylabel('Pr�dko�� k�towa [rad/s]')
 legend('Pr�dko�� k�towa', 'Ograniczenie')
 grid on;
% 
% % Rozruch ze znamionowym momentem biernym
 sim('symulacja1.slx')
% 
 w_wektor=ones(1,length(omega.time))*wn; % ograniczenie pr�dko�ci
 I_wektor=ones(1, length(Is.time))*In; % ograniczenie pr�du
 Un_wektor=ones(1,length(Us.time))*Un; % ograniczenie napi�cia
 
% Przebieg napi��ia
 f=figure;
 set(f,'name', 'Rozruch ze znamionowym momentem biernym')
 subplot(311)
 plot(Us.time(:,1), Us.signals.values(:,1), Us.time, Un_wektor,'--');
 title('Przebieg napi�cia')
 xlabel('Czas [s]')
 ylabel('Napi�cie [V]')
 legend('Napi�cie', 'Ograniczenie')
 grid on;
% 
% Przebieg pr�du twornika
 subplot(312)
 plot(Is.time(:,1), Is.signals.values(:,1), Is.time, I_wektor,'--');
 title('Przebieg pr�du twornika')
 xlabel('Czas [s]')
 ylabel('Pr�d twornika [A]')
 legend('Pr�d', 'Ograniczenie')
 grid on;
% 
% Przebieg pr�dko�ci k�towej
 subplot(313)
 plot(omega.time(:,1), omega.signals.values(:,1), omega.time, w_wektor,'--');
 title('Przebieg pr�dko�ci k�towej')
 xlabel('Czas [s]')
 ylabel('Pr�dko�� k�towa [rad/s]')
 legend('Pr�dko�� k�towa', 'Ograniczenie')
 grid on;
% 
% % Zapas fazy i modu�u
% 
%reg_momentu*przeksztaltnik*silnik*bezwladnosci
tau0=3.3e-3
RegPredkosc = tf([Kw*TR Kw],[TR,0]);        % transmitacja regulatora pr�dko�ci
RegMoment = tf([m 1],[V 0]);                % transmitancja regulatora momentu(pr�du)
PrzeksztTyryst = tf([Kp],[tau0 1]);         % transmitancja przekszta�tnika tyrystorowego
Silnik = tf([B 0],[Rt*B1*T1 Rt*(B1+T1) Rt]);% transmitancja silnik
Bezwlad = tf([1],[J 0]);                    % transmitancja momentu bezwladnosci
uklad = series(RegMoment,PrzeksztTyryst);     % liczymy szeregowo transmitancje regulatora predkosci i przekszta�tnika tyrystorowego
uklad = series(uklad,Silnik);                   % liczymy szeregowo transmitancj� poprzedniego uk�adu i silnika
uklad = feedback(uklad, Y);                     % liczymy sprz�enie zwrotne uzyskanej transmitancji oraz Y
uklad = series(RegPredkosc,uklad);              % liczymy szeregowo transmitancj� poprzedniego uk�adu i regulatora pr�dko�ci
disp('Transmitancja otwartego uk�adu regulacji: ');
uklad = series(uklad,Bezwlad)                  % liczymy szeregowo transmitancj� poprzedniego uk�adu i momentu bezwladnosci
% 
% Charakterystyki bodego otwartego uk�adu
 f = figure; set(f,'name','Bode plot dla uk�adu otwartego','numbertitle','off');
 bode(uklad);grid on;                          % charakterystyki Bodego
 
 %Charakterystyka nyquista otwartego uk�adu
 f = figure; set(f,'name','Nyquist plot dla uk�adu otwartego','numbertitle','off');
 subplot(2,1,1);
 nyquist(uklad);grid on;                       % charakterystyki nyquista
 subplot(2,1,2);
 nyquist(uklad); axis([-0.85 0 -0.1 0.1]);  
% 
% liczymy transmitancj� zamkni�tego uk�adu w celu wyznaczenia wg
 disp('Transmitancja zamkni�tego uk�adu regulacji: ');
 uklad = feedback(uklad, Kt)
 
 %Charakterystyka Bodego uk�adu zamkni�tego
 f = figure; set(f,'name','Bode plot dla uk�adu zamkni�tego','numbertitle','off');
 bode(uklad);grid on;                          % charakterystyki Bodego
 [Gm,Pm,Wgm,Wpm] = margin(uklad);
 
 %Charakterystyka nyquista uk�adu zamkni�tego oraz obliczanie zapasu fazy i wzmocnienia
 f = figure; set(f,'name','Nyquist plot dla uk�adu zamkni�tego','numbertitle','off');
 subplot(2,1,1);
 nyquist(uklad);grid on; title(sprintf('Zapas wzmocnienia %f[razy]\nZapas fazy %f[stopni]\nMaksymalne op�nienie %f[s]',Gm,Pm,Pm*pi/180/Wpm));
 subplot(2,1,2);
 nyquist(uklad); axis([-0.85 0 -0.1 0.1]); 
%            
% % Dokonac dyskretyzacji regulator�w dzia�ania ci �ag�ego i wykona � c symulacje dla r��znych czas�w pr�b- �
% wz = 1.25;                                  %[rad/s] odczytane z wykresu
% zwp = wz*10;
% fgp = 2*pi*zwp;
% Tp = 1/fgp;
% 
% N=20;
% Tp = Beta/N; 
 
% Nastawy dyskretnego regulatora pr�dko�ci obrotowej
% K1 = Kw;
% K2 = Kw*(Tp/TR-1);
% Nastawy dyskretnego regulatora pr�du
% K3 = m/V;
% K4 = (Tp-m)/V;
% Symulacja dla �le dobranego (zbyt du�ego) czasu pr�bkowania
% Tp = 0.5;
% sim('symulacja_dyskretny.slx')
% 
