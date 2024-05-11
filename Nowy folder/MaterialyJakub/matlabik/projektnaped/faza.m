%% Dane wejś›ciowe
Pn = 17e3; % Moc znamionowa
Un = 220;  % Znamionowe napiecie zasilania
In = 88;   % Prad znamionowy
nn = 1500; % Znamionowa predkosc obrotowa
Rt = 0.15; % Rezystancja uogolniona
Lt=0.01875; % Indukcyjnosc
Js = 0.275; % Moment bezwladnosci wirnika

L = Lt; % Indukcyjnosc calkowita
R = Rt; % Rezystancja uogolniona
J = 22*Js; % Moment bezwladnosci silnika i agregatu

wn = 2*pi*nn/60;

lambda_n = 2; % Stala ograniczajaca prad
p = 50; % Dopuszczalna krotnosc pradu
Beta = lambda_n/p; % Stala czasowa przebiegu pradu twornika

psi_e = (Un-In*R)/wn; % znamionowy strumieďż˝ skojarzony rotacyjnie z uzwojeniem tworika

B = J*R/(psi_e^2); % stała elektromechaniczna napędu
T = L/R; % elektromagnetyczna staďż˝a czasowa
Tm = 0; % stała czasowa obwodu wzbudzenia
Id= lambda_n*In; % dopuszczalny prąd twornika
Mn=psi_e*In;  % Moment znamionowy
Mm=Mn; % Moment mechaniczny


%Sprawdzanie warunku B > 4T
if B > 4*T
    disp('Ok')
else
    disp('Zle')
end

%% Obliczenie transmitancji
G_wU = tf([0 0 1/psi_e],[B*T B 1]);
G_wM = tf([0 Rt/(psi_e^2)*T Rt/(psi_e^2)],[B*T B 1]);
G_IU = tf([0 B/Rt 0],[B*T B 1]);
G_IM = tf([0 0 1/psi_e],[B*T B 1]);
%% Nastawy regulatora prądu

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
Uz0=0.937*lambda_n*In*Y*B1/(B1-Beta)
%Uz0=lambda_n*In*Y*B1/(B1-Beta)
% Nastawy regulatora prędkości PI
dI=(psi_e*V*Mn)/(psi_e^2*V+J*Kp*Y)
TR=4*Beta;
 Kw=J/(2*Kt*kz*Beta*psi_e);
 %% Nastawy regulatora prędkości PI
dI=(psi_e*V*Mn)/(psi_e^2*V+J*Kp*Y)
TR=4*Beta;
 Kw=J/(2*Kt*kz*Beta*psi_e);
%% Zapas fazy i modułu

% reg_momentu*przekształtnik*silnik*bezwładność
tau0=8.1e-3
RegPredkosc = tf([Kw*TR Kw],[TR 0]);        % transmitacja regulatora prędkości
RegMoment = tf([m 1],[V 0]);                % transmitancja regulatora momentu(prądu)
PrzeksztTyryst = tf([Kp],[tau0 1]);         % transmitancja przekształtnika tyrystorowego
Silnik = tf([1/R],[T 1]);                   % transmitancja silnik
Strumien = tf([psi_e]);
Bezwlad = tf([1],[J 0]);                    % transmitancja momentu bezwladnosci
odwrBezwlad = tf([J 0]);
prad = tf([Y]);
Ktt = tf([Kt]);
odwrStrumien = tf([1],[psi_e])
spr = tf([Y*J 0], [0 psi_e]);
uklad = series(Silnik,Strumien);
uklad = series(uklad,Bezwlad);
uklad = feedback(uklad,Strumien);
uklad = series(uklad,PrzeksztTyryst);
uklad = series(uklad,RegMoment);
uklad = feedback(uklad,spr);
disp('Transmitancja otwartego układu regulacji: ');
uklad = series(uklad,RegPredkosc)                  % liczymy szeregowo transmitancję poprzedniego układu i momentu bezwladnosci

% Charakterystyki bodego otwartego układu
f = figure; set(f,'name','Bode plot dla układu otwartego','numbertitle','off');
bode(uklad);grid on;                          % charakterystyki Bodego

% Charakterystyka nyquista otwartego układu
f = figure; set(f,'name','Nyquist plot dla układu otwartego','numbertitle','off');
subplot(2,1,1);
 nyquist(uklad);grid on; title(sprintf('Zapas wzmocnienia %f[dB]\nZapas fazy %f[°]\nMaksymalne opóźnienie %f[s]',Gm,Pm,Pm*pi/180/Wpm));                     % charakterystyki nyquista
subplot(2,1,2);
nyquist(uklad); 
axis([-1.1 0.1 -1 1]);  

% % wyznaczenie zapasu modulu i fazy
% f = figure; set(f,'name','Nyquist plot dla układu otwartego','numbertitle','off');
% [Gm,Pm,Wgm,Wpm] = margin(uklad)
% subplot(111); nyquist(uklad);grid on; title(sprintf('Zapas wzmocnienia %f[dB]\nZapas fazy %f[°]\nMaksymalne opóźnienie %f[s]',Gm,Pm,Pm*pi/180/Wpm));

% liczymy transmitancję zamkniętego układu w celu wyznaczenia wg
disp('Transmitancja zamkniętego układu regulacji: ');
uklad = feedback(uklad, Ktt)

% Charakterystyka Bodego układu zamkniętego
f = figure; set(f,'name','Bode plot dla układu zamkniętego','numbertitle','off');
bode(uklad);grid on;                          % charakterystyki Bodego
[Gm,Pm,Wgm,Wpm] = margin(uklad)

% Charakterystyka nyquista układu zamkniętego oraz obliczanie zapasu fazy i wzmocnienia
f = figure; set(f,'name','Nyquist plot dla układu zamkniętego','numbertitle','off');
subplot(2,1,1);
nyquist(uklad);grid on; title(sprintf('Zapas wzmocnienia %f[dB]\nZapas fazy %f[°]\nMaksymalne opóźnienie %f[s]',Gm,Pm,Pm*pi/180/Wpm));
subplot(2,1,2);
nyquist(uklad); 
axis([-1.1 0.1 -1 1]); 
