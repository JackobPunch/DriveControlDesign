%% 
%% 
clear all; close all
%% Dane wejœ›ciowe
Pn = 37000; % Moc znamionowa
Un = 230;  % Znamionowe napiecie zasilania
In = 185;   % Prad znamionowy
nn = 585; % Znamionowa predkosc obrotowa
Rt = 0.12; % Rezystancja uogolniona
Lt = 1.59e-3; % Indukcyjnosc
Js = 2.1; % Moment bezwladnosci wirnika

L = Lt; % Indukcyjnosc calkowita
R = Rt; % Rezystancja uogolniona
J = 8*Js; % Moment bezwladnosci silnika i agregatu

wn = 2*pi*nn/60;

lambda_n = 2; % Stala ograniczajaca prad
p = 50; % Dopuszczalna krotnosc pradu
Beta = lambda_n/p; % Stala czasowa przebiegu pradu twornika

psi_e = (Un-In*R)/wn; % znamionowy strumieï¿½ skojarzony rotacyjnie z uzwojeniem tworika

B = J*R/(psi_e^2); % staï¿½a elektromechaniczna napï¿½du
T = L/R; % elektromagnetyczna staï¿½a czasowa
Tm = 0; % sta³a czasowa obwodu wzbudzenia

Mn=psi_e*In;  % Moment znamionowy
Mm=Mn; % Moment mechaniczny
Mb=Mn*sign(wn);

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

sim('pkt1.slx')

wn_wektor=ones(1,length(dane.time))*wn;
figure; plot(dane.time, dane.signals.values(:,1),dane.time, wn_wektor)
title('Prêdkoœæ k¹towa'); xlabel('Czas [s]'); ylabel('\omega[rad/s]'); grid on;
legend('Prêdkoœæ k¹towa', 'Ograniczenie')


figure; plot(dane.time, dane.signals.values(:,2))
title('Przebieg pochodnej pr¹du twornika'); xlabel('Czas [s]'); ylabel('I[A]'); grid on;

In_wektor=ones(1, length(dane.time))*In;
figure; plot(dane.time, dane.signals.values(:,3),dane.time, In_wektor);
title('Przebieg pr¹du twornika'); xlabel('Czas [s]'); ylabel('I[A]'); grid on;
legend('Pr¹d twornika', 'Ograniczenie')


%% Nastawy regulatora pr¹du

%Beta = lambda_n/p;

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

%% Nastawy regulatora prêdkoœci PI

TR=4*Beta;
Kw=J/(2*Kt*kz*Beta*psi_e);

%% Rozruch bez momentu obci¹¿enia
sim('symulacja.slx')
figure;
subplot(311)
plot(Us.time(:,1), Us.signals.values(:,1));
title('Przebieg napiêcia')
xlabel('Czas [s]')
ylabel('Napiêcie [V]')

subplot(312)
plot(Is.time(:,1), Is.signals.values(:,1));
title('Przebieg pr¹du twornika')
xlabel('Czas [s]')
ylabel('Pr¹d twornika [A]')

subplot(313)
plot(omega.time(:,1), omega.signals.values(:,1));
title('Przebieg prêdkoœci k¹towej')
xlabel('Czas [s]')
ylabel('Prêdkoœæ k¹towa [rad/s]')




