%% Komputerowe uk³ady sterowania
% Projekt 2018  
% Autorzy: Tomasz Kielar, Tomasz Telesz  
% EAIiIB Elektrotechnika rok III


%% 1. Czyszczenie przestrzeni roboczej
 clc, clear all, close all

%% 2. Dane
Pn = 31000;  %moc znamionowa [W]
Un = 440;   %napiêcie znamionowe [V]
In = 77;    %pr¹d znamionowy [A]
Nn = 1400;  %prêdkoœæ znamionowa [obr/min]
R = 0.253;  %rezystancja twornika [Ohm]
L = 1.9*10^-3;   %indukcujnosc twornika [H]
Js = 0.75;  %moment bezwladnosci [kg*m^2]

%% 3. Prze³¹cznia uk³adu
regpred = 1;        %1 - regulator PI, 0 - reg P
przelacz3 = 0;      %1 - za³¹cza moment obci¹¿enia
delomega = 0;       %1 - deltaomega 2%, 0 - 5%

%% 4. Sprawdzenie warunków
fn = Nn/60;     %czêstotliwoœæ znamionowa
omegan = 2*pi*fn;   %pulsacja znamionowa
lambdan = 2;
p = 50;
T = L/R;    %elektromagnetyczna sta³a czasowa
J = 6*Js;   %zastêpszy moment bezw³adnoœci
psien = (Un-R*In)/omegan;   %znamionowy strumieñ
B = J*R/psien^2;    %elektromechaniczna sta³a czasowa

%% 5. Obliczenia
tau0 = 3.3e-3;  %[s] sta³a czasowa przekszta³tnika tyrystorowego
Kp = 1.5*Un/10;    %wzmocnienie wzmacniacza mocy
Y = 7.8/(2.5*In);     %zmienione z 10 na 7.8, bo nie spe³nia³o warunku z Uz0
Kt = 10/(1.2*omegan);   %wzmocnienie
Mn = psien*In;  %znamionowy moment obci¹¿enia
Top = 0;  %opóŸnienie za³¹czanie znamionowego momentu obci¹¿enia
if (delomega == 0)
    deltaomega = 0.05*omegan;   %statyzm regulacji
else
    deltaomega = 0.02*omegan;
end;
beta = lambdan/p;   %sta³a czasowa przebiegu pr¹du twornika
T1 = 0.5*B*(1-sqrt(1-4*T/B));
B1 = B-T1;
kz = (B1-beta)/(Y*B1);  %zastêpczy wspó³czynnik wzmocnienia
if (B > 4*T)
    m = T1;                         %parametry regulatora
    V = beta*Y*Kp*B/((B1-beta)*R);  %
else
    m = sqrt(B*T);
    V = Kp*Y/R*B*beta/(sqrt(B*T)-beta);
end
Uz0 = lambdan*In*Y*B1/(B1-beta);    %ograniczenie Uz0 < 10V

TR = 4*beta;  %parametr wg kryterium symetrycznego
if (regpred == 0)
    Kw = Mn/(psien*kz*Kt*deltaomega);    %regulator P
else
    Kw = J/(2*Kt*kz*beta*psien);    %regulator PI
end;
time = linspace(0,1,2);
Imax = In*lambdan+0*time;

%% 6. Odpowiedzi modelów w postaci transmitancji na skok impilsow
sim('Projekt_KUS_sim');
figure(1)
grid on
hold on
title('GIU');
xlabel('Czas [s]');
ylabel('Pr¹d [I]');
ylim([0 1800]);
plot(GUI_imp.time,GUI_imp.signals.values(:,1));
plot(time,Imax,'r--');

figure(2)
grid on
hold on
title('GIU pochodna');
xlabel('Czas [s]');
ylabel('Pochodna pr¹du [I/s]');
%ylim([-0.5, 28.e4]);
plot(GUI_der_imp.time,GUI_der_imp.signals.values(:,1));
plot(GUI_der_imp.time, p*In*ones(length(GUI_der_imp.time),1),'r--');
xlim([0 0.4]);

figure(3)
grid on
hold on
%title('GomegaU');
xlabel('Czas [s]');
ylabel('Prêdkoœæ k¹towa [rad/s]');
%ylim([0 60]);
plot(GomegaU_imp.time,GomegaU_imp.signals.values(:,1));
plot(GomegaM_imp.time, omegan*ones(length(GUI_der_imp.time),1),'r--');

figure(4)
grid on
hold on
title('GIM');
xlabel('Czas [s]');
ylabel('Pr¹d [I]');
%ylim([0 40]);
plot(GIM_imp.time,GIM_imp.signals.values(:,1));
plot(time,Imax,'r');

%% 7. Rysowanie przebiegów symulacyjnych - uk³ad analogowy
regpred = 1;    %regulator PI
przelacz3 = 0;
sim('projekt');

figure(5)
subplot(311)
grid on
hold on
title('Napiêcie silnika');
xlabel('Czas [s]');
ylabel('Napiêcie [V]');
plot(Us.time,Us.signals.values(:,1));
xlim([0 3]);

subplot(312)
grid on
hold on
title('Pr¹d silnika');
xlabel('Czas [s]');
ylabel('Pr¹d [A]');
plot(Is.time,Is.signals.values(:,1));
plot(Imax.time,Imax.signals.values(:,1),'r');
xlim([0 3]);

subplot(313)
grid on
hold on
title('Pulsacja silnika');
xlabel('Czas [s]');
ylabel('Pulsacja [rad/s]');
plot(omega.time,omega.signals.values(:,1));
xlim([0 3]);
saveas(gcf, 'PI_bez_mom.png');

regpred = 1;    %regulator PI
przelacz3 = 1;
sim('projekt');

figure(6)
subplot(311)
grid on
hold on
title('Napiêcie silnika');
xlabel('Czas [s]');
ylabel('Napiêcie [V]');
plot(Us.time,Us.signals.values(:,1));
xlim([0 3]);

subplot(312)
grid on
hold on
title('Pr¹d silnika');
xlabel('Czas [s]');
ylabel('Pr¹d [A]');
plot(Is.time,Is.signals.values(:,1));
plot(Imax.time,Imax.signals.values(:,1),'r');
xlim([0 3]);

subplot(313)
grid on
hold on
title('Pulsacja silnika');
xlabel('Czas [s]');
ylabel('Pulsacja [rad/s]');
plot(omega.time,omega.signals.values(:,1));
xlim([0 3]);
saveas(gcf, 'PI_z_mom.png');

regpred = 1;    %regulator PI
przelacz3 = 1;
Top = 2;
sim('projekt');

figure(7)
subplot(311)
grid on
hold on
title('Napiêcie silnika');
xlabel('Czas [s]');
ylabel('Napiêcie [V]');
plot(Us.time,Us.signals.values(:,1));
xlim([0 3]);

subplot(312)
grid on
hold on
title('Pr¹d silnika');
xlabel('Czas [s]');
ylabel('Pr¹d [A]');
plot(Is.time,Is.signals.values(:,1));
plot(Imax.time,Imax.signals.values(:,1),'r');
xlim([0 3]);

subplot(313)
grid on
hold on
title('Pulsacja silnika');
xlabel('Czas [s]');
ylabel('Pulsacja [rad/s]');
plot(omega.time,omega.signals.values(:,1));
xlim([0 3]);
saveas(gcf, 'PI_z_udarem.png');


regpred = 0;    %regulator P
Top = 0;
przelacz3 = 0;

sim('projekt');

figure(8)
subplot(311)
grid on
hold on
title('Napiêcie silnika');
xlabel('Czas [s]');
ylabel('Napiêcie [V]');
plot(Us.time,Us.signals.values(:,1));
xlim([0 3]);

subplot(312)
grid on
hold on
title('Pr¹d silnika');
xlabel('Czas [s]');
ylabel('Pr¹d [A]');
plot(Is.time,Is.signals.values(:,1));
plot(Imax.time,Imax.signals.values(:,1),'r');
xlim([0 3]);

subplot(313)
grid on
hold on
title('Pulsacja silnika');
xlabel('Czas [s]');
ylabel('Pulsacja [rad/s]');
plot(omega.time,omega.signals.values(:,1));
xlim([0 3]);
saveas(gcf, 'P_bez_mom.png');

regpred = 0;    %regulator P
Top = 0;
przelacz3 = 1;

sim('projekt');

figure(9)
subplot(311)
grid on
hold on
title('Napiêcie silnika');
xlabel('Czas [s]');
ylabel('Napiêcie [V]');
plot(Us.time,Us.signals.values(:,1));
xlim([0 3]);

subplot(312)
grid on
hold on
title('Pr¹d silnika');
xlabel('Czas [s]');
ylabel('Pr¹d [A]');
plot(Is.time,Is.signals.values(:,1));
plot(Imax.time,Imax.signals.values(:,1),'r');
xlim([0 3]);

subplot(313)
grid on
hold on
title('Pulsacja silnika');
xlabel('Czas [s]');
ylabel('Pulsacja [rad/s]');
plot(omega.time,omega.signals.values(:,1));
xlim([0 3]);
saveas(gcf, 'P_z_mom.png');

regpred = 0;    %regulator P
Top = 2;
przelacz3 = 1;

sim('projekt');

figure(10)
subplot(311)
grid on
hold on
title('Napiêcie silnika');
xlabel('Czas [s]');
ylabel('Napiêcie [V]');
plot(Us.time,Us.signals.values(:,1));
xlim([0 3]);

subplot(312)
grid on
hold on
title('Pr¹d silnika');
xlabel('Czas [s]');
ylabel('Pr¹d [A]');
plot(Is.time,Is.signals.values(:,1));
plot(Imax.time,Imax.signals.values(:,1),'r');
xlim([0 3]);

subplot(313)
grid on
hold on
title('Pulsacja silnika');
xlabel('Czas [s]');
ylabel('Pulsacja [rad/s]');
plot(omega.time,omega.signals.values(:,1));
xlim([0 3]);
saveas(gcf, 'P_z_udarem.png');
