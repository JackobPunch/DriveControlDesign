%% Komputerowe uk³ady sterowania
% Projekt 2018  
% Autorzy: Tomasz Kielar, Tomasz Telesz  
% EAIiIB Elektrotechnika rok III

%% Czyszczenie przestrzeni roboczej
 clc, clear all, close all

%% Dane
Pn = 31000;  %moc znamionowa [W]
Un = 440;   %napiêcie znamionowe [V]
In = 77;    %pr¹d znamionowy [A]p
nn = 1400;  %prêdkoœæ znamionowa [obr/min]
R = 0.253;  %rezystancja twornika [Ohm]
L = 1.9*10^-3;   %indukcujnosc twornika [H]
Js = 0.75*1;  %mowladnosci [kg*m^2]
lambdan = 2;
p = 50; 

%% Obliczenia
Y = 10/2.5/In;   % wzmocnienie toru pomiarowego pr¹du
wn = nn*2*pi/60;    % pulsacja uk³adu
Kt = 10/1.2/wn;     % wzmocnienie toru pomiarowego prêdkoœci
Kp = 1.5*Un/10;     % wzmocnienie wzmacniacza mocy
T = L/R;    %elektromagnetyczna sta³a czasowa
J = 6*Js;   %zastêpszy moment bezw³adnoœci
psien = (Un-R*In)/wn;   %znamionowy strumieñ
B = J*R/psien^2;    %elektromechaniczna sta³a czasowa
Mn = psien*In;  %znamionowy moment obci¹¿enia
time = linspace(0,1,2);
Imax = In*lambdan+0*time;
omegan = wn;
tau = 3.3e-3;
Top = 0;
beta = lambdan/p;
T1 = 0.5*B*(1-sqrt(1-4*T/B));
B1 = B - T1;
kz = (B1 - beta)/(Y*B1); % zastêpczy wspó³czynnik wzmocnienia

% sprawdzenie warunku B > 4T:
if (B > 4*T)
    m = T1;                         %parametry regulatora
    V = beta*Y*Kp*B/((B1-beta)*R);  %
else
    m = sqrt(B*T);
    V = Kp*Y/R*B*beta/(sqrt(B*T)-beta);
end
Uz0 = lambdan*In*Y*B1/(B1-beta);    %ograniczenie Uz0 < 10V

%% Prze³¹cznia uk³adu
regpred = 1;        %1 - regulator PI, 0 - reg P
przelacz3 = 0;      %1 - za³¹cza moment obci¹¿enia
delomega = 0;       %1 - deltaomega 2%, 0 - 5%

TR = 4*beta;  %parametr wg kryterium symetrycznego
if (regpred == 0)
    Kw = Mn/(psien*kz*Kt*deltaomega);    %regulator P
else
    Kw = J/(2*Kt*kz*beta*psien);    %regulator PI
end;
%% Transmitancje

GomegaU = tf([1/psien],[B*T B 1])
GomegaM = tf([R*T/((psien)^2) R/((psien)^2)],[B*T B 1])
GIU = tf([B/R 0],[B*T B 1])
GIM = tf([1/psien],[B*T B 1])
%Gpsi = tf([],[])
% warunki pocz¹tkowe - u nas niekonieczne BO TAK

sim('Projekt_KUS_sim');
figure(1)
grid on
hold on
title('Odpowiedzi modelu GIU w postaci transmitancji na skok impulsowy');
xlabel('Czas [s]');
ylabel('Pr¹d [I]');
ylim([0 1800]);
plot(GUI_imp.time,GUI_imp.signals.values(:,1));
plot(time,Imax,'r--');

figure(2)
grid on
hold on
title('Odpowiedzi modelu GIU pochodna w postaci transmitancji na skok impulsowy');
xlabel('Czas [s]');
ylabel('Pochodna pr¹du [I/s]');
%ylim([-0.5, 28.e4]);
plot(GUI_der_imp.time,GUI_der_imp.signals.values(:,1));
plot(GUI_der_imp.time, p*In*ones(length(GUI_der_imp.time),1),'r--');
xlim([0 0.4]);


figure(3)
grid on
hold on
title('Odpowiedzi modelu GomegaM w postaci transmitancji na skok impulsowy');
xlabel('Czas [s]');
ylabel('Prêdkoœæ k¹towa [rad/s]');
%ylim([0 60]);
plot(GomegaU_imp.time,GomegaU_imp.signals.values(:,1));
plot(GomegaM_imp.time, wn*ones(length(GUI_der_imp.time),1),'r--');

%% Regulator P - rozruch bez momentu obci¹¿enia

regpred = 0;
przelacz3 = 0;
sim projekt;

figure(4)
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


%% Regulator P - rozruch z momentem udarowym

regpred = 0;
Top = 2;
sim projekt;

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

%% Regulator P - rozruch z momentem czynnym

regpred = 0;
przelacz3 = 1;
sim projekt;

figure()
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

%% Regulator P - rozruch z momentem biernym


%% Regulator PI - rozruch bez momentu obci¹¿enia

regpred = 1;
przelacz3 = 0;
sim projekt;

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
%% Regulator PI - rozruch z momentem udarowym

regpred = 1;
Top = 2;
sim projekt;

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
%% Regulator PI - rozruch z momentem czynnym

regpred = 1;
przelacz3 = 1;
sim projekt;

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
%% Regulator PI - rozruch z momentem biernym

%% Zapas modu³u i fazy oraz okreœliæ dopuszczalne opóŸnienie

KR=m/V;
tau0=3.3*10^(-3);
reg_pred=tf([KR*TR KR], [TR 0]);        %regulacja predkosci
reg_mom=tf([m 1], [V 0]);               %regulacja momentu
przeksztyr=tf([Kp], [tau0 1]);          %przeksztaltnik tyrystorowy
silnik=tf([B 0], [R*B1*T1 R*(B1+T1) R]);    %moment silnika
bezwl=tf([1],[J 0]);                    %bezwladnosc

sys=series(reg_mom,przeksztyr); 
sys=series(sys,silnik);

sys=feedback(sys,Y);

sys=series(psien,sys);
sys=series(reg_pred,sys);
sys=series(sys,bezwl);

P = nyquistoptions;
P.ShowFullContour = 'off'; 

figure(10)
subplot(211); grid on;
bode(sys)
subplot(212)
nyquist(sys, P)

figure(11)
margin(sys)

%% Dyskretny regulator prêdkoœci/pr¹du

N = 10; %liczba próbek
Tp = beta/N;    %krok próbkowania
K1 = Kw;    %wzmocnienia regulatora prêdkoœci
K2 = Kw*(1-(Tp/TR));

K1I = m/V;  %wzmocnienia regulatora pr¹du
K2I = (Tp-m)/V;

deltaI = psien*V*Mn/(psien^2*V+J*Kp*Y);

%% Porównanie przebiegów dla prawid³owo dobranego czasu próbkowania i zbyt du¿ego

figure()
for x = 1:2
    if(x == 2)
        Tp = 4*Tp;
    end;
    sim('Projekt_KUS_cyfro_sim');
    subplot(311)
    hold on
    grid on
    if(x == 1)
        plot(us.time,us.signals.values(:,1),'b')
        xlim([0 5]);
    elseif (x == 2)
        plot(us.time,us.signals.values(:,1),'k')
        xlim([0 5]);
    end;
    title('Napiêcie silnika');
    xlabel('Czas [s]');
    ylabel('Napiêcie [V]');
    if(x == 2)
        legend('Prawid³owo dobrany','Nieprawid³owo dobrany (za du¿y)');
    end;

    subplot(312)
    hold on
    grid on
    if(x == 1)
        plot(dltpipi.time,dltpipi.signals.values(:,1),'b')
        xlim([0 5]);
    elseif (x == 2)
        plot(dltpipi.time,dltpipi.signals.values(:,1),'k')
        xlim([0 5]);
    end;
    plot(Imax.time,Imax.signals.values(:,1),'r');
    title('Pr¹d silnika');
    xlabel('Czas [s]');
    ylabel('Pr¹d [A]');
    if(x == 2)
        legend('Prawid³owo dobrany', 'Pr¹d maksymalny', 'Nieprawid³owo dobrany (za du¿y)');
    end;

    subplot(313)
    hold on
    grid on
    if(x == 1)
        plot(domegapipi.time,domegapipi.signals.values(:,1),'b')
        xlim([0 5]);
    elseif (x == 2)
        plot(domegapipi.time,domegapipi.signals.values(:,1),'k')
        xlim([0 5]);
    end;
    title('Pulsacja silnika');
    xlabel('Czas [s]');
    ylabel('Pulsacja [rad/s]');
    if(x == 2)
         legend('Prawid³owo dobrany','Nieprawid³owo dobrany (za du¿y)');
    end;
end;