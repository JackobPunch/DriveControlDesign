%% Komputerowe uk�ady sterowania
% projekt1 2016
% Autorzy: 
%  
%%
%%
%%% 1. Czyszczenie przestrzeni roboczej
close all, clc, clear all
%%
%%
%%% 2. Dane
Pn = 43000;  %moc znamionowa [W]
Un = 230;   %napi�cie znamionowe [V]
In = 215;    %pr�d znamionowy [A]
Nn = 585;  %pr�dko�� znamionowa [obr/min]
R = 0.112;  %rezystancja twornika [Ohm]
L = 1.59*10^-3;   %indukcujnosc twornika [H]
Js = 2.1;  %moment bezwladnosci [kg*m^2]
%%
%%
%%% 3. Prze��cznia uk�adu
regpred = 1;        %1 - regulator PI, 0 - reg P
przelacz3 = 0;      %1 - za��cza moment obci��enia
delomega = 0;       %1 - deltaomega 2%, 0 - 5%

%%
%%
%%% 4. Sprawdzenie warunk�w
fn = Nn/60;     %cz�stotliwo�� znamionowa
omegan = 2*pi*fn;   %pulsacja znamionowa
lambdan = 2;
p = 50;
T = L/R;    %elektromagnetyczna sta�a czasowa
J = 7*Js;   %zast�pszy moment bezw�adno�ci
psien = (Un-R*In)/omegan;   %znamionowy strumie�
B = J*R/psien^2;    %elektromechaniczna sta�a czasowa
%%
%%
%%% 5. Obliczenia
tau0 = 3.3e-3;  %[s] sta�a czasowa przekszta�tnika tyrystorowego
Kp = 1.5*Un/10;    %wzmocnienie wzmacniacza mocy
Y = 10/(2.5*In);     %zmienione z 10 na 7.8, bo nie spe�nia�o warunku z Uz0
Kt = 10/(1.2*omegan);   %wzmocnienie
Mn = psien*In;  %znamionowy moment obci��enia
Top = 0;  %op�nienie za��czanie znamionowego momentu obci��enia
if (delomega == 0)
    deltaomega = 0.05*omegan;   %statyzm regulacji
else
    deltaomega = 0.02*omegan;
end;
beta = lambdan/p;   %sta�a czasowa przebiegu pr�du twornika
T1 = 0.5*B*(1-sqrt(1-4*T/B));
B1 = B-T1;
kz = (B1-beta)/(Y*B1);  %zast�pczy wsp�czynnik wzmocnienia
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
%%
%%
%%% 6. Odpowiedzi model�w w postaci transmitancji na skok impilsowy
sim('projekt1_KUS_sim1');
figure()
grid on
hold on
title('GIU');
xlabel('Czas [s]');
ylabel('Pr�d [I]');
%%%ylim([0 1500]);
plot(GUI_imp.time,GUI_imp.signals.values(:,1));
plot(time,Imax,'r');

figure()
grid on
hold on
title('GIU pochodna');
xlabel('Czas [s]');
ylabel('Pochodna pr�du [I/s]');
%%%ylim([-0.5, 14.e4]);
plot(GUI_der_imp.time,GUI_der_imp.signals.values(:,1));

figure()
grid on
hold on
title('GomegaM');
xlabel('Czas [s]');
ylabel('Pulsacja [rad/s]');
%%%ylim([0 4]);
plot(GomegaM_imp.time,GomegaM_imp.signals.values(:,1));

figure()
grid on
hold on
title('GomegaU');
xlabel('Czas [s]');
ylabel('Pulsacja [rad/s]');
%%%ylim([0 200]);
plot(GomegaU_imp.time,GomegaU_imp.signals.values(:,1));

figure()
grid on
hold on
title('GIM');
xlabel('Czas [s]');
ylabel('Pr�d [I]');
%%%ylim([0 80]);
plot(GIM_imp.time,GIM_imp.signals.values(:,1));
plot(time,Imax,'r');

%%
%%
%%% 7. Rysowanie przebieg�w symulacyjnych - uk�ad analogowy
for x = 1:2
    for p = 1:3
        if(x == 1)
            regpred = 1;    %regulator PI
        elseif(x == 2)
            regpred = 0;    %regulator P
        end;
        if(p == 1)
            przelacz3 = 0;   %wy��cza moment obci��enia
        elseif(p == 2)
            przelacz3 = 1;  %za��cza moment obci��enia
        elseif (p == 3)
            Top = 5;    %op�nie za��czenie momentu (za��czenie udarowe)
        end;
        sim('projekt1');

        figure()
        subplot(311)
        grid on
        hold on
        title('Napi�cie silnika');
        xlabel('Czas [s]');
        ylabel('Napi�cie [V]');
        plot(Us.time,Us.signals.values(:,1));
        %%%xlim([0 2]);

        subplot(312)
        grid on
        hold on
        title('Pr�d silnika');
        xlabel('Czas [s]');
        ylabel('Pr�d [A]');
        plot(Is.time,Is.signals.values(:,1));
        plot(Imax.time,Imax.signals.values(:,1),'r');
        %%%xlim([0 2]);

        subplot(313)
        grid on
        hold on
        title('Pulsacja silnika');
        xlabel('Czas [s]');
        ylabel('Pulsacja [rad/s]');
        plot(omega.time,omega.signals.values(:,1));
        %%%xlim([0 2]);

    end;
end;
%%
%%
%%% 8. Przekszta�tnik tyrystrowy
Gp = tf(Kp,[tau0 1]);
Gi = tf(m+1, V);
Gs = tf(B/R,[B1 1]);
G1 = Kt*tf([TR 1],[TR 0]);
Gj = tf(1,[J 1]);

Gz = Gi*Gp*Gs/((1+Gi*Gp*Gs)*Y);
G0 = G1*Gz*psien*Gj;

figure()
margin(G0);

%Charakterystyki
figure()   %dla Gz
subplot(211)
title('Charakterystyki dla Gz');
hold on
bode(Gz);
grid on
subplot(212)
nyquist(Gz);

% figure()   %dla G0
% subplot(211)
% title('Charakterystyki dla G0');
% hold on
% bode(G0);
% grid on
% subplot(212)
% nyquist(G0);
% 
%%
%%
%%% 9. Dyskretny regulator pr�dko�ci/pr�du

N = 10; %liczba pr�bek
Tp = beta/N;    %krok pr�bkowania
K1 = Kw;    %wzmocnienia regulatora pr�dko�ci
K2 = Kw*(1-(Tp/TR));

K1I = m/V;  %wzmocnienia regulatora pr�du
K2I = (Tp-m)/V;

deltaI = psien*V*Mn/(psien^2*V+J*Kp*Y);
%%
%%
%%% 10. Por�wnanie przebieg�w dla prawid�owo dobranego czasu pr�bkowania i zbyt du�ego

figure()
for x = 1:2
    if(x == 2)
        Tp = 4*Tp;
    end;
    sim('projekt1_KUS_cyfro_sim1');
    subplot(311)
    hold on
    grid on
    if(x == 1)
        plot(us.time,us.signals.values(:,1),'b')
        %%%xlim([0 6]);
    elseif (x == 2)
        plot(us.time,us.signals.values(:,1),'k')
        %%%xlim([0 6]);
    end;
    title('Napi�cie silnika');
    xlabel('Czas [s]');
    ylabel('Napi�cie [V]');
    if(x == 2)
        legend('prawod�owy','za duzy');
    end;

    subplot(312)
    hold on
    grid on
    if(x == 1)
        plot(dltpipi.time,dltpipi.signals.values(:,1),'b')
        %%%xlim([0 7]);
    elseif (x == 2)
        plot(dltpipi.time,dltpipi.signals.values(:,1),'k')
        %%%xlim([0 7]);
    end;
    plot(Imax.time,Imax.signals.values(:,1),'r');
    title('Pr�d silnika');
    xlabel('Czas [s]');
    ylabel('Pr�d [A]');
    if(x == 2)
        legend('prawod�owy','pr�d maksymalny','za duzy');
    end;

    subplot(313)
    hold on
    grid on
    if(x == 1)
        plot(domegapipi.time,domegapipi.signals.values(:,1),'b')
        %%%xlim([0 6]);
    elseif (x == 2)
        plot(domegapipi.time,domegapipi.signals.values(:,1),'k')
        %%%xlim([0 6]);
    end;
    title('Pulsacja silnika');
    xlabel('Czas [s]');
    ylabel('Pulsacja [rad/s]');
    if(x == 2)
        legend('prawod�owy','za duzy');
    end;
end;

%%
%%
%%% 11. Rysowanie przebieg�w symulacyjnych - uk�ad cyfrowy
Tp = beta/N;
for x = 1:2
    for p = 1:3
        if(x == 1)
            regpred = 1;    %regulator PI
        elseif(x == 2)
            regpred = 0;    %regulator P
        end;
        Top = 0;
        if(p == 1)
            przelacz3 = 0;   %wy��cza moment obci��enia
        elseif(p == 2)
            przelacz3 = 1;  %za��cza moment obci��enia
        elseif (p == 3)
            Top = 5;    %op�nie za��czenie momentu (za��czenie udarowe)
        end;
        sim('projekt1_KUS_cyfro_sim1');
    
        figure()
        subplot(311)
        hold on
        grid on
        plot(us.time,us.signals.values(:,1))
        %%%xlim([0 2]);
        title('Napi�cie silnika');
        xlabel('Czas [s]');
        ylabel('Napi�cie [V]');

        subplot(312)
        hold on
        grid on
        plot(dltpipi.time,dltpipi.signals.values(:,1))
        %%%xlim([0 2]);
        plot(Imax.time,Imax.signals.values(:,1),'r');
        title('Pr�d silnika');
        xlabel('Czas [s]');
        ylabel('Pr�d [A]');

        subplot(313)
        hold on
        grid on
        plot(domegapipi.time,domegapipi.signals.values(:,1))
        %%%xlim([0 2]);
        title('Pulsacja silnika');
        xlabel('Czas [s]');
        ylabel('Pulsacja [rad/s]');
    end;
end;
regpred = 1;
Top = 0;
% %%
% %%
% %%%Legenda do figure
% %%
% %figure(1) - Odpowiedzi modelu GIU w postaci transmitancji na skok impilsowy
% %figure(2) - Odpowiedzi modelu GIU pochodna w postaci transmitancji na skok impilsowy
% %figure(3) - Odpowiedzi modelu GomegaM w postaci transmitancji na skok impilsowy
% %figure(4) - Odpowiedzi modelu GomegaU w postaci transmitancji na skok impilsowy
% %figure(5) - Odpowiedzi modelu GIM w postaci transmitancji na skok impilsowy
% %figure(6) - Uk�ad analogowy - regulator PI - brak momentu obci�zenia
% %figure(7) - Uk�ad analogowy - regulator PI - moment obci�zenia
% %figure(8) - Uk�ad analogowy - regulator PI - moment obci�zenia w 5s
% %figure(9) - Uk�ad analogowy - regulator P - brak momentu obci�zenia
% %figure(10) - Uk�ad analogowy - regulator P - moment obci�zenia
% %figure(11) - Uk�ad analogowy - regulator P - moment obci�zenia w 5s
% %%figure(12) - Przekszta�tnik tyrostorowy - charakterystyki dla Gz
% %figure(12) - margine G0
% %figure(13) - Przekszta�tnik tyrostorowy - charakterystyki dla G0
% %figure(14) - Uk�ad cyfrowy - por�wnanie przebieg�w dla r�nych czas�w Tp
% %figure(15) - Uk�ad cyfrowy - regulator PI - brak momentu obci�zenia
% %figure(16) - Uk�ad cyfrowy - regulator PI - moment obci�zenia
% %figure(17) - Uk�ad cyfrowy - regulator PI - moment obci�zenia w 5s
% %figure(18) - Uk�ad cyfrowy - regulator P - brak momentu obci�zenia
% %figure(19) - Uk�ad cyfrowy - regulator P - moment obci�zenia
% %figure(20) - Uk�ad cyfrowy - regulator P - moment obci�zenia w 5s
% %%
% %%