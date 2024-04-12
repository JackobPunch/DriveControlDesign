%% Komputerowe uk�ady sterowania
%
%% Projekt: Kaskadowa struktura regulacji nap�dem pr�du sta�ego
% Autor: Konrad Turbasa, Tomasz Kuder
% Data modyfikacji: 09.03.2016
% Kierunek: Elektrotechnika, Modu� C, Rok 3
% Termin zaj��: Wtorek godz. 9:30 - 11:00
%
%% Punkt 1
% Dane gr. 13

close all, clear all, clc

dane        %wczytanie danych projektowych
dane_sym    %wczytanie danych symulacji
warunek     %sprawdzenie warunk�w uk�adku

beta=lambdaN/p;
T1=0.5*B*(1-sqrt(1-4*T/B));
B1=B-T1;
Kz=(B1-beta)/(Y*B1);
m=T1;
V=beta*Y*Kp*B/((B1-beta)*Rt);
uz0=lambdaN*In*Y*B1/(B1-beta);

Kw=J/(2*Kt*Kz*beta*psien);
Tr=4*beta;

% Zad2
war=2;
switch war
% Rozruch bez momentu obci��enia
    case 1
        tm=0;
        Mn=0;
% Rozruch z momentem udarowym
    case 2
        tm=3;
        Mn=psien*In;
% Rozruch z momentem czynnym
    case 3
        tm=0;          
        Mn=psien*In;
    otherwise
        retun;
end
        
        deltaI=psien*V*Mn/(psien^2*V+J*Kp*Y);

% Zad 3
tau0=0.0033;
l1=[Kp*B/Rt];
m1=[V*B1*tau0 V*(tau0+B1) V+Kp*Y*B/Rt];
Gz=tf(l1,m1)


l2=[Kw*psien*Kp*B/Rt*Tr Kw*psien*Kp*B/Rt];
m2=[Tr*J*V*B1*tau0 Tr*J*V*(tau0+B1) Tr*J*V+Tr*J*Kp*Y*B/Rt 0 0];
Go=tf(l2,m2)
[GO,PO]=margin(l2,m2);
disp(['Zapas wzmocnienia: ' num2str(GO)]);
disp(['Zapas fazy: ' num2str(PO)]);

figure;
subplot(211)
bode(Gz)
title('Uk�ad zamkni�ty');
subplot(212)
bode(Go)
title('Uk�ad otwarty');
figure;
subplot(211)
nyquist(Gz)
title('Uk�ad zamkni�ty');
subplot(212)
nyquist(Go)
title('Uk�ad otwarty');


% Zad 4
% Uk�ad dyskretny
N=10;
Tp=beta/N;
Ti=Tr;
% Regulator pr�dko�ci
K2=Kw*(1-Tp/Ti);
K1=Kw;
% Regulator pr�du

q=2*N/(2^(16-1));
            