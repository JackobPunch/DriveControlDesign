%% Regulatory ciagłe
% Ważne: ze wzgledu na kolizje oznaczeń zakomentować to co nie "idzie" do
% simulinka

%% Parametry silnika i przekształtnika tyrystorowqego
B=J*Rt/(psi_en^2);
T=Lt/Rt;
beta=2*tau0;
kz=1/Y;

%% Przekształtnik tyrystorowy
Gp=tf(Kp, [tau0, 1])


%% Parametry regulatora prądu PI
KRi = T*Rt/(2*Kp*Y*tau0);
TRi = T;
Gpradu=tf([KRi*TRi KRi], [TRi 0])


%% Parametry regulatora predkości P
KomegaP=Mn/(psi_en*kz*Kt*statyzm);
GwP=tf([KomegaP], [1])

%% Parametry regulatora predkości PI
Komega=J/(2*Kt*kz*beta*psi_en)
Kw=Komega;
Tr_pred=4*beta;
GwPI=tf([Komega*Tr_pred Komega], [Tr_pred 0])

%% Filtr wejściowy
Gf=tf([1], [4*beta 1])

%% Sprawdzenie transmitancji regulatrora prądu
% wyznaczam Greg_pradu na podstawie Gp, Gpradu oraz GI
% nastepnie sprawdzam czy jest to w przybliżeniu równe Gc
 GI=tf([1/Rt], [T 1])
 
Greg_pradu=feedback(Gpradu*Gp*GI, Y)
Gc=tf([1/Y], [2*tau0 1])
Gc2=tf([1/Y], [2*tau0^2 2*tau0 1])
% bode(Greg_pradu)
% hold on
% bode(Gc, 'r')
% bode(Gc2, '.')
% 
% legend("G_{reg pradu}-na podstawie pełnego modelu str 161 'automaty 2020'", "Gc przybliżone", "Gc - pełne")

% Wnioski
% Dla transmitancji Gc dokładnej transmitancje są takie same, a wiec 
% obliczenia sa poprawne
% Dla transmitancji przybliżonej Gc podstawowy biegun jest podobny
