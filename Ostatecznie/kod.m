close all,  clear,  clc;

%dane dla silnika nr.14
P_N = 17000; % [W] 
U_N = 230; % [V]
I_N = 85; % [A]
n_N = 700; % [obr/min]
w_N = (2*pi*n_N)/60; % [rad/s]
R_t = 0.253; % [om]
L_t = 0.0019; % [H]
J_s = 0.75; % [kgm^2]

% wyjściowe dane (punkt 4.1 z instrukcji)
lambda_N = 1.8;
p = 50;
Y = (10)/(2.5*I_N); % [V] (wg instrukcji ma być: (10*V) - w gotowcu brak V, jest liczone niżej) 
K_T = (10)/(1.2*w_N); % [V] (wg instrukcji ma być: (10*V) - w gotowcu brak V, jest liczone niżej) 
K_P = (1.5*U_N)/10; % [V]

% wyznaczenie transmitancji
    % niezbędne parametry
    psi_e = (U_N-R_t*I_N)/w_N;  % [Vs]
    T = L_t/R_t;    % [s]
    J = J_s + 10*J_s;   % [kgm^2]
    B = J*R_t/(psi_e^2);    %[s]
    M_n = I_N*psi_e;    % [Nm]
    

    % ograniczenia
    I_max = lambda_N*I_N; % ograniczenie wartości prądu
    dI_max = p*I_N; % ograniczenie pochodnej prądu (p-dopuszczalna krotność prądu znamionowego w czasie 1 sekundy)
    w_max = w_N; % ograniczenie wartości prędkości obrotowej
    J_max = (4*T*psi_e^2)/R_t; %maksymalna wartość J (czy jest sens?? bez tego i tak działa warunek dla 11 J_s)

G_wU_s = tf([1/psi_e],[B*T B 1]);
G_wM_s = tf([R_t*T/psi_e^2 R_t/psi_e^2],[B*T B 1]);
G_IU_s = tf([B/R_t 0],[B*T B 1]);
G_IM_s = tf([1/psi_e],[B*T B 1]);

%% Nastawy regulatora prądu 

T1 = 0.5*B*(1-sqrt(1-4*T/B));
B1 = B - T1;
beta = lambda_N/p;
K_z = (B1-beta)/(Y*B1); %zastępczy współczynnik wzmocnienia
m = T1;
V = K_P*Y*B*beta/(R_t*(B1-beta));
u_z0 = lambda_N*I_N*Y*B1/(B1-beta); %ograniczenie

G_RI_s = tf([m 1],[V 0]);

%parametry uzyskane wg kryterium symetrycznego

K_w = J/(2*beta*psi_e*K_T*K_z);
Tr = 4*beta;

G_F_s = tf([1],[Tr 1]);

% transmitancje (AUX)
GRw=K_w*tf([Tr,1],[Tr 0]);
GRI=tf([m 1],[V 0]);
%Moment obciążenia
MN=P_N/w_N;

sim('projekt.slx')

%prędkość
figure
plot(GwU.time, GwU.signals.values)
title('Odpowiedz skokowa transmitancji G_{wU}')
xlabel('t[s]')
ylabel('w(t) [rad/s]')

%prąd
figure
plot(GIU.time, GIU.signals.values)
hold on
plot(GIU.time,ones(1,length(GIU.time)).*I_N*lambda_N,'r-')
hold on
title('Odpowiedz skokowa transmitancji G_{IU}')
hold on
xlabel('t[s]')
ylabel('I_t(t) [A]')

%pochodna prądu
figure
plot(GdIU.time, GdIU.signals.values)
hold on
plot(GdIU.time,ones(1,length(GdIU.time)).*I_N*p,'r-')
hold on
title('Odpowiedz skokowa transmitancji G_{dIU}')
hold on
xlabel('t[s]')
ylabel('dI_t/dt [A]')
