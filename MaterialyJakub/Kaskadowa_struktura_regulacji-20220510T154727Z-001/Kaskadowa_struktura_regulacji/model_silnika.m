%% Matematyczny model transmitancji silnika

% wypisane transmitancje ze strony 6
% transmitancje wypisuje w postaci symbolicznej, a nastepnie konwertuje do
% tf "control system"

B=J*Rt/(psi_en^2);
T=Lt/Rt;

% transmitancje, Mianownik wszedzie jest taki sam 
syms s

Gden = vpa(B*T*s^2+B*s+1, 3) % mianownik
GwUnum = vpa(1/psi_en, 3) % licznik 
GwU = syms2tf(GwUnum/Gden)

GwMnum = vpa(Rt*(T*s+1)/(psi_en^2), 3)
GwM = syms2tf(GwMnum/Gden)

GIUnum = vpa(B/Rt*s, 3) % licznik 
GIU = syms2tf(GIUnum/Gden)

GIMnum = vpa(1/psi_en, 3)
GIM = syms2tf(GIMnum/Gden)


%% Teraz te skoki jednostkowe
% Chyba chodzi o GwU i GwI, a pochodna prądu to GwI*s (s to idealny  człon
% różniczkujący)

% Uwaga: wyrażenia w latex'ie nie moga miec polskich znaków !!

close all

figure(1)
[y, t]=step(GIU*Un);
plot(t, y, t, In.*ones(1, length(t)), t, Imax.*ones(1, length(t)));
xlim([0, t(end)]);
title('Prad twornika $I$', 'Interpreter','latex');
legend('Prad twornika $I$', 'Prad znamionowy $I_N$', 'Ograniczenie pradu $I_{max}$','Interpreter','latex');

figure(2)
GIUderivate = syms2tf(GIUnum/Gden*s);
[y, t]=step(GIUderivate*Un);
plot(t, y, t, dIdtmax.*ones(1, length(t)));
xlim([0, t(end)]);
title('Pochodna pradu $\frac{dI}{dt}$', 'Interpreter','latex');
legend('Pochodna pradu $\frac{dI}{dt}$', 'Ograniczeni pochodnej pradu','Interpreter','latex');

figure(3)
[y, t]=step(GwU*Un); 
plot(t, y, t, omega_max.*ones(1, length(t)));
xlim([0, t(end)]);
title("Odpowiedz skokowa prędkości obrotowej");
legend("Prędkośc obrotowa \omega (t)", "Predkość znamionowa \omega_N");

