%% Stabilność kaskadowego układu regulacji %% statyzm 2%

close all
%% statyzm 2%
clear statyzm
statyzm=0.02*omega_n;
otwarty_model;

close all
figure(1)
% pełen Nyquist
nyquist(Go)
title("Wykres Nyquista - statyzm \Delta\omega=0.2\omega_N")
hold on
w=linspace(0, 2*pi, 100)
y=sin(w)
x=cos(w)
plot(y, x, 'r')
xlim([-2, 0])
ylim([-2, 2])

%Bode
handle=figure(2)
margin(Go)
[Gm,Pm,Wcg,Wcp] = margin(Go)
% Time Delay
text(12, -250, sprintf("Dopuszczalne opóźnienie \\tau_{max}=%.4fs", Pm*pi/180/Wcp))


% Phrase margin
figHandle=figure(3)
nyquist(Go)
hold on
w=linspace(0, 2*pi, 100)
y=sin(w)
x=cos(w)
plot(y, x, 'r')
xlim([-1.1, 0.1]);
ylim([-1.2, 1.2]);
plot([0 -0.0805*1.2], [0 -0.996*1.2], '-black')
addpath('circular_arrow/');
circular_arrow(figHandle, 0.5, [0 0], 180+Pm/2, Pm, 2, 'black', 5);
text(-0.4, -0.2, sprintf("P_{m}=%.2f^o", Pm))
title("Wykres Nyquista - statyzm \Delta\omega=0.2\omega_N")

% Gain margin
x_margin = -1/Gm
plot(x_margin, 0, 'o red')
plot([-1 -1], [-10 10], ':black')
plot([x_margin x_margin], [-10 10], ':black')
x = [0.22 0.82];
x2= flip(x)
y = [0.15 0.15];
    % adjust location of left arrow starting point (the sum of this with 'x' should not be negative)
annotation('textarrow',x,y, 'HeadLength', 3)
annotation('textarrow',x2,y, 'HeadLength', 3)
text(-0.9, -0.95, sprintf("G_{m}=%.2f=%.2fdB", Gm, 20*log10(Gm)))



