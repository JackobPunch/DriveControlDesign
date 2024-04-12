%% Stabilność kaskadowego układu regulacji PI

close all
Go=GwPI*Gc2*psi_en*Gs* Kt;

close all
figure(1)
% pełen Nyquist
nyquist(Go)
title("Wykres Nyquista - regulator PI")
hold on
w=linspace(0, 2*pi, 100)
y=sin(w)
x=cos(w)
plot(y, x, 'r')

%Bode
handle=figure(2)
margin(Go)
[Gm,Pm,Wcg,Wcp] = margin(Go)
% Time Delay
text(1.3, -250, sprintf("Dopuszczalne opóźnienie \\tau_{max}=%.4fs", Pm*pi/180/Wcp))


% Phrase margin
figHandle=figure(3)
nyquist(Go)
hold on
w=linspace(0, 2*pi, 100)
y=sin(w)
x=cos(w)
plot(y, x, 'r')
xlim([-1.1, 0.1]);
ylim([-1.5, 1.5]);
plot([0 -0.8418*1.2], [0 -0.539*1.2], '-black')
addpath('circular_arrow/');
circular_arrow(figHandle, 1, [0 0], 180+Pm/2, Pm, 2, 'black', 6);
text(-0.9, -0.1, sprintf("P_{m}=%.2f^o", Pm))
title("Wykres Nyquista - regulator PI")

% Gain margin
x_margin = -1/Gm
plot(x_margin, 0, 'o red')
plot([-1 -1], [-10 10], ':black')
plot([x_margin x_margin], [-10 10], ':black')
x = [0.22 0.62];
x2= flip(x)
y = [0.15 0.15];
    % adjust location of left arrow starting point (the sum of this with 'x' should not be negative)
annotation('textarrow',x,y, 'HeadLength', 4)
annotation('textarrow',x2,y, 'HeadLength', 4)
text(-0.9, -1.1, sprintf("G_{m}=%.2f=%.2fdB", Gm, 20*log10(Gm)))


