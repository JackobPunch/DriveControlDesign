% BADAM TRANSMITANCJE samego przekształtnika a ma byc cały model !!
% do poprawy, zostawiam na wzór

Gp=tf(Kp, [tau0, 1])

[Gm,Pm,Wcg,Wcp] = margin(Gp)
% BRAK ZAPASU AMPLITUDY

close 
figure(2)
bode(Gp)

% linie zapas amplitudy i fazy
% z tym że nie ma zapasu amplitudy, ch-ka fazowa nie przechodzi przez -180
% stopni
hold on
plot([Wcp Wcp],[-1e5 1e5], ':')
xlim([10, 5*10^5])
ylim([-100, 0])
% dodanie strzałek
% https://www.mathworks.com/matlabcentral/answers/160487-how-can-i-draw-a-line-with-arrow-head-between-2-data-points-in-a-plot
x = [0.7 0.64];
y = [0.2 0.15];
annotation('textarrow',x,y,'String','Zapas fazy')

%("test");