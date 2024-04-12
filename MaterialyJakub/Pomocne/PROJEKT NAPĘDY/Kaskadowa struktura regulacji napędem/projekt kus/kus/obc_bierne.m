figure; close all;
sim('model_moment_bierny')
subplot(311)
grid on
hold on
title('Napiêcie silnika');
xlabel('Czas [s]');
ylabel('Napiêcie [V]');
plot(wyniczki_bierne.time,wyniczki_bierne.signals.values(:,1));
xlim([0 4]);

subplot(312)
grid on
hold on
title('Pr¹d silnika');
xlabel('Czas [s]');
ylabel('Pr¹d [A]');
plot(Imax.time,Imax.signals.values(:,1),'r');
plot(wyniczki_bierne.time,wyniczki_bierne.signals.values(:,2));
xlim([0 4]);

subplot(313)
grid on
hold on
title('Pulsacja silnika');
xlabel('Czas [s]');
ylabel('Pulsacja [rad/s]');
plot(wyniczki_bierne.time,wyniczki_bierne.signals.values(:,4));
xlim([0 4]);
%saveas(gcf, 'PI_bez_mom.png');