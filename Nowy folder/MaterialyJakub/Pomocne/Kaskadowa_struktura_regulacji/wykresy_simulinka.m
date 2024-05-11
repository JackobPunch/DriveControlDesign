figure(1)

close all
plot(out.tout, out.yout(:,1))
hold on
plot(out.tout, out.yout(:,2))
plot(out.tout, out.yout(:,4))
plot(out.tout, out.yout(:,3))
legend("Prad twornika", "napiecie", "predkosc", "prad zadany")

