%% Cały model układu kaskadowej regulacji

regulatory

% Gf - filtr wejściowy
% GwP - regulator predkości P
% Gc - regulator prądu + przekształtnik
% Gc2 - regulator prądu + przekształtnik (dokładny)
% GI=tf([1/Rt], [T 1]) - transmitancja silnika GI=Itwornika\Uwejsciowe
Gs=tf([1], [J 0]) % - silnik
% psi_en - strumień znamionowy
% Y - wzmocnienie toru pomiarowego prądu
% Kt - wzmocnienie toru pomiarowego prędkości

% Transmitancje opisuje dokunując odpowiednich przekształceń
% robie to na kartce skan jest w pliku: 

Go=GwP*Gc2*psi_en*Gs* Kt;
