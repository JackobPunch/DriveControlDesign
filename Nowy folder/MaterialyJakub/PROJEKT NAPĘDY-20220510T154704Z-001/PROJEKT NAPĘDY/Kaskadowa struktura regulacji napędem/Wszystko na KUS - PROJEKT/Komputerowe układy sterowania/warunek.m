%% KUS projekt, wtorek 9.30



dane            % wczytanie danych do projektu



%% za�o�enia projektowe, wsp�czynniki w ograniczeniach

Jmr = 5*Js;     % moment bezw�adno�ci maszyny roboczej [ZA�O�ENIE PROJEKTOWE]

lambdaN = 2;     % pr�d dopuszczalny Id = 2*In

Id=lambdaN*In;

p = 50;         % dopuszczalna krotno�� pr�du z szybko�ci zmian pr�du



%dopuszczalna pr�dko�� silnika |w(t)| <= wd

T = Lt/Rt;              % elektromagnetyczna sta�a czasowa



wn = nn*2*pi/60;        % znamionowa pulsacja uk�adu

psien = (Un - Rt*In)/wn;  % <--- strumie� z charakterystyki mechanicznej

J = Js + Jmr;

B = J*Rt/psien^2;         % elektromechaniczna sta�a czasowa



if (B >= 4*T)

    disp('warunek ograniczenia pr�dko�ci obrotowej silnika zosta� spe�niony')    

else

    disp('warunek ograniczenia pr�dko�ci obrotowej silnika nie zosta� spe�niony')

end    