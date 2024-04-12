%% KUS projekt, wtorek 9.30



dane            % wczytanie danych do projektu



%% za³o¿enia projektowe, wspó³czynniki w ograniczeniach

Jmr = 5*Js;     % moment bezw³adnoœci maszyny roboczej [ZA£O¯ENIE PROJEKTOWE]

lambdaN = 2;     % pr¹d dopuszczalny Id = 2*In

Id=lambdaN*In;

p = 50;         % dopuszczalna krotnoœæ pr¹du z szybkoœci zmian pr¹du



%dopuszczalna prêdkoœæ silnika |w(t)| <= wd

T = Lt/Rt;              % elektromagnetyczna sta³a czasowa



wn = nn*2*pi/60;        % znamionowa pulsacja uk³adu

psien = (Un - Rt*In)/wn;  % <--- strumieñ z charakterystyki mechanicznej

J = Js + Jmr;

B = J*Rt/psien^2;         % elektromechaniczna sta³a czasowa



if (B >= 4*T)

    disp('warunek ograniczenia prêdkoœci obrotowej silnika zosta³ spe³niony')    

else

    disp('warunek ograniczenia prêdkoœci obrotowej silnika nie zosta³ spe³niony')

end    