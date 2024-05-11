clear; close all; clc;
%dane z tabelki
PN=54e3;
UN=440;
IN=130.7;
nN=1450; wN=nN*pi/30;
Rt=0.254;
Lt=1.63;
Js=0.97;

%obliczenia tych ograniczeń i innych
lambdaN=1.8;
p=50;
Y=10/(5*IN); %tak dobieramy ten mnożnik przed IN, aby uz0<10
KT=10/(1.2*wN);
Kp=1.5*UN/10;

fie=(UN-IN*Rt)/(wN);

T=0.001*Lt/Rt;
J=(Js+0*Js);
B=J*Rt/(fie^2);

%========================
%współczynniki (str. 8)

beta=lambdaN/p;

T1=0.5*B*(1-sqrt(1-4*T/B));
B1=B-T1;

%parametry reuglatora
m=T1;
V=beta*(Y*Kp*B)/((B1-beta)*Rt);

%zastępczy współczynnik wzmocnienia
kz=(B1-beta)/(Y*B1);

%ograniczenie
uz0=lambdaN*IN*(Y*B1)/(B1-beta)

%parametry uzyskane wg kryterium symetrycznego
TR=4*beta;
Kw=J/(2*KT*kz*beta*fie);

% transmitancje (AUX)
GRw=Kw*tf([TR,1],[TR 0]);
GRI=tf([m 1],[V 0]);
%Moment obciążenia
MN=PN/wN;

%prąd dodawany do lambda*IN podczas rozruchu napędu:
%deltaI=(fie*V*MN)/(fie^2*V+J*Kp*Y) 



