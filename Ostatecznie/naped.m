%dane:
close all, clear all, clc
%nasze dane
Pn=17;
Un=230;
In=85;
nn=700;
R1=0.253;
L1=1.9*10^(-3);
Js=0.75;

%wyliczenie parametrów:

wn=(2*pi*nn)/60;
psie=(Un-R1*In)/wn;
J=Js*4;
B=J*(R1/psie^2);
T=L1/R1;
Mn=In*psie;

%przyjêto:

lambdan=1.8;
p=50;
Y=10/(2.5*In);
KT=10/(1.2*wn);
Kp=(1.5*Un)/10;

%za³ozenia zmiennych stanu:

Id=lambdan*In;

pIn=p*In;

%wd????

%PUNKT 6

%nastawy regulatora pr¹du


T1=0.5*B*(1-sqrt(1-(4*T/B)));

beta=lambdan/p;

B1=B-T1;

Kz=(B1-beta)/(Y*B1);

m=T1;

V=((Kp*Y)/R1)*((B*beta)/(B1-beta));

uz0=lambdan*In*((Y*B1)/(B1-beta))



%regulator predkosci P:

KwP=Mn/(psie*Kz*KT*6.078)



%regulator predkosci PI:

KwPI=J/(2*beta*psie*KT*Kz)

TrPI=4*beta
