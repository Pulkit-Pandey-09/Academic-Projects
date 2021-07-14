function dydt = HH(t,y)
% Constants
ENa=55.17; % mv Na reversal potential
EK=-72.14; % mv K reversal potential
El=-49.42; % mv Leakage reversal potential

gbarNa=1.2; % mS/cm^2 Na conductance
gbarK=0.36; % mS/cm^2 K conductance
gbarl=0.003; % mS/cm^2 Leakage conductance

I = 0.1; %Applied Current

Cm = 0.01; %Membrane Capacitance

% Values set to equal input values
V = y(1);
n = y(2);
m = y(3);
h = y(4);

gNa=gbarNa*m^3*h;
gK=gbarK*n^4;
gl=gbarl;

INa=gNa*(V-ENa);
IK=gK*(V-EK);
Il=gl*(V-El);

%Hodgkin-Huxley Model Equation
dydt = [((1/Cm)*(I-(INa+IK+Il)));
         an(V)*(1-n)-bn(V)*n; 
         am(V)*(1-m)-bm(V)*m;
         ah(V)*(1-h)-bh(V)*h];