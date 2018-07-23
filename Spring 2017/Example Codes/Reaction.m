function main
%This is the matlab function to solve the 2 component gene network
% Need parameter, differential equation, and initial conditions.
% we will use Matlab ode45 built in ODE solver
 
%initial conditions 
Ca0 = 5;  %initial amount of p
Cb0 = 0;  %initial amount of q
Cc0=0;
T0=300;
Nw0=6.14;
time_begin = 0;
time_end = 360;

time = [time_begin, time_end];
initial_concentrations = [Ca0,Cb0,Cc0,T0,Nw0];
%Beginning of call to ODE solver
%Three steps- call the solver, tell it where the equations are to solve,
%capture the solution in T, Y

options = odeset('RelTol',1e-4);
[T,Y] = ode45(@rxn,time,initial_concentrations,options);

%plotting the solution
figure(1)
hold on
plot(T,Y(:,1),'-',T,Y(:,2),'.-',T,Y(:,3),'.-');
figure(2)
plot(T(:,1),Y(:,4));
hold on

end

function dy = rxn(t,y)
Ca =y(1);
Cb =y(2);
Cc =y(3);
T =y(4);
Nw =y(5);

v0=0.004;
Cb0=1;
UA=3000;
Ta=290;
cp = 75240;
T0=300;
dh=-7.9076e7;
Cw0=55;
k=0.39175*exp(5472.7*((1/273)-1/T));
Cd=Cc;
Vi=0.2;
Kc = 10^(3885.44/T);
cpa=170700;
V=Vi+v0*t;
Fb0 = Cb0*v0;
ra = -k*((Ca*Cb)-(Cc*Cd)/Kc);
Na = V*Ca;
Nb = V*Cb;
Nc = V*Cc;
rb = ra;
rc = -ra;
Nd = V*Cd;
NCp = cp*(Nb+Nc+Nd+Nw)+cpa*Na;

dCadt = ra-(v0*Ca)/V;
dCbdt = rb+(v0*(Cb0-Cb)/V);
dCcdt = rc-(Cc*v0)/V;
dTdt = (UA*(Ta-T)-Fb0*cp*(1+55)*(T-T0)+ra*V*dh)/NCp;
dNwdt = v0*Cw0;

dy = [dCadt;dCbdt;dCcdt;dTdt;dNwdt];

end