%% Alcohol Metabolism

%%
function main




% outputs time and y, requires function, time interval and initial
% condition

[t,y]=ode45(@metabolism, [0 7], [0;0]);

plot(t,y(:,2), '-o');
title('Alcohol Metabolism Over Time');
xlabel('Time t');
ylabel('B');
legend('y2');
hold on;


end

function dydt=metabolism(t,y)
%This is where we code in 21, 22
%y is concentrations

keMax=10.2;
a=0.00167;
D=11.2;
%11.2, 22.4, 33.6, 45
F=0.785;
%0.785,0.96,1,1
V=48;
Vm=0.202;
ka=25.1;
Km=0.0818;

ke=keMax/(1+a*D^2);

I=y(1); % The solver sends y=[y1,y2], redefine as P,Q
B=y(2); % Since these are what we actually have in the model

dIdt=ke*(F*D/V)*exp(-ke*t)-ka*I;
dBdt=ka*I-Vm*B/(Km+B);

dydt=[dIdt; dBdt];

end