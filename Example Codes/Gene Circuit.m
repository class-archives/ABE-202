

%%
function main
% MatLab code to solve for two component gene network
%
% Chen et al 2002 IEEE: 
% A model of periodic oscillation for genetic regulatory systems
%
% Solves equations (21) and (22)



% outputs time and y, requires function, time interval and initial
% condition

[t,y]=ode45(@gene_model, [0 15], [2;3]);

subplot(2,1,1);
plot(t,y(:,1), '-o',t,y(:,2),'-o');
title('Gene Expression Over Time');
xlabel('Time t');
ylabel('P,Q');
legend('y1', 'y2');


subplot(2,1,2);
plot(y(:,2), y(:,1), '-o');
title('P vs. Q Concentration');
ylabel('P');
xlabel('Q');
hold on;

end

function dydt=gene_model(t,y)
%This is where we code in 21, 22
%y is concentrations

kp=1;
kq=1;
k1=15;
k2=0.2;
k3=0.1;
k4=10;
epsilon=0.11;

P=y(1); % The solver sends y=[y1,y2], redefine as P,Q
Q=y(2); % Since these are what we actually have in the model

dPdt=-kp*P+k1/(Q+k2);
dQdt=(1/epsilon)*(-kq*Q+Q^2/(Q^2+k4)*P+k3);

dydt=[dPdt; dQdt];

end



