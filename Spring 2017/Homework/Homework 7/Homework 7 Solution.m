function kylelowryovereverything
for D=[140]

[t,y] = ode45(@alcohol_model,[0 24],[0; 0; 0;],[],D); %syntax for ODE for this problem
    figure(1)
    plot(t,(y(:,3)*.1),'-o')
    %ylim([0 5000])
    title('Blood acetaldehyde Content by Dose');
    xlabel('Time t');
    ylabel('blood acetaldehyde content');
    nameL=sprintf('dose %d',D);
    legend(nameL)



    figure(2)
    plot(t,(y(:,1)*.1),'-o')
    %ylim([0 5000])
    title('Blood alcohol Content by Dose');
    xlabel('Time t');
    ylabel('blood alcohol content');
    nameL=sprintf('dose %d',D);
    legend(nameL)
    hold on
end
end


function dydt = alcohol_model(t,y,D) 
%this is where we code in questions 21, 22
%
%define parameters
kemax = 10.2; %per hour
a = 0.00167; %1/g^2
ka=25.1; %per hour
V=44.100; %ml
Vm = 0.202;%mg/mlxhr
Km=0.0818;%mg/ml
VmaxB = .184;
VmaxA = .246;
KrevB = 1;
VrevB = 3.26;
KmB = 0.014; %mg/ml
KmA = 0.0000528; %mg/ml
if D==11.2
    F=0.785;
elseif D==22.4
    F=0.960;
else
    F=1;
        
end

ke=kemax/(1+a*D^2);

B = y(1); %the solver sends y = [y1,y2] we redefine these as P,Q
I = y(2); %Since these are what we actually have in the model
A = y(3);

dAdt = ((VmaxB*B - VrevB*A)/(KmB + B + KrevB*A)) - ((VmaxA*A)/(KmA+A)); % acetaldehyde in Body
dBdt = ka*I-  ((VmaxB*B - VrevB*A)/(KmB + B + KrevB*A)); % Alcohol in Body
dIdt = ke*(F*D/V)*exp(-ke*t)-ka*I; %Alcohol in Intestine
dydt = [dBdt; dIdt; dAdt;];

end

