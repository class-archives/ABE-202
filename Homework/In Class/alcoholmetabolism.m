% Kathryn Atherton
% Barbara McAnulty

function main
   %solving eqtn 21 and 22 
   %[t,y] =ode45(function, time interval, initial conditions)
   %VDP is Van Der Pols ODE 
   D = [11.2, 22.4, 33.6, 45];
   Dparam = 11.2;
   [t, y] = ode45(@metmodel1, [0 15], [2;3], [],Dparam);
   plot(t, y)
   hold on 
   title('')
   xlabel('Time t')
   ylabel('Alcohol Conc y')
   legend('')
   Dparam = 22.4;
   [t, y] = ode45(@metmodel1, [0 15], [2;3], [],Dparam);
   plot(t, y)
   hold on 
   Dparam = 33.6;
   [t, y] = ode45(@metmodel1, [0 15], [2;3], [],Dparam);
   plot(t, y)
   hold on 
   Dparam = 45;
   [t, y] = ode45(@metmodel1, [0 15], [2;3], [],Dparam);
   plot(t, y)
   hold on 
end


function dydt = metmodel1(t, y, param)
    %y is a vector with two components
    %MATLAB file help has dydt = [y(2); 1-y*1(^2)*y(2)-y(1)];
    %code in eq 21 and 22 
    kemax = 10.2; %hr-1
    ka = 25.1; %hr-1
    a = 0.00167; %g-2
    V = 44.1; %L volume of water in bodu
    Vm = 0.202; %mg/mL
    Km = 0.0818; %mg/mL
    %Ffif= 0.785; %15 mL fraction of dose absorbed at lose dose
    %Fthir=0.960; %30 mL 
    I = y(1);
    B = y(2);
    D = param;
    ke = kemax/(1+a*D^2); %rate constant
    %if D<= 11.2
     %   F = Ffif;
      %  dIdt = (ke * (F*D/V)*exp(-ke*t)) - ka*I;
       % dBdt = ka*I - ((Vm*B) /(Km+B)); %what is I and B 
    %absorption rate into body minus liver's metabolism
%     elseif D>11.2 && D<=22.4
%         F = Fthir;
%         dIdt = ke * (F*D/V)*exp(-ke*t) - ka*I;
%         dBdt = ka*I - (Vm*B) /(Km+B);
%     else 
        F = 1;
        dIdt = ke * (F*D/V)*exp(-ke*t) - ka*I;
        dBdt = ka*I - (Vm*B) /(Km+B);
%     end
    dydt = [dIdt; dBdt];
end 