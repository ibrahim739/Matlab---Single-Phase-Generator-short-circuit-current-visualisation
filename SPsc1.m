%% This script reprents a Short circuit current calculation of an unloaded single phase generator.

% The Generator is modeled as a voltage source and a
% series RL branch which represents the internal parameters of the
% Generator

%% Brief explanation of the Single Phase short circuit calculation: 

% The total fault current also known as the Asymmetrical fault current is
% composed of two components: The Ac fault current also named as
% symmetrical current or steady-state fault current + the dc offset current
% which decays exponentially with time constant taw
% In this script, we assume Maximum Dc offset which represents the worst
% case.

%% The inputs to the function are:
% V which is the Rms voltage in Volts
% L is the internal inductance in Henry
% R is the internal resistance in ohm
% f is the operating frequency in Hertz

V = input('Please enter the Rms voltage in Volt: ');
L = input('Please enter the internal inductance in Henry: ');
R = input('Please enter the internal Resistance in ohm: ');
f = input('Please enter the operating frequency in Hertz: ');

%% Example (Typical input problems)
% A bolted short circuit occurs in a single phase generator with the maximum dc offset :
% V =20 kV, L = 0.0212 H, R = 0.8 Ω , F = 60 Hz.

% In this case, the inputs to the function are :

% V = 20000;
% L = 0.0212;
% R = 0.8;
% F = 60 ;
%--------------------------------------------------------------------------------------------------------

% Done by Ibrahim Al-Salloum (Electrical power and machines engineer)
% if you have any other specification that you want to add to this project,
% i would be glad to do it, send me a Message on WhatsApp on: +96176537146
% or via email: ibrahimmsalloum12@gmail.com

% -----------------------------------------------------------------------------------------------------
t = 0:0.01:3; % You can change the time limit 
              % number whenever you want, just change the 
              % time end value ( 3 in this
              % example)

w = 2*pi*f;
X = L*w;
taw = L/R;
Iac_rms = V/sqrt(X^2+R^2);

Ish = zeros(1,length(t));
Irms = zeros(1,length(t));
Idc = zeros(1,length(t));
for i= 1:length(t)
    Ish(i) = abs(sqrt(2) * Iac_rms * sind(w*t(i)-180/2) - sqrt(2) * Iac_rms * sind(-180/2) * exp(-t(i)/taw)) ;
    Idc(i) = sqrt(2) * Iac_rms * exp(-t(i)/taw);
    Irms(i) = sqrt(Iac_rms^2 + Idc(i)^2);
end
Iac_rms = Iac_rms *ones(1,length(t));
plot(t,Ish,'m')
hold on
plot(t,Irms,'k')
plot(t,Iac_rms,'r')
plot(t,Idc,'b')
pause(0.1)
yy = Irms(1) + Irms(1)/10;
xx = t(end) + 0.1;
axis([-0.08 xx 0 yy])
xlabel('Time in seconds')
ylabel('Short circuit current in Ampere')
title('Singe Phase Short circuit current')
legend('Iassymetrical(total)-instantaneous','Iassymetrical(total)-rms','Iac-rms','Idc')
grid on
hold off









    

