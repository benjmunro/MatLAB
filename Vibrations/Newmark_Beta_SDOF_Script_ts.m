% Newmark Beta SDOF script
clc
clear

%load('kobe.mat');
%constants
m = 10;
d = 4;
k = 1700;

%Initial Conditions
IC=[0;0;0];

%# Loop to Calculate response over various dt
for dt = 0.025:-0.005:0.001

%Forcing Function
[time,F]=interpclass(dt);
n=length(time);
t=max(time);
%Choose integrating Function
Nwmark=@Newmark_B_SDOF;

% Execute program
[x,xdot,xdotdot]=Nwmark( m, k, d, F, t, dt, IC );

%# Calculate the maximum displacement and display
maxdisp=max(abs(x));
fprintf('the maximum absolute displacement time step %4.3f s  is %4.4f m \n  ',dt, maxdisp)

hold all
%plot result
plot(time,x)
title('Newmark Beta SDOF System');
xlabel('Time (s)');
ylabel('Displacement (m)');
hold off
figure(1)
end

