%-------------- Euler time stepping-SDOF Script - ENME402 ----------------%
%
% Written by Ben Munro
% LAST MODIFIED:10/04/2013

%# Clear Screen
clc
clear

%# Optional Load Data
%load('kobe.mat');

%# Enter constants
m = 10;
d = 4;
k =1700;

%# Initial Conditions
IC=[0;0;0];

%# Choose integrating Function
euler=@Euler;

%# Loop to Calculate response over various dt
for dt = 0.025:-0.005:0.001

[time,F]=interpclass(dt);


n=length(time);
t=max(time);

%# Execute program
[x,xdot,xdotdot]=euler( m, d, k, t, dt, IC, F );

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


