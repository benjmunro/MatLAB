%------------------- Euler-SDOF Script - ENME402 ------------------------%
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

%# Time Values
%t=60;
dt=0.001;

%# Initial Conditions
IC=[0;0;0];

%#Forcing Function
[time,F]=interpclass(dt);

%# Choose integrating Function
euler=@Euler;

%# Calculated values
n=length(time);
t=max(time);

%# Execute program
[x,xdot,xdotdot]=euler( m, d, k, t, dt, IC, F );

%# Calculate the maximum displacement and display
maxdisp=max(abs(x));
fprintf('the maximum absolute displacement time step %4.3f s  is %4.4f m \n  ',dt, maxdisp)

%# Plot result

plot(time,x)
title('Newmark Beta SDOF System');
xlabel('Time (s)');
ylabel('Displacement (m)');
str = fprintf('%s',dt);
legend('str')

 
figure(1)
