%-------------------- Modal Domain Analysis - ENME402 --------------------%
%
%
% Written by Ben Munro
% Developed - 18/04/2013
% Last updated 19/04/2013

%# Clear screen, data
clc
clear
close all

%# Load Data, Optional
load('kobe.mat')

%# Known Variables
dt = 0.02;      %# Time step (s)
N = 3000;       %# Number of time steps
dof = 3;         %# Degrees of fredom of system

% Mass per story
m1=10000; m2=m1; m3=m2;

% Stiffness per story
k1=1600000; k2=k1; k3=k2;

% Damping per story
d1=13000; d2=d1; d3=d2;

% Mass Matrix
M = [m1 0 0; 0 m2 0; 0 0 m3];

% Damping Matrix
D = [d1+d2 -d2 0; -d2 d2+d3 -d3; 0 -d3 d3];

% Stiffness Matrix
K = [k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3];

%Force Matrix
Force(1:3,1:N)=-M*[1;1;1]*acc;

[X_hat, lambda] = eig(K, M);

%# plot mode shapes
figure(1)
hold on
plot([0; X_hat(:,1)],0:3,'s-m','MarkerSize',6);
plot([0; X_hat(:,2)],0:3,'p-b','MarkerSize',6);
plot([0; X_hat(:,3)],0:3,'o-r','MarkerSize',6);
%# Labels are erased, so generate them manually,add a legend in the upper left:
title('Mode Shapes and Natural Vibration Modes')
legend('Mode 1','Mode 2','Mode 3','Location','NorthEast')
hold off

%# Decoupling - Create modal coordinate EOM 
Mstar = eye(dof);
Kstar = lambda;
Dstar = X_hat'*D*X_hat;
Fstar = X_hat'*Force;

% Numerical Intergration
Nwmark=@Newmark_B_MDOF;     %Choose integrating Function
IC=[0 0 0];
[phi,phidot,phidotdot]=Nwmark( Mstar, Kstar, Dstar, Fstar, max(t), dt, IC, 3 );    % Execute program

%# Plot Modal Response
figure(2)
subplot(3,1,1); plot(t(1:2999),phi(1,:),'r')
 title('Mode 1 displacement');
 xlabel('Time (sec)');
 ylabel({'Modal Displacement' ;'(Scaled to mode shape)'});

 subplot(3,1,2); plot(t(1:2999),phi(2,:),'k')
 title('Mode 2 displacement');
 xlabel('Time (sec)');
 ylabel({'Modal Displacement' ;'(Scaled to mode shape)'});
 
 subplot(3,1,3); plot(t(1:2999),phi(3,:))
 title('Mode 3 displacement');
 xlabel('Time (sec)');
 ylabel({'Modal Displacement' ;'(Scaled to mode shape)'});
 
 %# Reconstruction using modal transform
 x = X_hat * phi;
 
 %# Plot Deformation in time domain
 figure(5)
 subplot(3,1,1); plot(t(1:2999),x(1,:),'r','linewidth',2)
 title('First Story');
 xlabel('Time (sec)');
 ylabel('Disp (m)');

 subplot(3,1,2); plot(t(1:2999),x(2,:),'k','linewidth',2)
 title('Second Story');
 xlabel('Time (sec)');
 ylabel('Disp (m)');
 
 subplot(3,1,3); plot(t(1:2999),x(3,:),'b','linewidth',2)
 title('Third Story');
 xlabel('Time (sec)');
 ylabel('Disp (m)');
