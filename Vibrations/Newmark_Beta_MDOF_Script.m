
% Newmark Beta MDOF script

load('Kobe.mat')

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


%Initial Conditions
IC=zeros(3,3);

%Time Values
t=60;
dt=0.02;
n=t/dt;

%Forcing Function
F=zeros(3,n+1);
F(:,2:n+1)=-M*[1;1;1]*acc;

%Choose integrating Function
Nwmark=@Newmark_B_MDOF;


% Execute program
[x,xdot,xdotdot]=Nwmark( M, K, D, F, t, dt, IC, 3 );
 
 
 %Time series for plot
 timeseries=0:dt:t;
 
 % Plot result
 subplot(3,1,1); plot(timeseries,x(1,:))
 title('First Story');
 xlabel('Time (sec)');
 ylabel('Disp (m)');

 subplot(3,1,2); plot(timeseries,x(2,:))
 title('Second Story');
 xlabel('Time (sec)');
 ylabel('Disp (m)');
 
 subplot(3,1,3); plot(timeseries,x(3,:))
 title('Third Story');
 xlabel('Time (sec)');
 ylabel('Disp (m)');
 