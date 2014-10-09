%-------------- Freqency Domain Analysis-MDOF - ENME402 ------------------%
%
%
% This script takes an input force, transfers it to the freqency domian 
% using fft. The script calculates the trasfer function which is then 
% mutiplied by the fft input force to get the response in the freqency
% domian which is then transfered back to the time domain using ifft.
%
%
% Written by Ben Munro
% Last updated 16/04/2013

%# Clear screen & data, close windows
clc
clear
close all

%# Opltional load data
load('Kobe.mat')

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

%# Calulated Variables
df = 1 / (N*dt);    %# Frequency step (Hz)    
fmax = 1 / dt;      %# Max Freqency (Hz)
Fnyst = fmax/2;     %# Nyquist Freqency (Hz)

%# Calculating Nyquist Frequency
if rem(N,2)==0
    Nnyst=N/2;
else 
    Nnyst=(N+1)/2;
end

%Forcing Function
Force=zeros(3,N);
Force(:,1:N)=-M*[1;1;1]*acc;

%# Fast Fourier Transform
F=fft(Force,[],2);

%# Calculate Transfer Function
omega=(0:df:Fnyst-df)*2*pi; %# Freqency range for transfer function (Rad/s)

%# Calculating deformation in freqency domain
X=zeros(dof,N);
for i=1:Nnyst 
    X(:,i)=(((K-omega(i)^2*M)+(1i*omega(i)*D)))\F(:,i);
end

%# Inverse Fourier transform
x=ifft(X,[],2,'symmetric');

%# Plot deformation
subplot(3,1,1); plot(t,x(1,:))
 title('First Story');
 xlabel('Time (sec)');
 ylabel('Disp (m)');

 subplot(3,1,2); plot(t,x(2,:))
 title('Second Story');
 xlabel('Time (sec)');
 ylabel('Disp (m)');
 
 subplot(3,1,3); plot(t,x(3,:))
 title('Third Story');
 xlabel('Time (sec)');
 ylabel('Disp (m)');
