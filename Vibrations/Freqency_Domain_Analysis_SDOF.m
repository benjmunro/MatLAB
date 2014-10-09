%-------------- Freqency Domain Analysis-SDOF - ENME402 ------------------%
%
% This script takes an input force, transfers it to the freqency domian 
% using fft. The script calculates the trasfer function which is then 
% mutiplied by the fft input force to get the response in the freqency
% domian which is then transfered back to the time domain using ifft.
%
% Written by Ben Munro
% Last updated 16/04/2013


%# Clear screen, data
clc
clear

%# Load Data, Optional
load('kobe.mat')

%# Known Variables
m = 10000;      %# Mass (kg)
k = 1600000;    %# Stiffness Coeficient (N/m)
d = 5060;       %# Damping Coefficent (Ns/m)
s = 0.02;       %# Damping
dt = 0.02;      %# Time step (s)
N = 3000;       %# Number of time steps

%# Calulated Variable
df = 1 / (N*dt);    %# Frequency step (Hz)    
fmax = 1 / dt;      %# Max Freqency (Hz)
Fnyst = fmax/2;     %# Nyquist Freqency (Hz)

%# Calculating Nyquist Frequency
if rem(N,2)==0
    Nnyst=N/2;
else
    Nnyst=(N+1)/2;
end

%# Force Vector
Force=-acc*m;

%# Fast Fourier Transform
F=fft(Force);

%# Calculate Transfer Function
omega=(0:df:Fnyst-df)*2*pi; %# Freqency range for transfer function (Rad/s)

H=(((k-(omega.^2.*m))+(1i.*omega.*d)).^(-1));

%# Calculating deformation in freqency domain
X =H.*(F(1,1:Nnyst));

%#Padding
X_vector= padarray(X,[0 Nnyst],'post');

%# Inverse Fourier transform
x=ifft(X_vector,'symmetric');

%# Plot deformation
plot(t,real(x))
title('Frequency Domain Analysis SDOF System');
xlabel('Time (s)');
ylabel('Displacement (m)');
