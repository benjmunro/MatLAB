clc, clear

%# known constants

load('Christchurch.mat')

k1 	= 12000;        %kN/m	
k2 	= 180000;       %kN/m
c1 	= 1800;      %Ns/m	
c2 	= 5500;      %Ns/m
m1	= 10000; 	 %kg
m2	= m1;

N=length(t);
fmax=1/dt;
fnyst=fmax/2;
df=fmax/N;

%# Calculating Nyquist Frequency
if rem(N,2)==0
    Nnyst=N/2;
else
    Nnyst=(N+1)/2;
end

omeg = (0:df:fnyst)*2*pi;
omeghz=(0:df:fnyst);

input=fft(disp');

H1=((c1*1i*omeg)+k1)./((-m1*omeg.^2)+(c1*1i*omeg)+k1 );

H2=(c2*omeg*1i+k2)./(-m2*omeg.^2+c2*omeg*1i+k2 );

%H=(((k-(omega.^2.*m))+(1i.*omega.*d)).^(-1));

%# Calculating deformation in freqency domain
X =H1.*(input(1,1:Nnyst));

%#Padding
X_vector= padarray(X,[0 Nnyst],'post');

%# Inverse Fourier transform
x=ifft(X_vector,'symmetric');

%# Plot deformation
plot(t,real(x(:,1:4329)))
title('Frequency Domain Analysis SDOF System');
xlabel('Time (s)');
ylabel('Displacement (m)');