 function [F,Re, Im]=DFT(input,N,dt)

%-------------- Discrete Fourier Transform - ENME402 ---------------------%
%
% 
% Returns the real and imaginary parts of the F in the frequency domain
% 
%
% Written by Ben Munro
% LAST MODIFIED:11/04/2013

%----------------------------- Input -------------------------------------%
%
%       X        = Input function
%       N        = Number of Steps        [DOF,DOF]
%       dt       = Time step
%  
%----------------------------- Output ------------------------------------%
%
%       Re        = Real part of the DFT of F 
%       Im        = Imaginary part of the DFT of F
%
% Notes
% Will display the graph of Real and Imaginary against frequency 


df = 1 / (N*dt);

fmax = 1 / dt;

xk=zeros(N,1);

for i = 1:N
   
   k=i-1;
   sum=0;
   
   for j = 1:N
       
       n=j-1;
       term=2*pi*n*k/N;
       sum=sum+input(j)*(cos(term)-1i*sin(term));
       
   end 
   
   xk(i)=sum;
    
end
F=xk'/N;
Re=real(xk)/N;
Im=imag(xk)/N;

%subplot(2,1,1);plot(0:df:fmax-df,Re)
%title('Real');
%xlabel('Freqency (Hz)');
%ylabel('Real(Fn)');

%subplot(2,1,2); plot(0:df:fmax-df,Im)
%title('Imaginary');
%xlabel('Freqency (Hz)');
%ylabel('Imag(Fn) (m)');
plot(0:df:fmax-df,abs(F))
end