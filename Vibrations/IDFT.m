 function [x, Re, Im]=IDFT(F,N,dt)

%--------------Inverse Discrete Fourier Transform - ENME402 ---------------------%
%
% 
% Returns the time domain form from the real and imaginary parts of the F 
% in the frequency domain
%
% Written by Ben Munro
% LAST MODIFIED:14/04/2013

%----------------------------- Input -------------------------------------%
%
%       F        = Input function
%       N        = Number of Steps        [DOF,DOF]
%       dt       = Time step
%  
%----------------------------- Output ------------------------------------%
%
%       Re        = Real part of the IDFT of X 
%       Im        = Imaginary part of the IDFT of X
%
% Notes
% Will display the graph of Real and Imaginary formas in t domain 


df = 1 / (N*dt);

fmax = 1 / dt;

xk=zeros(N,1);

for i = 1:N
   
   k=i-1;
   sum=0;
   
   for j = 1:N
       
       n=j-1;
       term=2*pi*n*k/N;
       sum=sum+F(j)*(cos(term)+1i*sin(term));
       
   end 
   
   xk(i)=sum;
    
end
x=xk;
Re=real(xk)/N;
Im=imag(xk)/N;

%subplot(2,1,1);plot(0:dt:(N-1)*dt,Re)
%title('Real');
%xlabel('Freqency (Hz)');
%ylabel('Real(Fn)');

%subplot(2,1,2); plot(0:dt:(N-1)*dt,Im)
%title('Imaginary');
%xlabel('Freqency (Hz)');
%ylabel('Imag(Fn) (m)');

end