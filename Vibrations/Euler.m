function [ x, xdot, xdotdot ] = Euler( m, d, k, t, dt, IC, f )

%-------------------- Euler integration - ENME402 ------------------------%
% SINGLE degree of freedom system
% Written by Ben Munro
% LAST MODIFIED:10/04/2013
%------------------------------Inputs-------------------------------------%
%   m   = mass
%   d   = damping
%   k   = spring constant
%   dt  = time step
%   t   = total time
%   IC  = inital conditions 
%   f   = forcing function
%-----------------------------Outputs-------------------------------------%
%   x       = diplacement
%   xdot    = velocity
%   xdotdot = acceleration

n =ceil(t/dt);  %# number of steps

x=zeros(1,n);
xdot=zeros(1,n);
xdotdot=zeros(1,n);


%# Intial Conditions

x(1) = IC(1);
xdot(1) = IC(2);
xdotdot(1) = IC(3);

for i= 1:n
    
    x(i+1) = .5*xdotdot(i)*dt^2 + xdot(i)*dt + x(i);
    
    xdot(i+1) = xdotdot(i)*dt + xdot(i);
    
    xdotdot(i+1)= (1/m)*(f(i)-(d+k*dt)*xdot(i)-k*x(i)-(d*dt+(k*(dt)^2)/2)*xdotdot(i));
end


%# Plot option
%plot(0:dt:t,x)
%title('Euler integration, SDOF System');
%xlabel('Time (s)');
%ylabel('Displacement (m)');


end

