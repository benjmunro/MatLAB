function [  x,xdot,xdotdot ] = Newmark_B_SDOF( m, k, d,F,t, dt, IC,varargin )
%-------------- Newmrk Beta integration - ENME402 ------------------------%
%
% Integrates a 1-DOF system with mass "m", spring stiffness "k" and damping
% coeffiecient "d", when subjected to an external load f(t).
% Returns the displacement, velocity and acceleration of the system with
% respect to an inertial frame of reference.
%
% Written by Ben Munro
% LAST MODIFIED:10/04/2013

%----------------------------- Input -------------------------------------%
%
%       m        = System Mass             [1,1]
%       k        = System Stiffness        [1,1]
%       d        = System Damping          [1,1]
%       f        = Externally Applied Load [n,1]
%       dt       = Time Step               {1,1]
%       t        = Totlal Time             [1,1]
%       IC       = Initial Conditions      [3,1]
%       varargin = Options (see below)
%       
%   
%---------------------------- Outputs -------------------------------------%
%
%       x         = Displacemente Response  [n,1]
%       xdot      = Velocity                [n,1]
%       xdotdot   = Acceleration            [n,1]
%
%  N = number of time steps
%
% The options include changing the value of the "gamma" and "beta"
% coefficient which appear in the formulation of the method. By default
% these values are set to gamma = 1/2 and alpha = 1/4.
% NOTE  - gamma = 1/2 and alpha = 1/4 is the constant acceleration method
%       - gamma = 1/2 and alpha = 1/6 is the linearly varying method
% EXAMPLE
% To change nemark's coefficients, say to gamma = 1/3 and alpha = 1/5, 
% the syntax is:
% [x, xdot, xdotdot] = newmark_AB( m, k, d,f,t, dt, IC, 1/3, 1/5 )

if nargin == 7
    disp('Using default values, Constant average acceleration method :');
    disp('    alpha = 1/2');
    disp('    beta  = 1/4');
    alpha = 1/2;
    beta = 1/4;
else
    if nargin == 9
        alpha = varargin{1};
        beta = varargin{2};
    disp('Using input values, linearly varying method (or your in trouble):');  
    else
        error('Incorrect number of imput arguments');
    end
end


%# Number of Steps
n =ceil( t/dt);  

%# Preallocate Array
x=zeros(1,n);
xdot=zeros(1,n);
xdotdot=zeros(1,n);


%# Specify Intial Conditions 
x(1)=IC(1);
xdot(1)=IC(2);
xdotdot(1)=IC(3);


%# Calculate Intergration Constants
b1 = 1 / ( beta * ( dt^2 ) );

b2 = 1 / ( beta * dt );

b3 = (1/( 2 * beta ) - 1 );

b4 = alpha / (beta * dt);

b5 = alpha / beta - 1; 

b6 = dt * (alpha / ( 2 * beta ) - 1);

b7 = m*b1 + d*b4 + k;

for i = 1:n-1

    % Solve for node dosplacement at time t
    x(i+1) = (F(i+1) + m*(b1*x(i) + b2*xdot(i) + b3*xdotdot(i))...
        + d*(b4*x(i) + b5*xdot(i) + b6*xdotdot(i))) / b7;

   
    % Calculate node velocitys and accelerations at time t
    xdot(i+1) = b4*(x(i+1) - x(i)) - b5*xdot(i) - b6*xdotdot(i);
    
    xdotdot(i+1) = b1*(x(i+1) - x(i)) - b2*xdot(i) - b3*xdotdot(i);

end

end