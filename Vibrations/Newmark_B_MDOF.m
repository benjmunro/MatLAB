
function [  x,xdot,xdotdot ] = Newmark_B_MDOF( M, K, D, F, t, dt, IC,DOF ,varargin )
%-------------- Newmrk Beta integration - ENME402 ------------------------%
%
% Integrates a M-DOF system with mass Matrix "M", stiffness Matrix "K" and 
% damping coeffiecient Matrix "D", when subjected to an external load f(t).
% Returns the displacement, velocity and acceleration of the system with
% respect to an inertial frame of reference.
%
% Written by Ben Munro
% LAST MODIFIED:10/04/2013

%----------------------------- Input -------------------------------------%
%
%       M        = System Mass             [DOF,DOF]
%       K        = System Stiffness        [DOF,DOF]
%       D        = System Damping          [DOF,DOF]
%       F        = Externally Applied Load [n,1]
%       dt       = Time Step               {1,1]
%       t        = Totlal Time             [1,1]
%       IC       = Initial Conditions      [3,1]
%       varargin = Options (see below)
%       
%   
%---------------------------- Outputs -------------------------------------%
%
%       x         = Displacemente Response  [n,DOF]
%       xdot      = Velocity                [n,DOF]
%       xdotdot   = Acceleration            [n,DOF]
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
% [x, xdot, xdotdot] = Newmark_B_MDOF( ( M, K, D, F, t, dt, IC,DOF, 1/3, 1/5 )


if nargin == 8
    disp('Using default values, Constant average acceleration method :');
    disp('    alpha = 1/2');
    disp('    beta  = 1/4');
    alpha = 1/2;
    beta = 1/4;
else
    if nargin == 10
        alpha = varargin{1};
        beta = varargin{2};
        
    else
        error('Incorrect number of imput arguments');
    end
end

% Number of Steps
n = t/dt;  

% Preallocate Array
x=zeros(DOF,n);
xdot=zeros(DOF,n);
xdotdot=zeros(DOF,n);

% Specify Intial Conditions 
x(:,1)=IC(:,1);
xdot(:,1)=IC(:,2);
xdotdot(:,1)=IC(:,3);

% Calculate Intergration Constants
b1 = 1 / ( beta * ( dt^2 ) );

b2 = 1 / ( beta * dt );

b3 = (1/( 2 * beta ) - 1 );

b4 = alpha / (beta * dt);

b5 = alpha / beta - 1; 

b6 = dt * (alpha / ( 2 * beta ) - 1);

Khat = M*b1 + D*b4 + K;
 
% LU Decomposition
[L,U,P]=lu(Khat,'vector');

for i = 1:n-1

    % Solve for node displacement at time t
    Fhat = (F(:,i+1) + M*(b1*x(:,i) + b2*xdot(:,i) + b3*xdotdot(:,i))...
        + D*(b4*x(:,i) + b5*xdot(:,i) + b6*xdotdot(:,i))) ;

 
    x(:,i+1) = U\(L\(Fhat(P,:))); 
   
    % Calculate node velocitys and accelerations at time t
    xdot(:,i+1) = b4*(x(:,i+1) - x(:,i)) - b5*xdot(:,i) - b6*xdotdot(:,i);
    
    xdotdot(:,i+1) = b1*(x(:,i+1) - x(:,i)) - b2*xdot(:,i) - b3*xdotdot(:,i);

end

