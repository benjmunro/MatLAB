%-------------------- ENGR401 - Unsteady Heat Transfer -------------------%
%
%
% Ben Munro
% Developed : 28/04/2013
% Last Modified : 28/04/2013
%
% V & M Example 8.1
% Discription : Solve for temperature along a thin plate wih tranisent
% 1D heat conduction,no internal source or heat loss and specified end
% temperatures.

clc
clear

%%%%##### WHAT METHOD ####%%%%
theta = 0; %# Fully explict
%# theta = 0.5; %# Crank Nickelson
%# theta = 1; %# Fully implicit

%# Known Values
W = 0.02;
dx = 0.004;
k = 10;
pc=10e6;
t = 120;
%dt = pc * (dx ^2)/(2*k);
dt=2;
Nt = t/dt;

%# Create Mesh
x = dx/2:dx:W;
N=length(x);

%# Inital Condition
T = zeros(N,Nt);
T(:,1) = 200;

%# Boundary Conditions
%# dT/dx = 0 at x = 0 and t > 0
%# T = 0 at x = L and t > 0
TB = 0;
%# Caculated Constants
aw = k / dx;
ae = k / dx;
ap0 = pc*(dx/dt);
ap = ap0 + theta*(aw + ae);

%# Create A Matrix
A = zeros(N,N);
b = zeros(N,1);

for i = 1:N
    
    if x(i) == min(x)                           %# Boundary Condition 1
        
        A(i,i) = ap0 + theta*ae;
        A(i,i+1) = -theta*ae;
        
    elseif x(i) == max(x)                       %# Boundary Point 2
        
        A(i,i) = ap0 + theta*(2*aw + aw);
        A(i,i-1) = -theta*aw;
        
    else                                        %# Internal Nodes
        
        A(i,i) = ap;
        A(i,i+1) = -theta*ae;
        A(i,i-1) = -theta*aw;
        
    end
    
end



for i = 2:Nt+1
    
    b(1) = ae*(1-theta)*T(2,i-1)+( ap0 - (1-theta)*ae)*T(1,i-1) ;
    
    b(N) = aw*(1-theta)*T(N-1,i-1)+2*aw*(theta*TB + (1-theta)*TB) ...
        +( ap0 - (1-theta)*2*ae -(1-theta)*aw )*T(N,i-1) ;
    
   
    for j = 2:N-1
        
        b(j) = (1-theta)*(aw*T(j-1,i-1) + ae*T(j+1,i-1)) + ...
            ( ap0 - (1-theta)*ae -(1-theta)*aw )*T(j,i-1);  
           
    end
   
      [L,U]=lu(A);
       y=L\b;           % Forward Substitution
      T(:,i)=U\y;      % Backward Substitution
    
end

plot(T)

