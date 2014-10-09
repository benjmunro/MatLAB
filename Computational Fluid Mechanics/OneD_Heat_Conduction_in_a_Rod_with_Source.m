%----------- ENGR401 - 1D Heat Conduction in a Rod with Source -----------%
%
%
% Ben Munro
% Developed : 27/04/2013
% Last Modified : 28/04/2013
%
%V & M Example 4.2
% Discription : Solve for temperature along a rod wih 1D heat conduction,
% with source, no heat loss and specified end temperatures.

clc
clear

%# Known Values
L = 0.02;
dx = 0.004;
k = 0.5;
q = 1000000;
TA = 100;
TB = 200;

%# Create Mesh
x = dx/2:dx:L;
N=length(x);

%# Calculated Constants
Area = 1;
aw = ((k/dx)*Area);
ae = ((k/dx)*Area);

%# Create Matrix
A = zeros(N,N);
b = zeros(N,1);

for i = 1:N
    
    if x(i) == min(x)                           %# Boundary Point 1
        
        A(i,i) = +ae + ((k*2/dx)*Area);
        A(i,i+1) = -ae;
        b(i,1) = ((k*2/dx)*Area) * TA + q*Area*dx;
        
    elseif x(i) == max(x)                       %# Boundary Point 2
        
        A(i,i) = aw + ((k*2/dx)*Area);
        A(i,i-1) = -aw;
        b(i,1) = ((k*2/dx)*Area) * TB + q*Area*dx;
        
    else                                        %# Internal Nodes 
        
        A(i,i) = ae + aw;
        A(i,i+1) = -ae;
        A(i,i-1) = -aw;
        b(i,1) = q*Area*dx;
        
    end
    
end

%# Solve linear system of equations For Temp along bar 
T = A\b;

disp('Temperature along rod')
disp(T)

%# Plot results
plot([0,x,L],[TA;T;TB], 'x-')
title('1D Heat Conduction in a Rod')
xlabel('meters (m)')
ylabel('Temperature (C)')


