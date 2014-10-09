%---------------- ENGR401 - 1D Heat Conduction in a Rod ------------------%
%
%
% Ben Munro
% Developed : 27/04/2013
% Last Modified : 28/04/2013
%
%V & M Example 4.1
% Discription : Solve for temperature along a rod wih 1D heat conduction,
% no internal source or heat loss and specified end temperatures.

clc
clear

%# Known Values
L = 0.5;
dx = 0.1;
k = 1000;
TA = 100;
TB = 500;

%# Create Mesh
x = dx/2:dx:L;
N=length(x);

%# Calculated Constants
Area = dx * dx;
aw = ((k/dx)*Area);
ae = ((k/dx)*Area);

%# Create Matrix
A = zeros(N,N);
b = zeros(N,1);

for i = 1:N
    
    if x(i) == min(x)                           %# Boundary Point 1
        
        A(i,i) = +ae + ((k*2/dx)*Area);
        A(i,i+1) = -ae;
        b(i,1) = ((k*2/dx)*Area) * TA;
        
    elseif x(i) == max(x)                       %# Boundary Point 2
        
        A(i,i) = aw + ((k*2/dx)*Area);
        A(i,i-1) = -aw;
        b(i,1) = ((k*2/dx)*Area) * TB;
        
    else                                        %# Internal Nodes 
        
        A(i,i) = ae + aw;
        A(i,i+1) = -ae;
        A(i,i-1) = -aw;
        
    end
    
end

%# Solve linear system of equations For Temp along bar 
T = A\b;

disp('Temperature along rod')
disp(T)

%# Plot results
plot([0,x,L],[TA;T;TB])
title('1D Heat Conduction in a Rod')
xlabel('meters (m)')
ylabel('Temperature (C)')


