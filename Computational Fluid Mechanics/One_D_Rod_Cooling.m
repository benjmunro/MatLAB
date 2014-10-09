%------------------ ENGR401 - 1D cooling in a Rod  -----------------------%
%
%
% Ben Munro
% Developed : 29/04/2013
% Last Modified : 29/04/2013
%
%V & M Example 4.3
% Discription : Cooling of a cirular fin by convection along its length

clc
clear

%# Known Values
L = .15;
dx = .015;
k = 398;
h = 100;
TA = 25;
TB = 100;
Area = .005^2*pi/4;
p = pi * 5/1000;

n= h*p/k/Area;

%# Create Mesh
x = dx/2:dx:L;
N=length(x);

%# Calculated Constants
Area = 1;
aw = 1/dx;
ae = 1/dx;

%# Create Matrix
A = zeros(N,N);
b = zeros(N,1);

for i = 1:N
    
    if x(i) == min(x)                           %# Boundary Point 1
        
        A(i,i) = ae + (n * dx + (2/dx));
        A(i,i+1) = -ae;
        b(i,1) = (n * dx *TA) + (2/dx)*TB;
        
    elseif x(i) == max(x)                       %# Boundary Point 2
        
        A(i,i) = ae+  n *dx;
        A(i,i-1) = -aw;
        b(i,1) = (n * dx *TA); 
        
    else                                        %# Internal Nodes 
        
        A(i,i) = ae + aw + n * dx;
        A(i,i+1) = -ae;
        A(i,i-1) = -aw;
        b(i,1) = (n * dx *TA) ;
        
    end
    
end

%# Solve linear system of equations For Temp along fin 
T = A\b;

disp('Temperature along rod')
disp(T)

%# Plot results
plot([0,x],[TB;T], 'o-','markersize',10)
title('Cooling of a Cirular Fin')
xlabel('meters (m)')
ylabel('Temperature (C)')


