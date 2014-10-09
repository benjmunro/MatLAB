%------- ENGR401 - Property Transported by Convetion and Diffusion -------%
%
%
% Ben Munro
% Developed : 28/04/2013
% Last Modified : 28/04/2013
%
%
% Uses HYBRID DIFFERENCING SCHEME
%
% V & M Example 5.1
% Discription : Solve for temperature along a rod wih 1D heat conduction,
% no internal source or heat loss and specified end temperatures.

clc
clear

%# Known Values
L = 1;
dx = 0.02;
rho = 1;
PhiA = 1;
PhiB = 0;
u = 1;
Gamma = 0.1;

%# Create Mesh
x = dx/2:dx:L;
N=length(x);

%# Calculated Constants
Dw = Gamma/dx;      De = Gamma/dx;
Fw = rho * u;       Fe = rho * u;
Area = dx * dx;

%# Create Matrix
A = zeros(N,N);
b = zeros(N,1);
ae = max([-Fe,(De-Fe/2),0]);
aw = max([Fw,Dw+Fw/2,0]);
for i = 1:N                         
    
    if x(i) == min(x)                       %# Boundary Point 1
        
        A(i,i) = ae + (2*De+Fe);
        A(i,i+1) = -ae;
        b(i,1) = (2*De+Fe)*PhiA;
        
    elseif x(i) == max(x)                   %# Boundary Point 2
        
        A(i,i) = aw + (2*Dw);
        A(i,i-1) = -aw;
        b(i,1) = (2*Dw)*PhiB;
        
    else                                    %# Internal Nodes
        A(i,i) = ae + aw + (Fe - Fw);
        A(i,i+1) = -ae;
        A(i,i-1) = -aw;
        
    end
    
end

%# Solve linear system of equations For Temp along bar
Phi = A\b;
disp('Phi along rod')
disp(Phi)

%# Plot Results
plot([0,x,L],[PhiA;Phi;PhiB],'o-b')
title('Property Transported by Convetion and Diffusion')
xlabel('meters (m)')
ylabel('Phi')


