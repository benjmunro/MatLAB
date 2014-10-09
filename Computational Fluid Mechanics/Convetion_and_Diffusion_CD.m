%------- ENGR401 - Property Transported by Convetion and Diffusion -------%
%
%
% Ben Munro
% Developed : 28/04/2013
% Last Modified : 28/04/2013
%
% Uses CENTRAL DIFFERENCING SCHEME
%
% V & M Example 5.1
% Discription : Solve for temperature along a rod wih 1D heat conduction,
% no internal source or heat loss and specified end temperatures.

clc
clear

%# Known Values
L = 1;
dx = 1/20;
rho = 884;
PhiA = 100;
PhiB = 20;
Q=1.7e-12;
k =0.144;
c=1910;
Gamma = 0.1;

%# Create Mesh
x = dx/2:dx:L;
N=length(x);

%# Calculated Constants

Area =linspace((0.01/2)^2*pi,(0.005/2)^2*pi,N);
u = Q./Area;



%# Create Matrix
A = zeros(N,N);
b = zeros(N,1);

for i = 1:N
    
    if x(i) == min(x)
        
        
        Dw = k*Area(i)/dx;      De = k*Area(i+1)/dx;
        Fw = Q*rho*c;       Fe = Q*rho*c;
        aw = Dw + (Fw/2);   ae = De - (Fe/2);
        
        A(i,i) = ae + (2*De+Fe);
        A(i,i+1) = -ae;
        b(i,1) = (2*De+Fe)*PhiA;
        
    elseif x(i) == max(x)
        
        Dw = k*Area(i-1)/dx;      De = k*Area(N)/dx;
        Fw = Q*rho*c;       Fe = Q*rho*c;
        aw = Dw + (Fw/2);   ae = De - (Fe/2);
        
        A(i,i) = aw + (2*Dw-Fw);
        A(i,i-1) = -aw;
        b(i,1) = (2*Dw-Fw)*PhiB;
        
    else
        Dw = k*Area(i-1)/dx;      De = k*Area(i+1)/dx;
        Fw = Q*rho*c;       Fe = Q*rho*c;
        aw = Dw + (Fw/2);   ae = De - (Fe/2);
        A(i,i) = ae + aw + Fe - Fw;
        A(i,i+1) = -ae;
        A(i,i-1) = -aw;
        
    end
    
end

%# Solve linear system of equations For Phi along bar
Phi = A\b;
disp('Phi along rod')
disp(Phi)

%# Plot Results
plot([0,x,L],[PhiA;Phi;PhiB],'o--k')
title('Property Transported by Convetion and Diffusion')
xlabel('meters (m)')
ylabel('Phi')


