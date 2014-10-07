%--------FVM Simulation of a Potential Flow around a 90° Corner ----------%
%
%
% Written by Ben Munro
% LAST MODIFIED:24/04/2013

clc
clear

%# Mesh Size
dx=1;
dy=1;

%# Geometry
ymax = 16;
xmax = 15;
ymid = 10;
xmid = 10;

Ui=10;
Uo=20;

%# Show Mesh
MeshMatrix=@Meshing;
[ NodeMaxrix, N, xnodesb,xnodess,ynodes,ynodess]= Meshing( dx,dy );
x = NodeMaxrix(:,1);
y = NodeMaxrix(:,2);

A=zeros(N,N);
b= zeros(1,N);
aw = dx*dy/dx; ae = aw; as = aw; an = aw;

for i = 1:N
    %# Are you one of my corner nodes???
    
    if y(i) == dy && x(i) == dx                               %# Corner One
        A(i,i) =aw+ae+as+an;
        A(i,i+xnodesb) = -an;
        A(i,i+1) = -ae;
        b(i) = (dy*Ui)*aw;
        
    elseif y(i) == dy && x(i) == xmax-dx                      %# Corner Two
        A(i,i) =aw+ae+as+an;
        A(i,i-1) = -aw;
        A(i,i+xnodesb) = -an;
        
    elseif y(i) == ymax-dy && x(i) == xmax-dx               %# Corner Three
        A(i,i) = aw+ae+as+an;
        A(i,i-1) = -aw;
        A(i,i-xnodess) = -as;
        b(i) = (dx*Uo)*an;
        
    elseif y(i) == ymax-dy && x(i) == xmid+dx                %# Corner Four
        A(i,i) = aw+ae+as+an;
        A(i,i+1) = -ae;
        A(i,i-xnodess) = -as;
        b(i) = 100*aw + ((5-dx) *Uo)*an;
        
        elseif y(i) == ymid-dy && x(i) == xmax-dx            %#Corner six
        
            A(i,i) =aw+ae+as+an;
            A(i,i-1) = -aw;
            A(i,i-xnodesb) = -as;
            A(i,i+xnodess) = -an;
       
        
        
    elseif y(i) == ymid-dy && x(i) == dx                       %#Corner six
        
        A(i,i) = aw+ae+as+an;
        A(i,i-xnodesb) = -as;
        A(i,i+1) = -ae;
        b(i)=(y(i)*10)*aw +100*aw;
        
        
        
    elseif  y(i) == dy                %# Now along the bottom boundary
        A(i,i) =aw+ae+as+an;
        A(i,i-1) = -aw;
        A(i,i+1) = -ae;
        A(i,i+xnodesb) = -an;
        
        
    elseif x(i) == xmax-dx                  %# Now along the Right boundary
        
        if y(i) < ymid
            A(i,i) =aw+ae+as+an;
            A(i,i-1) = -aw;
            A(i,i-xnodesb) = -as;
            A(i,i+xnodesb) = -an;
            
        else
            A(i,i) =aw+ae+as+an;
            A(i,i-1) = -aw;
            A(i,i-xnodess) = -as;
            A(i,i+xnodess) = -an;
        end
        
        
    elseif y(i) == 16-dy                                    %# Top Boundary
        
        A(i,i) = aw+ae+as+an;
        A(i,i+1) = -ae;
        A(i,i-xnodess) = -as;
        A(i,i-1) = -aw;
        b(i)=((15-x(i))*20)*an;
        
        
         
    elseif x(i) == xmid+dx && y(i) >= ymid             %# Left Mid Boundary
        
        A(i,i) = aw+ae+as+an;
        A(i,i+xnodess) = -an;
        A(i,i-xnodess) = -as;
        A(i,i+1) = -ae;
        b(i)=(100)*aw;
        
    elseif y(i) == ymid-dy && x(i) <= xmid              %# Top Mid Boundary
        
        A(i,i) = aw+ae+as+an;
        A(i,i-xnodesb) = -as;
        A(i,i+1) = -ae;
        A(i,i-1) = -aw;
        b(i)=(100)*aw;
        
    elseif x(i) == dx                                     %# left Boundary
        
        A(i,i) = aw+ae+as+an;
        A(i,i+xnodesb) = -an;
        A(i,i-xnodesb) = -as;
        A(i,i+1) = -ae;
        b(i)=(y(i)*10)*aw;
    
    else                                                   %#Internal Nodes
       
        if y(i) >= ymid
            
        A(i,i) = aw+ae+as+an;
        A(i,i-1) = -aw;
        A(i,i+1) = -ae;
        A(i,i+xnodess) = -an;
        A(i,i-xnodess) = -as;
        
        elseif y(i) == ymid-dy
            
            A(i,i) = aw+ae+as+an;
            A(i,i-1) = -aw;
            A(i,i+1) = -ae;
            A(i,i+xnodess) = -an;
            A(i,i-xnodesb) = -as;
        
        else
            
            A(i,i) = aw+ae+as+an;
            A(i,i-1) = -aw;
            A(i,i+1) = -ae;
            A(i,i+xnodesb) = -an;
            A(i,i-xnodesb) = -as;
            
        end
  
    end
         
end
A=sparse(A);
Phi=A\b';

%# Calculate Velocity at Node 5 using Backwards diffrence
Vel_Nodes_5 = find(y==5);
Vel_Nodes_5_dy = find(y==5-dy);
Phi_Vel_5=[50;Phi(Vel_Nodes_5);0];
Phi_Vel__dy=[(5+dy)*10;Phi(Vel_Nodes_5_dy);0];
Velocity = abs(Phi_Vel__dy-Phi_Vel_5)/(dy) ;
figure(2)
plot([0;x(Vel_Nodes_5);15],Velocity)

% Ploting Stream Lines
k=0;
for i=2:ynodes+1
    
    for j=2:xnodesb+1
        k=k+1;
        xx(i,j)=x(k);
        yy(i,j)=y(k);
        Ans(i,j)=Phi(k);
        
    end
    
end

for i=ynodes+2:ynodess+1
    
    for j=xnodesb-xnodess+2:xnodesb+1
        k=k+1;
        xx(i,j)=x(k);
        yy(i,j)=y(k);
        Ans(i,j)=Phi(k);
        
    end
end

Ans(Ans==0)=NaN;
yy(xx==0)=NaN;
xx(yy==0)=NaN;
 
contour(xx,yy,Ans,30,'LineWidth',1)




    
    
    
    
    
    
    
