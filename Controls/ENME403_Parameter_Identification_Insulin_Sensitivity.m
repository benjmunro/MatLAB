%--------- ENME403 Parameter Identification - Insulin Sensitivity---------%
%
% Written by Ben Munro
% Developed - 17/04/2013
% Last updated - 17/04/2013

%# Reset
clc
clear

%# A-priori Prameters
nT = 0.2;       %# Clearance rate (1/min)
nI = 0.05;      %# Plasma to interstitium transport rate (1/min)
nC = 0.05;      %# Cell metabolism of insulin (1/min)
xL = 0.7;       %# First pass extraction ()
pG = 0.004; 
VP = 4.3;  
time = 40;      %# Time (min)
t = 0:.1:time; %# Time vector (min)

%# Initial Condtitions
I0 = (50*(1-xL)/((VP*(nT+(nI/2)))));
Q0 = I0/2;

%# Initilize vectors
I=I0*ones(size(t));
Q=Q0*ones(size(t));
Ux=zeros(size(t));
Un=zeros(size(t));
Px=zeros(size(t));


for i = 1:length(t)
   
    if t(i)>=15 && t(i)<=16
        Ux(i) = 1000;
    end
        
    if t(i) < 6
        Un(i)=50;
    elseif t(i)>=6 && t(i)<=7
        Un(i)=500;
    else
        Un(i)=70;
    end
    
     if t(i)>=5 && t(i)<=6
        Px(i) = 10;
    end
    
end

%---------------- Forward Simulation of Insulin Model---------------------%

for j = 1:50;
    
    I = I0 + cumtrapz(t,-1*nT*I + nI*(Q-I) + ((Ux + (1-xL)*Un)/VP));
    Q = Q0 + cumtrapz(t,nI*(I-Q) - nC*Q);
    
end
hold all
figure(1)
plot(t,Q,'Color','yellow','LineWidth',2)
plot(t,I,'Color','blue','LineWidth',2)
legend('Interstitial Insulin','Plasma Insulin')
title('Forward Simulation of Insulin')
xlabel('Time (minutes) ')
ylabel('Insulin Concentration (mU.L-^1) ')


%--------------------------- Glucose Data --------------------------------%

%# Original Glucose data
G_T1 = [4.2 4.2 8.36 7.80 6.45 5.09 4.23 3.73 3.43];
G_T2 = [5.3 5.3 8.86 8.70 8.33 7.89 7.54 7.27 7.08];

%# Add noise to data
%G_T1N = G_T1+G_T1.*(0.05*randn(1,length(G_T1)));
%G_T2N = G_T2+G_T2.*(0.05*randn(1,length(G_T2)));

%# Ensuring a steady basal pre-bolus period
%G_T1N(2) = G_T1N(1);
%G_T2N(2) = G_T2N(1);

G_T1N = [4.3978 4.3978 8.2041 8.1180 6.7076 5.1206 4.3508 3.8070 3.2607];
G_T2N = [5.5013 5.5013 8.5925 8.7770 8.2019 7.8380 7.7644 7.6505 7.0099];

%# Time of samples
T_samples = 0:5:40;

%------------- Objective Surface Contours of Vg and Si -------------------%

  N=61;
  SI = linspace(0,2e-3,N);
  VG = linspace(10,20,N);
  
  
  for i = 1:length(T_samples)
      index(i)= find(t==T_samples(i));
  end
  
for simNo = 1:2

    %# Set the data to be the desired set
    if simNo == 1
        Grun = G_T1;
    else
        Grun = G_T2;
    end
    
    G0=Grun(1);
    Gsim=G0*ones(size(t));
    Gsim_e = zeros(size(T_samples));
    error = zeros (N,N);
    
for i = 1:N
    
   for j = 1:N
       
       for k = 1:50
       Gsim = G0 + cumtrapz(t,-pG*(Gsim-G0)-SI(i)*(Gsim.*Q-G0*Q0)+(1000*Px/(180*VG(j)))); 
       end
       
       Gsim_e = Gsim(index);
       error(i,j) = norm(Grun-Gsim_e);     
   end
end

figure(2)
subplot(1,2,simNo)
contour(VG,SI,error,100);
legend('Contours')
hold all
title(sprintf('Objective Function Contour for Data Set %d',simNo))

end

%----------------- IIM Method to Identify SI and VI ----------------------%

for simNo = 1:2
   
    if simNo == 1
        Grun = G_T1N;
   else
        Grun = G_T2N;
    end 
    
  
    
    %# Initial Values 
    Si=0.001;   %# What am I
    VG=15;      %# What am I
    G0=Grun(1);
    
    
    %# Initialize Matrix
    Gsim = G0*ones(size(t));
    Mat_1=zeros(length(T_samples)-1,2);
    Mat_2=zeros(length(T_samples)-1,1);
    SI_1=zeros(1,30);
    VG_1=zeros(1,30);
    
    
    %Loop for IIM
    for i = 1:30
        
        for j = 1:50
            Gsim = G0 + cumtrapz(t,-1*pG*(Gsim-G0) - Si*(Gsim.*Q - G0*Q0) + (1000*Px/(180*VG)));
        end
        
        %Loop to fill out IIM matrices for system of equations
        
        for k = 1:length(Grun)-1
            INDEX = find(t==T_samples(k+1));
            Mat_1(k,1)=trapz(t(1:INDEX),-1*(Gsim(1:INDEX).*Q(1:INDEX)-G0*Q0));
            Mat_1(k,2)=trapz(t(1:INDEX),(1000*Px(1:INDEX)/180));
            Mat_2(k)=Grun(k+1) - G0 + pG*trapz(t(1:INDEX),Gsim(1:INDEX)-G0);
        end
        
        %Solve for SI and VG
        Mat_3 = Mat_1\Mat_2;
        Si=Mat_3(1);
        VG=(1/Mat_3(2));
        SI_1(i)=Si;
        VG_1(i)=VG;
    end
    
    SIConv(simNo)=Si;
    VGConv(simNo)=VG;
    
    figure(2)
    figure(3)
    
    hold all
    subplot(1,2,1)
    plot(1:30,SI_1);
    xlabel('Iteration');
    ylabel('SI');
    
    title('Convergence of Residuals for SI')
    if simNo==2
        legend('Data Set 1','Data Set 2')
        axis([0 30 0 1.8e-3])
    end
    
    hold all
    subplot(1,2,2)
    plot(1:30,VG_1);
    xlabel('Iteration');
    ylabel('VG');
    title('Convergence of Residuals for VG')
    if simNo==2
        legend('Data Set 1','Data Set 2')
         axis([0 30 10 20])
    end
end

%------------------- Simulate Glucose Profile SI and VI ------------------%

for simNo=1:2
    
    if simNo==1
        Si=SIConv(1);
        VG=VGConv(1);
        Grun= G_T1N;;
    else
        Si=SIConv(2);
        VG=VGConv(2);
        Grun= G_T2N;;
    end
    
    G0=Grun(1);
    Gsim = G0*ones(size(t));
    for i = 1:30
        Gsim = G0 + cumtrapz(t,-1*pG*(Gsim-G0) - Si*(Gsim.*Q - G0*Q0) + (1000*Px/(180*VG)));
    end
    
    figure(4)
    subplot(1,2,simNo)
    plot(t,Gsim,T_samples,Grun,'k+')
    title(sprintf('Simulated Glucose Profile for Data Set %d with Optimum Parameter Values',simNo))
    legend('Model Glucose Values','Actual Glucose Values' )
    xlabel('Time (Minutes) ')
    ylabel('Glucose Concentration (mmol.L^-^1) ')
end

%----------- Levenberg Marquardt Method to Identify SI and VI ------------%

for simNo = 1:2
    
    if simNo == 1
        Grun = G_T1N;
    else
        Grun = G_T2N;
    end
    

    OPTIONS = optimset('Algorithm','levenberg-marquardt','TolX',1e-8);
    X0 = [0.001,15];
    
    %x=lsqnonlin(@fit_simp,X0,[],[],OPTIONS,T,G,t,Q,Q0,Px);
    x=lsqnonlin(@My_lm_fun,X0,[],[],OPTIONS,index,Grun,t,Q,Q0,Px,pG);
    
    SI_LM(simNo) = x(1);
    VG_LM(simNo) = x(2);
    
    figure(2)
    subplot(1,2,simNo)
    plot(VG_LM(simNo),SI_LM(simNo),'k+','markersize',14,'linewidth',3)
    subplot(1,2,simNo)
    plot(VGConv(simNo),SIConv(simNo),'rx','markersize',15,'linewidth',2)
    legend('Contours','Levenberg-Marquardt Parameter''IIM Optimal Value')
end

%------------------------ Montecarlo Investigation -----------------------%

for simNo = 1:1
   
%# Initial Values
Si=0.001;   %# What am I
VG=15;      %# What am I

%# Initialize Matrix

Mat_1=zeros(length(T_samples)-1,2);
Mat_2=zeros(length(T_samples)-1,1);
Montecarlo_MAT = zeros(100,2);

for y = 1:300
    
    if simNo == 1
    %# Add noise to data
    G_T1N = G_T1+G_T1.*(0.05*randn(1,length(G_T1)));
    
    %# Ensuring a steady basal pre-bolus period
    G_T1N(2) = G_T1N(1);
    Grun = G_T1N;
    G0=Grun(1);
    Gsim = G0*ones(size(t));
    else
         %# Add noise to data
    G_T2N = G_T2+G_T2.*(0.05*randn(1,length(G_T1)));
    
    %# Ensuring a steady basal pre-bolus period
    G_T2N(2) = G_T2N(1);
    Grun = G_T2N;
    G0=Grun(1);
    Gsim = G0*ones(size(t));
    end
    
    %# Now run IIM method as previous
    
    %Loop for IIM
    for i = 1:30
        
        for j = 1:50
            Gsim = G0 + cumtrapz(t,-1*pG*(Gsim-G0) - Si*(Gsim.*Q - G0*Q0) + (1000*Px/(180*VG)));
        end
        
        %Loop to fill out IIM matrices for system of equations
        
        for k = 1:length(Grun)-1
            
            %# Index T_samples in t
            INDEX = find(t==T_samples(k+1));
            
            
            Mat_1(k,1)=trapz(t(1:INDEX),-1*(Gsim(1:INDEX).*Q(1:INDEX)-G0*Q0));
            Mat_1(k,2)=trapz(t(1:INDEX),(1000*Px(1:INDEX)/180));
            Mat_2(k)=Grun(k+1) - G0 + pG*trapz(t(1:INDEX),Gsim(1:INDEX)-G0);
        end
        
        %# Solve for SI and VG
        Mat_3 = Mat_1\Mat_2;
        Si=Mat_3(1);
        VG=(1/Mat_3(2));
        
    end
    
    %# Store Si an Vi Values
    Montecarlo_MAT(y,:) = [Si , VG];
    
    
end
figure(2)
subplot(1,2,simNo)
plot(Montecarlo_MAT(:,2),Montecarlo_MAT(:,1),'kx');

end

%# Statistical Analysis

mean_SI=mean(Montecarlo_MAT(:,1));
mean_VG=mean(Montecarlo_MAT(:,2));
[COEFF, SCORE, LATENT, TSQUARED] = princomp(Montecarlo_MAT);

Vari=std(Montecarlo_MAT)./mean(Montecarlo_MAT);
 

