function [ NodeMaxrix, N, xnodesb,xnodess,ynodes,ynodess] = Meshing( dx,dy ) 

%# Grid Sizing
ymax = 16;
xmax = 15;
%# Number of Nodes
%# Count
i = 0;
for y= dy:dy:ymax-dy/2
    if y <= 10-dy/2
      for x = dx:dx:xmax-dx/2
        i=i+1;
        NodeMaxrix(i,1:2) =[x,y]; 
      end      
    else  
        for x = 10+dx:dx:xmax-dx/2
        i=i+1;
        NodeMaxrix(i,1:2) =[x,y]; 
        end
    end 
end


ynodess=length(dy:dy:ymax-dy/2);
ynodes=length(dy:dy:10-dy/2);
xnodess=length(10+dx:dx:xmax-dx/2);
xnodesb=length(dx:dx:xmax-dx/2);

N=length(NodeMaxrix);
hold all
plot(NodeMaxrix(:,1),NodeMaxrix(:,2),'+k')
plot([0 15 15] , [0 0 16] ,[10 10 0 ],[16 10 10 ],'color','r','linewidth',2 ) 
axis([-5 20 -5 20])
end