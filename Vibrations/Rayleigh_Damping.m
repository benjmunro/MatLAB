



z1 = [0.5  ;0.015  ];

ws=[16.7432 91.7515].^.5; %# if from lamda
w1=ws(1);
w2=ws(2);

w=.5*[ 1/w1 w1;1/w2 w2];

ab=w\z1;

%# get alpha and beta
alp = ab(1);
bet = ab(2);

%# calculate matracies
D=alp*M+bet*K;
Dstar=alp+bet*lambda;