%# Caughey Damping





[X_hat, lambda] = eig(K, M);

%# Zi values
zi = [0.1 0.05 0.03];

%# Frequencys
w = [30.9574 196.7946  382.2480  ].^0.5;

Dstardiag= 2*zi.*w;

Dstar=diag(Dstardiag);

M;
X_hat;

D=M*X_hat*Dstar*X_hat'*M;