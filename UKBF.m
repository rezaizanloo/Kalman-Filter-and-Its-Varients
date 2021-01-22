% unscented kalman filter by reza izanloo
function [xhat,P]=UKBF(fstate,xhat,P,hmeas,y,Q,R)

L=numel(xhat);                                 %numer of states
m=numel(y);                                 %numer of measurements
alpha = 1e-1;                                 %default, tunable
ki=0;
beta=2;                                     %default, tunable
lambda=alpha^2*(L+ki)-L;                    %scaling factor
a = L+lambda;                                 %scaling factor
Wm=[lambda/a 0.5/a+zeros(1,2*L)];           %weights for means
Wc=Wm;
Wc(1)=(lambda)/(a+(1-alpha^2+beta));               %weights for covariance
% Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
w_m = Wm.';
W = (eye(2*L+1) - repmat(w_m,1,2*L+1)) * diag(Wc) * (eye(2*L+1) - repmat(w_m,1,2*L+1))';                  %scaling factor
% %% formula 26
% c = alpha^2*(L+ki);
% cholp = chol(P);
% X = repmat(xhat,1,2*L+1) + sqrt(c) * [zeros(L,1) cholp -cholp];
% X_hat = ut(fstate,X,L,Q);
% xhat = X_hat * w_m;
% P = X_hat * W * X_hat' + Q ;
% %% formula 27
% cholp = chol(P);
% X = X_hat + sqrt(c) * [zeros(L,1) cholp -cholp];
% Y=ut(hmeas,X,m,R);       %unscented transformation of measurments
% u_k = Y * w_m ;
% S = Y * W * Y' + R ;
% C = X_hat * W * Y' ;
% %% formula 28
% K = C * inv(S);
% xhat = xhat + K * (y - u_k);
% P = P - K * S * K';
%% formula 29 continious form
LL = eye(L);
V = eye(m);
Q_c = Q;
R_c = R;
c = alpha^2*(L+ki);
cholp = chol(P);
X = repmat(xhat,1,2*L+1) + sqrt(c) * [zeros(L,1) cholp -cholp];
X_hat = ut(fstate,X,L,Q);
Y=ut(hmeas,X,m,R);       %unscented transformation of measurments
u_k = Y * w_m ;
K = X * W * Y' * inv(V * R_c * V);
xhat = xhat + K * (y - u_k);
P = X * W * X_hat' + X_hat * W * X' + LL * Q_c * LL' - K * V * R_c * V' * K';

function [Y]=ut(f,X,n,R)
%Unscented Transformation
L=size(X,2);
Y=zeros(n,L);
for k=1:L                   
    Y(:,k)=f(X(:,k));              
end
  