% unscented kalman filter by reza izanloo
function [xhat,S]=SR_UKF(fstate,xhat,S,hmeas,y,Q,R)

L=numel(xhat);                                 %numer of states
m=numel(y);                                 %numer of measurements
alpha = 0.9;                                 %default, tunable
ki=2;
beta=2;                                     %default, tunable
lambda=alpha^2*(L+ki)-L;                    %scaling factor
a = L+lambda; 
gama = sqrt(a);
Wm=[lambda/a 0.5/a+zeros(1,2*L)];           %weights for means
Wc=Wm;
Wc(1)=(lambda)/(a+(1-alpha^2+beta));               %weights for covariance               %weights for covariance
% Wm(1)=Wc(1)+(1-alpha^2+beta);
Wm(1)=lambda/(L+lambda);
%% calculate sigma point and time update
X=sigmas(xhat,S,gama); 
[x1,X_star,P1,X2]=ut(fstate,X,Wm,Wc,L,Q);          %unscented transformation of process
[S,O] = qr([Wm(2)*(X2(:,2:end)) sqrt(Q)]); 
S = cholupdate(S,sqrt(Wm(1))*(X_star(:,1)-x1),'+'); 
X=sigmas(x1,S,gama);
[y1,Y1,P2,Y2]=ut(hmeas,X,Wm,Wc,m,R);       %unscented transformation of measurments
%% Measurement update equation
[S_y,O] = qr([(Wm(2)*(Y2(:,2:end))) sqrt(R)]); 
S_y = cholupdate(S_y,sqrt(Wm(1))*(Y1(:,1)-y1),'+'); 
P12=X2*diag(Wc)*Y2';                        %transformed cross-covariance
K = (P12/S_y')/S_y;
xhat=x1+K*(y-y1);
U = K * P12';
% U = K * S_y;
S = chol(S'*S+U); 
% S = chol(S'*S + Wm(1)*Y2(:,1)*Y2(:,1)');

function [y,Y,P,Y1]=ut(f,X,Wm,Wc,n,R)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations
L=size(X,2);
y=zeros(n,1);
Y=zeros(n,L);
for k=1:L                   
    Y(:,k)=f(X(:,k));       
    y=y+Wm(k)*Y(:,k);       
end
Y1=Y-y(:,ones(1,L));
P=Y1*diag(Wc)*Y1'+R;          

function X=sigmas(x,S,gama)
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points
Y = x(:,ones(1,numel(x)));
X = [x Y+(gama*S) Y-(gama*S)]; 
 