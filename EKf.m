%EKF by reza izanloo
%y: measurment vector R:covariance matrix of mesurment noise  
%Q:covariance matrix of state noise   

xhat = [ 0 ; 0];
P = 0.1 * ones(2);

for k=2:N
xhat(1,k)=xhat(2,k-1);
xhat(2,k) = xhat(1,k-1)^2 + xhat(2,k-1) ;
xhat(:,k) = xhat(:,k-1) ;
A = [0 1; 2*xhat(1,k)  1];
P = A * P * A' + Q;

rk= y(1,k) - (xhat(1,k)+xhat(2,k)^3); 
C = [1 3*xhat(2,k)^2];
Gain = P*C'/(C*P*C'+R);
xhat(:,k) = xhat(:,k) + Gain * rk;
P = (eye(2)-Gain*C)*P;
end

%rmseKalman=sqrt(sum(eKalman));
hold on;
title('EKF Measurement:blue , Prediction:red , State:green')
t=1:1:100;
plot(t,xhat(1,t),'r');
hold on
plot(t,x(1,t),'g');
figure,
plot(t,xhat(2,t),'r');
hold on
plot(t,x(2,t),'g');

hold off;