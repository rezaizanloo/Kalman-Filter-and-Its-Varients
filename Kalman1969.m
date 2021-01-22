% standard kalman filter by reza izanloo
%y: measurment vector R:covariance matrix of mesurment noise  
%Q:covariance matrix of state noise   
%y(k)=C*x(k)+v(k)
%x(k+1)=A*x(k)+w(k);
n=length(y);
xu=ones(2,n); 
pu=ones(2,2);
xp=ones(2,n);
pp=ones(2,2);
xu(:,1)=[5,3];
pu=[2 1.5;2 2];
for k=2:1:n
[xp(:,k),pp]=predictionK(xu(:,k-1),pu,Q,A,k);
[xu(:,k),pu]=updateK(xp(:,k),pp,y(:,k),R,C);
end

%rmseKalman=sqrt(sum(eKalman));
hold on;
title('Kalman filter 3-dimension Measurement:blue , Prediction:red , State:green')
t=1:1:N;
plot3(t,y(1,t),y(2,t),'b');
hold on
plot3(t,xu(1,t),xu(2,t),'r');
hold on
plot3(t,x(1,t),x(2,t),'g');

hold off;