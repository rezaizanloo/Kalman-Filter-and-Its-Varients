% standard kalman filter by reza izanloo
%x(k+1)=A*x(k)+w(k)
%y(k)=C*x(k)+v(k)
clear all;
N=1000;
A=[1 0 
   0 1 ];
C=[1 0 
    0 1];
w = wgn(2,N,-50);
v = wgn(2,N,-50);

Q=cov(w');%covw
R=cov(v');%covv
x(1:2,1)=0;
y(1:2,1)=0;
for k=1:N
x(:,k+1)=A*x(:,k)+w(:,k);
y(:,k)=C*x(:,k)+v(:,k); 
end
t=1:1:N+1;
plot(t,x,'B')
hold on
plot(t,y,'R')