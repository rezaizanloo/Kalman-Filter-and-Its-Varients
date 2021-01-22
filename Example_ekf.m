% Extended kalman filter by reza izanloo
%x(k+1)=f(x(k))+w(k)      f  : x(1,k+1)=x(2,k)  and x(2,k+1) = x(1,k)^2 + x(2,k)
%y(k)=h(x(k))+v(k)   h : y(k) = x(1,k)+x(2,k)^3
clear all;
N=1000;
w = wgn(2,1000,-50);
v = wgn(1,1000,-20);
x = zeros(2,1000);
y = zeros(1,1000);

Q=cov(w');%covw
R=cov(v');%covv
x(1:2,1)=0;

for k=1:N
x(1,k+1)=x(2,k);
x(2,k+1) = x(1,k)^2 + x(2,k) ;
x(:,k+1) = x(:,k+1) + w(:,k);

y(1,k)=x(1,k)+x(2,k)^3+v(:,k); 
end
y(1,N)=x(1,N)+x(2,N)^3+v(:,N); 