% unscented kalman filter by reza izanloo
%% Example 1
% num_vec = 3; % dimension of state vector
% num_meas = 2; % number of measurment
% Q=0.2*eye(num_vec); % covariance of process
% R=0.2*eye(num_meas);        % covariance of measurement 
% x=[0;0;1]; 
% xhat=[0;0;1]; % initial state
% P = eye(num_vec);
% S = chol(P);        % initial state covraiance
% f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
% h=@(x)[x(1); x(3)];                               % measurement equation
%%
teta = pi/3;
T = 3 ;
num_vec = 4; % dimension of state vector
num_meas = 2; % number of measurment
Q = diag([0.1,0.1,0.1,0.1]); % cov of Noise1 process
R = diag([0.1,0.1]);% cov of Noise1 measurment
A =[1 0 T 0  % process Equation
    0 1 0 T
    0 0 1 0
    0 0 0 1];
B =[1 0 0 0  % meaurment equation
    0 1 0 0];
x = [1;1;0;0]; % initial of state vector X
xhat = [1;1;0;0]; % initial of state vector X
P = diag([4,4,3,3]);% initial of cov matrix P
S = chol(P);
f = @(x)[x(1)+T*x(3);x(2)+T*x(4);x(3);x(4)];  % nonlinear state equations for UKF
h = @(x)[x(1);x(2)];                               % measurement equation for UKF
%%
N=50;                                     % total dynamic steps
xV = zeros(num_vec,N);          %estmate        % allocate memory
sV = zeros(num_vec,N);          %actual
zV = zeros(num_meas,N);
for k=1:N
      x = f(x) + sqrt(Q)*randn(num_vec,1);                % update process 
  z =h(x) + sqrt(R)*randn(num_meas,1);                     % measurments
  sV(:,k)= x;                             % save actual state
  zV(:,k)  = z;                             % save measurment
%   [xhat, P] = ukf(f,xhat,P,h,z,Q,R);            % ekf 
%    [x, P] = UKBF(f,x,P,h,z,Q,R); 
%       [x,S] = SR_UKF(f,x,S,h,z,Q,R); 

  xV(:,k) = x;                            % save estimate
end

for k=1:3                                 % plot results
  subplot(3,1,k)
  plot(1:N, sV(k,:), 'b', 1:N, xV(k,:), 'g')
end



