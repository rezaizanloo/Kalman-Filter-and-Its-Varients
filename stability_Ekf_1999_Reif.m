% article 2000 by k. Reif
% statbility of EKF
clc
clear all
close all
%%%%%%%%%%%%%%%%%%%5
%% define variables
Nstate = 2 ; Nmeas = 1 ; Nsample = 10; T =0.001;
PP = zeros(Nsample/T,Nstate,Nstate);
FF = zeros(Nsample/T,Nstate,Nstate);
zhat = zeros(2,Nsample/T);
z = zeros(2,Nsample/T);
y = zeros(1,Nsample/T);
P= eye(Nstate);
PP(1,:,:) = P;
z(:,1) = [ 0.8 ; 0.2 ];
Q = T*eye(Nstate);%Gt*Gt';
R = 1/T; %Dt*Dt';

% % example 1  (bounded)
% Gt = 0.1 * eye(Nstate);
% Dt = 0.1 ;
% zhat(:,1) = [0.5,0.5];

% Gt = 0 * eye(Nstate);
% Dt = sqrt(10^3) ;

Gt =sqrt(0.00001) * eye(Nstate);
Dt = sqrt(10) ;

zhat(:,1) =[0.5,0.5];
norm_z0 = norm(zhat(:,1)-z(:,1)); 

% % example 3  (divergent)
% Gt = 0.1 * eye(Nstate);
% Dt = 0.1 * eye(Nmeas);
% zhat(:,1) = [1.5;1.5];

% Q = Gt*Gt';
% R = Dt*Dt';

C = [1 0];

%% system simulation
for t = 2 : Nsample/T
    z(1,t) = z(1,t-1)+T*z(2,t-1) ;
    z(2,t) = z(2,t-1)+T*(- z(1,t-1) + (z(1,t-1)^2 + z(2,t-1)^2 - 1) * z(2,t-1));
    temp =  Gt * randn(Nstate,1);
    z(:,t) = z(:,t) +  temp;
    
    y(:,t) = C*z(:,t) +  Dt * randn(Nmeas,1);
end
PPP =P;

%% EKF
for t = 2 : Nsample/T
    %% prediction
    temp1 = zhat(1,t-1)+T*zhat(2,t-1) ;
    temp2 =zhat(2,t-1) + T * (- zhat(1,t-1) + (zhat(1,t-1)^2 + zhat(2,t-1)^2 - 1) * zhat(2,t-1));
    zhat(:,t) = [temp1;temp2];
    F = [1   T;(-1+2*zhat(1,t-1)*zhat(2,t-1))*T    1+(zhat(1,t-1)^2+3*zhat(2,t-1)^2-1)*T ];
    FF(t,:,:) = F;
    P = F * P * F' + Q;
    %% update
    S = C * P * C' + R;
    Gain = P  * C' * inv(S);
    zhat(:,t) = zhat(:,t) + Gain * (y(:,t) - C *zhat(:,t));
    P= (eye(Nstate) - Gain*C)*P;
    PPP =[PPP P] ;    
end
for i=1:2:Nstate*(Nsample/T)
    PP(i,:,:) = PPP(:,i:i+1);
end

figure,
hold on
axis([0 10000 -1 10])
t = 1:1 : Nsample/T;
c = 2; %  first element of state vector
plot(t,z(c,t) ,'blue');
plot(t,zhat(c,t),'green');
legend('exat solution','EKF');
xlabel('seconds');
ylabel('position');

%% stability analysis
size_seq = size(z,2);

%% assumption 1 (earn  :  c_bar   p_down_bar    p_up_bar   q_bar  r_bar)
temp = zeros(1,size_seq);
for i=1: size_seq
    s=svd(reshape(FF(i,:,:),Nstate,Nstate));
     temp(1,i) = max(s);
end
a_bar = max(temp); % (28)
c_bar = max(svd(C)); %  (29)

s = [];
temp = zeros(1,size_seq);
for i=1:2: size_seq
%     temp(1,i) = norm(reshape(PP(i,:,:),Nstate,Nstate));
s= [s svd(reshape(PP(i,:,:),Nstate,Nstate))];
% s =[s eig(reshape(PP(i,:,:),Nstate,Nstate))];
end

p_down_bar = min(min(s)) % (30)
p_up_bar = max(max(s)) %  (30)

% q_bar = min(svd(Gt*Gt')) %  (31)
% r_bar = min(svd(Dt*Dt')) %   (32)

q_bar = min(svd(Q)); %  (31)
r_bar = min(svd(R)); %   (32)

%% assumption 2 (earn  :  upper_norm_z   epsilon_phi   epsilon_x   K_phi  K_x)
upper_norm_z = zeros(1,size_seq);
temp  = zeros(1,size_seq);
 for i= 1:size_seq
upper_norm_z(1,i)=norm(z(:,i),2);
  temp(1,i) = norm(z(:,i)-zhat(:,i));
 end
%  disp('for all z  we have   ||z|| < up_norm_z ');
up_norm_z = max(upper_norm_z);
% disp('for all z  we have   ||z-zhat|| < epsilon_phi    and    ||z-zhat|| < epsilon_x ' );
epsilon_phi = max(temp);
epsilon_x = max(temp);
epsilon_hat = min(temp);

K = [];
for i=-10 : 0.01:10
for j= -10 : 0.01 :10
    if norm([i j])<=10
        K = [K [i;j]];
    end
end
end
size_seq = size(K,2)

for i=1:Nstate
    temp2 = zeros(1,size_seq);
    for j=1:size_seq
        if i==1
        temp1 = zeros(2);
        else
            temp1 =  [2*T*K(2,j)  2*T*K(1,j) ;2*T*K(1,j)  6*T*K(2,j) ];
        end
%         s  = max(eig(temp1));
 s = norm(temp1,2);
%  s = max(svd(temp1));
    temp2(j) = s ; % The spectral norm of a matrix A is the largest singular value of A
    end
    sup_hess(i) =max(temp2);
end

% K_phi = max(sup_hess)
K_phi = max(sup_hess)
K_x = 0 % because H is linear

temp = (p_up_bar)*(a_bar+(a_bar*p_up_bar*c_bar^2)/r_bar)^2;
alpha = 1 - 1/ ( 1 + (q_bar/temp) ) ;
k_prim = K_phi + a_bar*p_up_bar*c_bar*(1/r_bar)*K_x;
epsilon_prim = min(epsilon_phi,epsilon_x);
k_nonl =k_prim*(1/p_down_bar)*(2*(a_bar+a_bar*p_up_bar*c_bar^2/r_bar)+(k_prim*epsilon_prim));
epsilon = min(epsilon_prim,alpha/(2*k_nonl*p_up_bar))   

k_noise = (Nstate/p_down_bar) + ((a_bar^2*p_up_bar^2 * c_bar^2*Nmeas )/(r_bar^2 * p_down_bar));
delta = (alpha*epsilon_hat^2)/(2*p_up_bar*k_noise) 

norm(zhat(:,1)-z(:,1))
norm(zhat(:,1)-z(:,1)) <= epsilon
max(svd(Gt*Gt'))
max(svd(Gt*Gt')) <= delta
max(svd(Dt*Dt')) 
max(svd(Dt*Dt')) <= delta
