% Z : sequence of main signal (size_seq * Nstate)
% ZHAT : sequence of estimated signal (size_seq * Nstate)
% PP : sequence of error covariance matrix (size_seq * Nstate * Nstate)
% AA : sequence of  process equation matrix (size_seq * Nstate * Nstate) 

Nstate = 2;
Nmeas = 1;
z = z.';
zhat = zhat.';
size_seq = size(z,1);
% PP2 = zeros(size_seq,12,12);
% for i=1: size_seq
%     temp = P(i,:);
% PP2(i,:,:) = reshape(temp,12,12);
% end
% disp(['|| z0-z0hat|| <=', num2str(norm(ZHAT(1,:)-Z(1,:))),'']);
% disp(['Qfilter <= ', num2str(max(svd(Qfilter))),'']);
% disp(['Rfilter <= ',num2str(max(svd(R))),'']);

%% assumption 1 (earn  :  c_bar   p_down_bar    p_up_bar   q_bar  r_bar)

c_bar = max(svd(H)) %  (30)

s = [];
for i=1:2: size_seq
% s=[s eig(reshape(PP(i,:,:),Nstate,Nstate))];
s =[s svd(reshape(PP(i,:,:),Nstate,Nstate))];
end

p_down_bar = min(min(s)) % (31)
p_up_bar = max(max(s)) %  (31)


q_bar = min(svd(Q)) %  (32)
r_bar = min(svd(R)) %   (33)
%% assumption 2 (earn  :  upper_norm_z   epsilon_phi   epsilon_x   K_phi  K_x)
upper_norm_z = zeros(1,size_seq);
temp  = zeros(1,size_seq);
 for i= 1:size_seq
upper_norm_z(1,i) = norm(z(i,:));
  temp(1,i) = norm(z(i,:)-zhat(i,:));
 end
 disp('for all z  we have   ||z|| < up_norm_z ');
up_norm_z = max(upper_norm_z)
% disp('for all z  we have   ||z-zhat|| < epsilon_phi    and    ||z-zhat|| < epsilon_x ' );
epsilon_phi = max(temp)
epsilon_x = max(temp)
epsilon_hat = min(temp);

% K = sub_set_K()';  % dimension of K :  n*2
% K = Z;
set_K;
size_seq = size(K,2);

for i=1:Nstate
    temp2 = zeros(1,size_seq);
%     temp2(1) = norm(F1,2);
    for j=2:size_seq
        temp1 = CalcF_result(i,K(:,j));
         s = norm(temp1,2);
%         s  = max(svd(temp1));
%         s  = max(eig(temp1));
    temp2(j) = s ; % The spectral norm of a matrix A is the largest singular value of A
    end
    sup_hess(i) =max(temp2);
end

K_phi = max(sup_hess)
% K_phi = max(sup_hess(1:4))
K_x = 0 % because H is linear

k_nonl = (2*K_phi/p_down_bar)+(2*c_bar*K_x/r_bar);
epsilon_prim = min(epsilon_phi,epsilon_x);
epsilon = min(epsilon_prim,q_bar/(2*k_nonl*p_up_bar^2))   %50

k_noise = (12/p_down_bar) + ((p_up_bar^2 * c_bar^2 * 4)/(r_bar^2 * p_down_bar));
delta = (q_bar*p_down_bar*epsilon_hat^2)/(2*p_up_bar^3*k_noise)   %52









