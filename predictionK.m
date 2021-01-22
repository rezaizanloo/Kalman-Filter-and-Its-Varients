% standard kalman filter by reza izanloo
function [xp,pp] = predictionK(xu2,pu,Q,A,k)
xp=A*xu2;
pp=A*pu*A'+Q;
end