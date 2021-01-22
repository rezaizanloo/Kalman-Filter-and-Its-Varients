% standard kalman filter by reza izanloo
function [xu,pu] = updateK(xp1,pp,y,R,C)
Gain=pp*C'/(C*pp*C'+R);
xu=xp1+Gain*(y-C*xp1);
pu=(eye(1)-Gain*C)*pp;
end