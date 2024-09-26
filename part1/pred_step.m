function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time

% Assign Orientation(roll,pitch,yaw) of prev state
r = uPrev(4,1);
p = uPrev(5,1); 
y = uPrev(6,1); 
% Assign Linear velocity values of the previous state
vx = uPrev(7,1);
vy = uPrev(8,1);
vz = uPrev(9,1);
% Assign Gyroscope values from previous state
wx = angVel(1,1);
wy = angVel(2,1);
wz = angVel(3,1);
% Gyroscope bias bg
bgx = uPrev(10,1);
bgy = uPrev(11,1);
bgz = uPrev(12,1);
% Gyroscope noise
ngx = 0;
ngy = 0;
ngz = 0;
% Assign Accelerometer values from previous state
ax = acc(1,1);
ay = acc(2,1);
az = acc(3,1);
% Accelerometer bias ba
bax = uPrev(13,1);
bay = uPrev(14,1);
baz = uPrev(15,1);
% Accelerometer noise
nax = 0;
nay = 0;
naz = 0; 

%Combining data into single variable for better manipulation
lvel = [vx;vy;vz];
bg = [bgx;bgy;bgz];
ba = [bax;bay;baz];

% Calculate inverse of Euler Rate Parameterisation (ZYX) for f(x,u,n)
G_inv = [(cos(y)*sin(p))/(cos(p)*cos(y)^2 + cos(p)*sin(y)^2), (sin(p)*sin(y))/(cos(p)*cos(y)^2 + cos(p)*sin(y)^2), 1;
                               -sin(y)/(cos(y)^2 + sin(y)^2),                        cos(y)/(cos(y)^2 + sin(y)^2), 0;
                  cos(y)/(cos(p)*cos(y)^2 + cos(p)*sin(y)^2),          sin(y)/(cos(p)*cos(y)^2 + cos(p)*sin(y)^2), 0];
 
% Calculate Rotation matrix (ZYX) for f(x,u,n)
R=[cos(p)*cos(y),  cos(y)*sin(p)*sin(r) - cos(r)*sin(y),  sin(r)*sin(y) + cos(r)*cos(y)*sin(p)
   cos(p)*sin(y),  cos(r)*cos(y) + sin(p)*sin(r)*sin(y),  cos(r)*sin(p)*sin(y) - cos(y)*sin(r)
         -sin(p),                         cos(p)*sin(r),                         cos(p)*cos(r)];

% Concatenate all elements of f(x,u,n)
fxun = vertcat(lvel,G_inv * (angVel - bg),[0;0;-9.8] + R *(acc - ba),zeros(6,1));

%At=jacobian(fxun,[pos,ori,lvel,bg,ba]); was used to calculate the State
%transition matrix
At= [0, 0, 0,                                                                                                              0,                                                                                                              0,                                                                                                                                                  0, 1, 0, 0,  0,       0,              0,              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                                              0,                                                                                                                                                  0, 0, 1, 0,  0,       0,              0,              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                                              0,                                                                                                                                                  0, 0, 0, 1,  0,       0,              0,              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                -cos(p)*cos(y)*(bgz + ngz - wz),                                                                                           cos(y)*(bgy + ngy - wy) + sin(p)*sin(y)*(bgz + ngz - wz), 0, 0, 0,  0,  sin(y), -cos(y)*sin(p),              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                 sin(p)*sin(y)*(bgz + ngz - wz),                                                                                           sin(y)*(bgy + ngy - wy) - cos(p)*cos(y)*(bgz + ngz - wz), 0, 0, 0,  0, -cos(y), -cos(p)*sin(y),              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                        sin(p)*(bgz + ngz - wz),                                                                                                                                                  0, 0, 0, 0, -1,       0,        -cos(p),              0,                                      0,                                      0
0, 0, 0, - (sin(r)*sin(y) + cos(r)*cos(y)*sin(p))*(bay - ay + nay) - (cos(r)*sin(y) - cos(y)*sin(p)*sin(r))*(baz - az + naz), cos(y)*sin(p)*(bax - ax + nax) - cos(p)*cos(r)*cos(y)*(baz - az + naz) - cos(p)*cos(y)*sin(r)*(bay - ay + nay), (cos(r)*cos(y) + sin(p)*sin(r)*sin(y))*(bay - ay + nay) - (cos(y)*sin(r) - cos(r)*sin(p)*sin(y))*(baz - az + naz) + cos(p)*sin(y)*(bax - ax + nax), 0, 0, 0,  0,       0,              0, -cos(p)*cos(y),   cos(r)*sin(y) - cos(y)*sin(p)*sin(r), - sin(r)*sin(y) - cos(r)*cos(y)*sin(p)
0, 0, 0,   (cos(y)*sin(r) - cos(r)*sin(p)*sin(y))*(bay - ay + nay) + (cos(r)*cos(y) + sin(p)*sin(r)*sin(y))*(baz - az + naz), sin(p)*sin(y)*(bax - ax + nax) - cos(p)*cos(r)*sin(y)*(baz - az + naz) - cos(p)*sin(r)*sin(y)*(bay - ay + nay), (cos(r)*sin(y) - cos(y)*sin(p)*sin(r))*(bay - ay + nay) - (sin(r)*sin(y) + cos(r)*cos(y)*sin(p))*(baz - az + naz) - cos(p)*cos(y)*(bax - ax + nax), 0, 0, 0,  0,       0,              0, -cos(p)*sin(y), - cos(r)*cos(y) - sin(p)*sin(r)*sin(y),   cos(y)*sin(r) - cos(r)*sin(p)*sin(y)
0, 0, 0,                                                     cos(p)*sin(r)*(baz - az + naz) - cos(p)*cos(r)*(bay - ay + nay),                      cos(p)*(bax - ax + nax) + cos(r)*sin(p)*(baz - az + naz) + sin(p)*sin(r)*(bay - ay + nay),                                                                                                                                                  0, 0, 0, 0,  0,       0,              0,         sin(p),                         -cos(p)*sin(r),                         -cos(p)*cos(r)
0, 0, 0,                                                                                                                   0,                                                                                                              0,                                                                                                                                                  0, 0, 0, 0,  0,       0,              0,              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                                              0,                                                                                                                                                  0, 0, 0, 0,  0,       0,              0,              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                                              0,                                                                                                                                                  0, 0, 0, 0,  0,       0,              0,              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                                              0,                                                                                                                                                  0, 0, 0, 0,  0,       0,              0,              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                                              0,                                                                                                                                                  0, 0, 0, 0,  0,       0,              0,              0,                                      0,                                      0
0, 0, 0,                                                                                                                   0,                                                                                                              0,                                                                                                                                                  0, 0, 0, 0,  0,       0,              0,              0,                                      0,                                      0];
 
%Ut=jacobian(fxun,[ng,na,nba,nbg]); was used to calculate the Control Input matrix
Ut = [0,                                                    0,  0,              0,                                      0,                                      0, 0, 0, 0, 0, 0, 0
      0,                                                    0,  0,              0,                                      0,                                      0, 0, 0, 0, 0, 0, 0
      0,                                                    0,  0,              0,                                      0,                                      0, 0, 0, 0, 0, 0, 0
          (cos(p)*cos(y))/(sin(p)*cos(y)^2 + cos(p)*sin(y)^2),    (cos(p)*sin(y))/(sin(p)*cos(y)^2 + cos(p)*sin(y)^2), -1,                                  0,0,0, 0, 0, 0, 0, 0, 0
          (cos(p)*sin(y))/(sin(p)*cos(y)^2 + cos(p)*sin(y)^2),   -(cos(y)*sin(p))/(sin(p)*cos(y)^2 + cos(p)*sin(y)^2),  0,                                  0,0,0, 0, 0, 0, 0, 0, 0
                  -cos(y)/(sin(p)*cos(y)^2 + cos(p)*sin(y)^2),            -sin(y)/(sin(p)*cos(y)^2 + cos(p)*sin(y)^2),  0,                                  0,0,0, 0, 0, 0, 0, 0, 0
      0,                                                    0,  0, -cos(p)*cos(y),   cos(r)*sin(y) - cos(y)*sin(p)*sin(r), - sin(r)*sin(y) - cos(r)*cos(y)*sin(p), 0, 0, 0, 0, 0, 0
      0,                                                    0,  0, -cos(p)*sin(y), - cos(r)*cos(y) - sin(p)*sin(r)*sin(y),   cos(y)*sin(r) - cos(r)*sin(p)*sin(y), 0, 0, 0, 0, 0, 0
      0,                                                    0,  0,         sin(p),                         -cos(p)*sin(r),                         -cos(p)*cos(r), 0, 0, 0, 0, 0, 0
      0,                                                    0,  0,              0,                                      0,                                      0, 1, 0, 0, 0, 0, 0
      0,                                                    0,  0,              0,                                      0,                                      0, 0, 1, 0, 0, 0, 0
      0,                                                    0,  0,              0,                                      0,                                      0, 0, 0, 1, 0, 0, 0
      0,                                                    0,  0,              0,                                      0,                                      0, 0, 0, 0, 1, 0, 0
      0,                                                    0,  0,              0,                                      0,                                      0, 0, 0, 0, 0, 1, 0
      0,                                                    0,  0,              0,                                      0,                                      0, 0, 0, 0, 0, 0, 1];

% Initialize Ft, Vt, Q, Qd(covariance matrix of process noise) for the discretization step
Ft = eye(15) + dt*At;
Vt = Ut;
Q = 0.01;
Qd = dt*(eye(12) * Q);

% Compute the predicted state
uEst = uPrev + (dt*fxun);

% Calcualte predicted covariance
covarEst = (Ft * covarPrev * Ft') + (Vt * Qd * Vt');
end

