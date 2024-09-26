function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively

%Initializing the C matrix for linearization and is not computed since the
%measurement model z_t is available
C = horzcat(zeros(3,6), eye(3), zeros(3,6)); %Dimension 3x15
% 6x6 matrix R for introducing Sensor noise for linear update
R = eye(3) * 0.00009;
% Kalman Gain Calculation
KalG = (covarEst * C') / ((C * covarEst * C') + R);
% Current Covariance matrix
covar_curr = covarEst - (KalG * C * covarEst);
% Current State
uCurr = uEst + (KalG * (z_t - (C*uEst)));

end