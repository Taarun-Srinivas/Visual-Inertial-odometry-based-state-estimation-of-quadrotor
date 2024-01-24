function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state

    Ct = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];

R = diag([1e-3 0.01 0.1 0.01 0.1 0.01]);

Kt = covarEst*transpose(Ct)*inv(Ct*covarEst*transpose(Ct) + R);

uCurr = uEst + Kt*(z_t - Ct*uEst);
covar_curr = covarEst - Kt*Ct*covarEst;
    
end

