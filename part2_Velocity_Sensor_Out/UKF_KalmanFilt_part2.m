clear all; % Clear variables
addpath('../data')
datasetNum = 1; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime, proj2Data] = init(datasetNum);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.01*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0;
vel = proj2Data.linearVel;
angVel2 = proj2Data.angVel;
Zn = vertcat(vel',angVel2');

%% Calculate Kalmann Filter
for i = 1:length(sampledTime)

    if (i == 1)
        dt = sampledTime(1);
        [covarEst, uEst] = pred_step(uPrev, covarPrev, sampledData(1).omg, sampledData(1).acc, dt);
        [uCurr, covar_curr] = upd_step(Zn(:,1), covarEst, uEst);
        
    else
        dt = sampledTime(i) - sampledTime(i-1);
        [covarEst, uEst] = pred_step(uCurr, covar_curr, sampledData(i).omg, sampledData(i).acc, dt);
        [uCurr, covar_curr] = upd_step(Zn(:,i), covarEst, uEst);
        
    end

    
    uPrev = uCurr;
    covarPrev = covar_curr;
    savedStates(:,i) = uCurr;
end

plotData(savedStates, sampledTime, sampledVicon, 2, datasetNum);