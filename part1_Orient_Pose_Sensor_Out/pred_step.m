function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 

    alpha = 0.001;
    k = 1;
    beta = 2;
    Q = diag([1e-3 1e-3 1e-3 1e-3 1e-2 1e-3 1e3 1e3 1e3 1e-2 1e-2 1e-2]);

    uAug = [uPrev;zeros(12,1)];
    n_prime = length(uAug);
    Xt = zeros(15,1);
    final_mean = zeros(15,1);
    final_covar = zeros(15,15);

 
    L = alpha^2 *(n_prime + k) - n_prime;
    P_Aug = [covarPrev zeros(15,12);
             zeros(12,15)        Q];

    covarAug = chol(P_Aug,"lower");
    
    sigma_point(:,1) = uAug;

    for i = 1:n_prime
        xAug_odd = uAug + sqrt(n_prime + L)*covarAug(:,i);
        xAug_even= uAug - sqrt(n_prime + L)*covarAug(:,i);
        xAug = [xAug_odd xAug_even];
        sigma_point = horzcat(sigma_point,xAug);    
    end

    for i = 1:length(sigma_point)
        pos = sigma_point((1:3),i);
        q = sigma_point((4:6),i);
        lv = sigma_point((7:9),i);
        gb = sigma_point((10:12),i);
        ab = sigma_point((13:15),i);

        gn = sigma_point((16:18),i);
        an = sigma_point((19:21),i);
        gbn = sigma_point((22:24),i);
        abn = sigma_point((25:27),i);

        %System dyanmics
        pd = sigma_point((7:9),i);

        R = [cos(q(2))*cos(q(3)) cos(q(3))*sin(q(1))*sin(q(2))-cos(q(1))*sin(q(3)) sin(q(1))*sin(q(3))+cos(q(1))*cos(q(3))*sin(q(2));
             cos(q(2))*sin(q(3)) cos(q(1))*cos(q(3))+sin(q(1))*sin(q(2))*sin(q(3)) cos(q(1))*sin(q(2))*sin(q(3))-cos(q(3))*sin(q(1));
             -sin(q(2)) cos(q(2))*sin(q(1)) cos(q(1))*cos(q(2))];

        G = [0 -sin(q(3)) cos(q(2))*cos(q(3));
             0  cos(q(3)) cos(q(2))*sin(q(3));
             1 0 -sin(q(2))];

        qd = flip(R*inv(G)*(angVel - gb - gn));
        lvd = [0;0;-9.81] + R*(acc - ab - an);

        x = [pos; q; lv; gb; ab];
        sd = [pd*dt; qd*dt; lvd*dt; gbn; abn];
        nlf = x + sd;

        if (i == 1)
            Xt(:,i) = nlf;
        else
            Xt = horzcat(Xt,nlf);
        end
    end

    for x = 1:length(Xt)
        disp(x);
        
        if (x == 1)
            W = L / (n_prime + L);
        else
            W = 1/(2*(n_prime + L));
        end

        mean = W*Xt(:,x);
        final_mean = final_mean + mean;

    end

    for y = 1:length(Xt)
        
        if(y == 1)
            W = (L/(n_prime +L)) + (1-alpha^2+beta);
        else
            W = 1/(2*(n_prime + L));
        end
        covar = W * (Xt(:,y) - uPrev) * transpose(Xt(:,y) - uPrev);
        final_covar = final_covar + covar;
    end
    
    uEst = final_mean;
    covarEst = final_covar;

end

