function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state

    alpha = 0.001;
    k = 1;
    beta = 2;
    R = diag([1e-2 1e-2 1e-2]);
    
    Rbc = [ 0.7071 -0.7071   0
           -0.7071 -0.7071   0
                 0       0  -1];
    
    T = [-0.04; 0; -0.03];
    S = [    0  0.03    0;
         -0.03     0 0.04;
             0 -0.04    0];

    Zt = zeros(3,1);
    final_zut = zeros(3,1);
    final_Ct = zeros(15,3);
    final_St = zeros(3,3);

    uAug = uEst;
    n_prime = length(uAug);
    L = alpha^2 *(n_prime + k) - n_prime;
    
    covarAug = chol(covarEst,"lower");
    sigma_point(:,1) = uAug;

    for i = 1:n_prime
        xAug_odd = uAug + (sqrt(n_prime + L))*covarAug(:,i);
        xAug_even= uAug - (sqrt(n_prime + L))*covarAug(:,i);
        xAug = [xAug_odd xAug_even];
       
        sigma_point = horzcat(sigma_point,xAug); 
    end

    for j = 1: length(sigma_point)
        q = sigma_point((4:6),j);
        lv =sigma_point((7:9),j);

        R_bw = eul2rotm(transpose(q), "ZYX");

        l = inv(R_bw) * lv;
        nlf = Rbc*l - Rbc*S * inv(Rbc)*z_t(4:6);

        if (j == 1)
            Zt(:,1) = nlf;
        else 
            Zt = horzcat(Zt,nlf);
        end 
    end

     for x = 1:length(Zt)
        
         if (x == 1)
            Wm = L / (n_prime + L);
        else
            Wm = 1/(2*(n_prime + L));
         end

        z_ut = Wm*Zt(:,x);
        final_zut = final_zut + z_ut;
     end
    


     for y = 1:length(Zt)
        
        if(y == 1)
            Wc = (L/(n_prime +L)) + (1 - (alpha^2) + beta);
        else
            Wc = 1/(2*(n_prime + L));
        end

        Ct = Wc*(sigma_point(:,y) - uEst)*(transpose(Zt(:,y) - final_zut));
        final_Ct = final_Ct + Ct;
     end

     

     for z = 1:length(Zt)
        
        if(z == 1)
            Wc = (L/(n_prime +L)) + (1- (alpha^2) + beta);
        else
            Wc = 1/(2*(n_prime + L));
        end

        St = Wc*(Zt(:,z) - final_zut)*(transpose(Zt(:,z) - final_zut));
        final_St = final_St + St;
     end

     

     Kt = final_Ct * inv(final_St + R);
     uCurr = uEst + Kt*(z_t(1:3) - final_zut);
     covar_curr = covarEst - Kt*final_St*transpose(Kt);
     
  
end

