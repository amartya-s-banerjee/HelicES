function [lambda0_mnk, mprime_vect, nprime_vect] = MATVEC_prereq(S)

lambda0_mnk = zeros(S.Num_Basis, 1);
mprime_vect = zeros(S.Num_Basis, 1);
nprime_vect = zeros(S.Num_Basis, 1);

inc = 1;

for m = -S.M_max:S.M_max
    for n = -S.N_max:S.N_max        
        abs_val_sq = abs(2*pi*(m-S.alpha*n*S.N)*(1/S.Tau))^2;       
        for k = 1:S.K_max
            
            mprime_vect(inc) = 1i*2*pi*m;
            nprime_vect(inc) = 1i*2*pi*n*S.N;
%             nprime_vect(inc) = n*S.N;
            
            lambda0_mnk(inc) = (Get_Bessel_Roots(n, k, S)/S.R)^2 + ...
                abs_val_sq;
            
            inc = inc + 1;
        end
    end
end

end