function [num_k_pts, weights_k_pts, kptlist] = ...
    Create_List_and_Weights_of_k_points(S)

% List of k-points
h_nkpt = S.Input_num_k_pts; % number of k-points in the theta1 direction
num_k_pts = S.N*h_nkpt; % total number of k-points
kptlist = zeros(num_k_pts,2); % list of k-points

count = 1;

for i=1:S.N
    for j=1:h_nkpt
        kptlist(count,:) = [(2*j-h_nkpt-1)/(2*h_nkpt), i-1];
        count = count + 1;
    end
end

weights_k_pts = ones(num_k_pts,1)/num_k_pts; % weights for the k-points
cgo = S.N;

% Fold the nu values to the negative side : needed for the spectral scheme
nu_lim = fix(S.N/2);
for ii = 1:num_k_pts
    nu = kptlist(ii,2);
    if(nu > nu_lim)
        nu = nu_lim - nu;
    end
    kptlist(ii,2) = nu;
end



% reduction by time-reversal symmetry
if S.Flag_Time_Reversal_Sym == 1
    Iw = zeros(num_k_pts,1);
    Irem = zeros(num_k_pts,1);
    for p=1:num_k_pts
        for q=p+1:num_k_pts
            % % if abs(kptlist(qq,1)+kptlist(pp,1))<1e-6 && ...
            % % abs(kptlist(qq,2)+kptlist(pp,2)-cgo+1)<1e-6 
            if abs(kptlist(q,1)+kptlist(p,1))<1e-6 && ...
                    (abs(kptlist(q,2)+kptlist(p,2)-cgo)<1e-6 || ...
                    abs(kptlist(q,2)+kptlist(p,2))<1e-6)
                Iw(p) = 1;
                Irem(q) = 1;
            end
        end
    end
    
    Iw = Iw > 0.5;
    Irem = Irem > 0.5;
    weights_k_pts(Iw) = 2*weights_k_pts(Iw);
    kptlist = kptlist(~Irem,:);
    weights_k_pts = weights_k_pts(~Irem,:);
    num_k_pts = size(kptlist,1);
end

end