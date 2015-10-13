function [Gamma,Alpha,s,ch_Alpha,ch_s] = vMF(Z,C,Alpha,kappa,s,P,K,N)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = C./repmat(sqrt(sum(C.^2,2)),1,K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Posterior (E-step)
Gamma = repmat(repmat(s,[1,3]).*Alpha,[1,1,N]); % posterior mean of factor geolocations in time

for i=1:P    
    Mi = size(Z{i},2);    
    Zi = reshape(Z{i}(1:3,:),[1,3,Mi]);
    Ci = repmat(kappa(i)*C(i,:)',[1,3,Mi]);    
%     Ci = repmat(C(i,:)',[1,3,Mi]);    
    gamma = Ci.*repmat(Zi,[K,1,1]);
    hr_ch = [1,find(diff(Z{i}(4,:))>0)+1,Mi];
    n_hr = length(hr_ch)-1;
    for j=1:n_hr
        Gamma(:,:,Z{i}(4,hr_ch(j))) = Gamma(:,:,Z{i}(4,hr_ch(j))) + sum(gamma(:,:,hr_ch(j):hr_ch(j+1)-1),3);
    end
end

G_norm = sqrt(sum(Gamma.^2,2));
Gamma = Gamma./repmat(G_norm,[1,3,1]);

% Parameter Updates (M-step)
Alpha_old = Alpha;
s_old = s;

Alpha = sum(Gamma,3);   % factor geolocations averaged over time
a_norm = sqrt(sum(Alpha.^2,2));
Alpha = Alpha./repmat(a_norm,[1,3]);

r = a_norm/N;
s = (3*r-r.^3)./(1-r.^2);

ch_Alpha = norm(Alpha-Alpha_old,'fro')/norm(Alpha_old,'fro');
ch_s = norm(s-s_old)/norm(s_old);
