% function [C,d,mu,Sigma,kappa,Alpha,s,Q,iter] = MultiModalFactorAnalysis(Y,Z,H,K,M)
function [C,Alpha,s,Q,iter] = MultiModalFactorAnalysis_vMF_Mult(Z,H,K,N,kappa)

% Inputs:
%   Y:      Hashtag counts - Poisson (PxN matrix of P hashtags in N hours)
%   Z:      Spherical coordinates (x,y,z) - von Mises-Fisher (Px1 cell. Each cell is a M_ix3 matrix.)
%   H:      Word counts - Multinomial (DxP matrix of P hashtags for D words)
%   K:      number of factors
%   M:      maximum number of coordinates among P hashtags

%% Initialize parameters
P = length(Z);
D = size(H,2);
D = D-1;    % set the last word as pivot
H(:,end) = [];
[C,Alpha,s,Q,Eta] = initialize_MMFA_vMF_Mult(P,K,D);
thr = 1e-2;
change = 1;
iter = 0;

%% EM Algorithm
while change && iter<100
    iter = iter+1;
    fprintf('ITERATION %d \n',iter)
    % vMF
    tic        
    [Gamma,Alpha,s,ch_Alpha,ch_s] = vMF(Z,C,Alpha,kappa,s,P,K,N);
    fprintf('vMF: %f sec, ch_Alpha=%f, ch_s=%f \n',toc,ch_Alpha,ch_s) 
    % Multinomial
    tic
    [U,Theta,Q,Eta,ch_Q,ch_Eta] = Multinomial(H,C,Q,Eta,P,K,D);
    fprintf('Multinomial: %f sec, ch_Q=%f, ch_Eta=%f \n',toc,ch_Q,ch_Eta)
    % Update Factor Loadings    
    tic
    [C,ch_C] = Factor_Loadings_vMF_Mult(Z,C,kappa,Gamma,P,K,N,U,Theta);
    fprintf('Factor Loadings: %f sec, ch_C=%f \n',toc,ch_C)
    change = ch_Alpha>thr | ch_s>thr | ch_Q>thr | ch_Eta>thr*10 | ch_C>thr*2;
end
