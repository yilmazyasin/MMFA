% function [C,d,mu,Sigma,kappa,alpha,s,Q,Eta] = initialize_MMFA(P,K,D)
function [C,Alpha,s,Q,Eta] = initialize_MMFA_vMF_Mult(P,K,D)

C = randn(P,K);     % factor loadings
% C = C./repmat(sqrt(sum(C.^2,2)),1,K);
% C = repmat(1/K,P,K);

Alpha = randn(K,3);     % vMF parameters
Alpha = Alpha./repmat(sqrt(sum(Alpha.^2,2)),1,3);
s = rand(K,1);

Q = eye(K);         % Multinomial parameters
% Taylor series expansion point for the log-sum-exp function
Eta = randn(P,D)*100;   
