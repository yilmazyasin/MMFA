function [U,Theta,Q,Eta,ch_Q,ch_Eta] = Multinomial(H,C,Q,Eta,P,K,D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = C./repmat(sqrt(sum(C.^2,2)),1,K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Posterior (E-step)
% Posterior Covariance of Latent Factor Scores (over words)
while rank(Q)<K
    Q = Q + 1e-10*diag(diag(Q));
end
Qi = Q\eye(K);
Cc = C'*C;
tDp1 = 2*(D+1);

CcQi = Cc/2+Qi;
if sum(sum((isnan(CcQi)|isinf(CcQi))))>0
    keyboard
end
while rank(CcQi)<K
    CcQi = CcQi + 1e-10*diag(diag(CcQi));
end
PHI_ti = CcQi\eye(K);

CcQi = Cc/tDp1+Qi;
if sum(sum((isnan(CcQi)|isinf(CcQi))))>0
    keyboard
end
while rank(CcQi)<K
    CcQi = CcQi + 1e-10*diag(diag(CcQi));
end
Delta = -( eye(P) + (D/tDp1)*C*( CcQi\eye(K) )*C' )/tDp1;   % P in notes
Delta = PHI_ti*C'*Delta*C*PHI_ti;

% Posterior Mean of Latent Factor Scores
exp_nor = 0;
oflow = max(max(Eta));
if oflow>=700    % to avoid exp(710)=inf
%     exp_nor = floor(oflow/700)*700;
    exp_nor = oflow;
    Eta = Eta-exp_nor;
end
Exp_eta = exp(Eta);
P_eta = Exp_eta./repmat(sum(Exp_eta,2)/exp(exp_nor),1,D);
Prob_err = H./repmat(sum(H,2),1,D) - P_eta;
% Ht = 2*(Prob_err+repmat(sum(Prob_err),D,1)) + Eta;
% Ht_s = sum(Ht);
% G = Ht'/2 - repmat(Ht_s'/tDp1,1,D);
G = Prob_err + Eta/2 - repmat(sum(Eta,2)/tDp1,1,D); % PxD (G' in notes)
Phi = PHI_ti*C'*G + repmat(Delta*C'*sum(G,2),1,D); % columns are mean vectors for each word over factors

% Parameter Updates (M-step)
Q_old = Q;
Eta_old = Eta;
% R = D*(Delta+PHI_ti) + Phi*Phi';
% Q = R/D;
Q = PHI_ti + Delta + Phi*Phi'/D;
Eta = C*Phi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Eta = Eta./repmat(sqrt(sum(Eta.^2,2)),1,D)*50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eta = randn(D,P);
ch_Q = norm(Q-Q_old,'fro')/norm(Q_old,'fro');
ch_Eta = norm(Eta-Eta_old,'fro')/norm(Eta_old,'fro');

% for updating Factor Coefficients
Phi_s = sum(Phi,2);
% U = R/2 - (D^2*Delta+D*PHI_ti + Phi_s*Phi_s')/tDp1;
U = (D^2/tDp1)*PHI_ti + (D/tDp1)*Delta + Phi*Phi'/2 - Phi_s*Phi_s'/tDp1;
% Theta = Phi*repmat(Ht_s,D,1);
Theta = Phi*G';
