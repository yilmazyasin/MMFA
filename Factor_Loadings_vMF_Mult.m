function [C,ch_C] = Factor_Loadings_vMF_Mult(Z,C,kappa,Gamma,P,K,N,U,Theta)

C_old = C;
Gamma = mat2cell(Gamma,K,3,ones(1,N));
Gamma = reshape(Gamma,1,N);
% lam = 1e-3;
lamq = 1e0;
lamV = 1e0;
% lamM = max(max(U.^-1));
lamM = 1e0;
Hess = lamM*U+lamq; % -Hessian  
for i=1:P     
        Mi = size(Z{i},2);        
        Zi = mat2cell(Z{i}(1:3,:),3,ones(1,Mi)); 
        Gi = Gamma(Z{i}(4,:));
        Gz = cellfun(@mtimes,Gi,Zi,'UniformOutput',false);
        clear Zi Gi
        Gz = sum(squeeze(cell2mat(Gz)),2);        
        Grad = lamV*kappa(i)*Gz + lamM*Theta(:,i);% - lam*sign(C(i,:))';  % Grad-Hess          
        
        if sum(sum((isnan(Hess)|isinf(Hess))))>0
            keyboard
        end
        while rank(Hess)<K
            Hess = Hess + 1e-10*diag(diag(Hess));
        end
        C(i,:) = Hess\Grad;     
end

% C = C./repmat(sqrt(sum(C.^2,2)),1,K);
ch_C = norm(C-C_old,'fro')/norm(C_old,'fro');
% nC = norm(C,'fro');
% if nC>P*K*1e6
%     C = C/nC*P*K*1e5;
% end
