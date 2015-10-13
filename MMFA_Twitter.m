load Twitter_Multimodal_Aug14

K = 20;   % number of factors

[C,mu,Sigma,kappa,Alpha,s,Q,iter] = MultiModalFactorAnalysis(Y,Z,H,K,M);

nC = sqrt(sum(C.^2,2));
Cn = C./repmat(nC,1,K);
[Cns,Is] = sort(Cn,1,'descend');
Hs = Htag(Is);    % shows the dominant hashtags in the factors

% Optional Clustering
% L = 2*K;
% Idx = kmeans(Cn,L,'Distance','cosine','Replicates',10,'EmptyAction','drop','MaxIter',1000);
% Hc = cell(P,L);
% Ii = cell(L,1);
% for i=1:L
%     Ii{i} = find(Idx==i);
%     Hc(1:length(Ii{i}),i) = Htag(Ii{i});
% end

% Geolocations of factors
Alpha_sp = zeros(K,2);
Alpha_sp(:,1) = asin(Alpha(:,3));
As = asin(Alpha(:,2)./cos(Alpha_sp(:,1)));
At = atan(Alpha(:,2)./Alpha(:,1))/pi*180;
same = As./At<0;
Alpha_sp(same,2) = At(same);
Alpha_sp(~same,2) = At(~same)-180*sign(At(~same));
Alpha_sp(:,1) = Alpha_sp(:,1)/pi*180;
