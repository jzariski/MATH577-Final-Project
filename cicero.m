%----------------------------------------------------------------
% File:     cicero.m
%----------------------------------------------------------------
%
% Author:   Marek Rychlik (rychlik@arizona.edu)
% Date:     Sat Apr 22 12:31:27 2023
% Copying:  (C) Marek Rychlik, 2020. All rights reserved.
% 
%----------------------------------------------------------------
% Basic marginalization estimates for the Markov model

f=fopen('cicero.txt','r');
X=fread(f,inf,'uint8');
% Recode as an ASCII string
txt=char(X)';
fclose(f);

L = numel(txt);
k = 4;

% States are blocks of length k, so record them
A = zeros(L-k, k, 'uint8');
for j=1:L-k+1
    A(j,:) = X(j:j+k-1);
end
% ... also returns index vectors IA and IC such that C = A(IA,:) and A = C(IC,:)
[C,IA,IC] = unique(A, 'rows');

%tbl=histcounts(IC,1:size(C,1)+1);
%bar(tbl);

% The number of states in the Markov chain
N = size(C,1);

% Estimate counts of transitions
S = sparse(N,N);

% Initial distribution
I = zeros(N,1);
% Record transitions from block of length k+1
for j=1:L-k;
    from = X(j:j+k-1)';
    to = X(j+1:j+k)';
    [lia, locb] = ismember([from; to], C, 'rows');
    S(locb(1),locb(2)) = S(locb(1),locb(2)) + 1;
    % Emissions
    E(locb(1)) = to(end);
    if(j==1)
        % Record the initial state
        I(locb(1)) = 1;
    end
end

r = amd(S);

% Transition probability matrix
esttr = S./(eps+sum(S,2));
estemit = E./sum(E);

%... [V,D,W] = EIG(A) also produces a full matrix W whose columns are the
%    corresponding left eigenvectors so that W'*A = D*W'.

% Use sparse version of eig, which yields 6 top eigenvalues/vectors
% Thus, we must transpose.
[V,D] = eigs(esttr',15);
%[V,D] = eig(full(esttr)');
% Modulus of eigenvalues
D=diag(D);
Da = abs(D);
spectral_gap=1-max(Da(Da<.99));

[~,idx] = max(Da);
mu = V(:,idx)';                            % L^2 unit vector
mu = mu./sum(mu);                       % Make it L^1 unit vector;

% Find the entropy
H = sum(-mu*(esttr.*log2(eps+esttr)));

% Compute the initial state
idx = find(I);
initstate = char(C(idx,:));

t = tiledlayout('flow');

nexttile(t);
spy(S,'*');
title('Transition matrix sparsity pattern')

d = symrcm(S);
nexttile(t);
spy(S(d,d),'ro')
title('S(d,d) After Cuthill-McKee Ordering')

%----------------------------------------------------------------
%
% Reorder states to get convergence of eigenvalues
%
%----------------------------------------------------------------

nexttile(t);
r = amd(S);
spy(S(r,r))
title('S(r,r) After Minimum Degree Ordering')
% Now the sparse eigenvalues algorithm converges
P=S(r,r);
P=P./(eps+sum(P,2));
[V1,D1] = eigs(P',20);
if ~any(isnan(D1),'all')
    disp('The sparse eigenvalues algorithm converges')
end
D1 = diag(D1);
D1 = D1(~isnan(D1));
crcl = exp(i*linspace(0,2*pi,30))';
nexttile(t);
plot(D1,'bs');
hold on;
plot(crcl,'g.-');
hold off;