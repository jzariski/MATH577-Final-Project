f=fopen('coriolan.txt','r');
X=fread(f,inf,'uint8');
% Recode as an ASCII string
txt=char(X)';
fclose(f);

L = numel(txt);
k = 1;

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
S = ones(N,N); %% Needs to be non-sparse matrix

% Initial distribution
I = zeros(N,1);
% Record transitions from block of length k+1
L-k
for j=1:L-k
    from = X(j:j+k-1)';
    to = X(j+1:j+k)';
    [lia, locb] = ismember([from; to], C, 'rows');
    S(locb(1),locb(2)) = S(locb(1),locb(2)) + 1;
    % Emissions
    % CHANGED THIS PART
    rnum = rand;
    E(locb(1), from(1)) = rnum;
    if(j==1)
        % Record the initial state
        I(locb(1)) = 1;
    end
    j
end

r = amd(S);
P = S;
% Transition probability matrix
S = 100*S+1; % Need to start with non-zero matrix
esttr = S./(eps+sum(S,2));
E = E+1e-3; % Need to start with non-zero matrix
estemit = E./sum(E,2);

% ADDED SEQ HERE
seq = A(:,1)';

%%
P0 = P./sum(P,2);

s = RandStream('mlfg6331_64');


n = 100000;
x = zeros(1,n);
x(1) = 1;

%% Testing

% Generate a sequence
for j = 1:(n-1)
    x(j+1) = datasample(s, 1:N, 1 ,'Weights',esttr(x(j),:));
end
[C, ia, ic] = unique(x);
%tbl = histcounts(ic,1:size(C,2)+1);

x_last = x(end);
n = length(x);
y={};

for u=1:N
    if u == x_last
        idx = find(x(1:n-1)==u)+1;
    else
        idx = find(x==u)+1;
    end
    y{u} = x(idx);
    code{u} = arithenco(y{u}, P(u,:));
end

y_dec={};
for u=1:N
    count = numel(y{u});
    if count==0 % Added to make sure arithendeco works
        count=1;
    end
    y_dec{u} = arithdeco(code{u}, P(u,:), count);
end

% Decode the original sequence
cnt = ones(1,N);
xx = zeros(1,n);
xx(1) = x(1);                           % Must be sent
for j=1:(n-1)
    xx(j+1) = y_dec{xx(j)}(cnt(xx(j))); cnt(xx(j)) = cnt(xx(j)) + 1;
end

% Validate decoding
assert(all(x==xx));

% bits per symbol
num_code_bits = numel([code{:}]);

% Heuristic entropy
H = num_code_bits./numel(x)

% Exact entropy of the Markov Chain
[V,D,W] = eig(esttr);
Q = -esttr .* log2(esttr);
Q(isnan(Q))=0;
mu = W(:,1)';                            % L^2 unit vector
mu = mu./sum(mu);                       % Make it L^1 unit vector;
Hexact = sum(mu*Q)


