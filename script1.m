
%----------------------------------------------------------------
% File:     script1.m
%----------------------------------------------------------------
%
% Author:   Marek Rychlik (rychlik@arizona.edu)
% Date:     Thu Apr 13 11:38:17 2023
% Copying:  (C) Marek Rychlik, 2020. All rights reserved.
% 
%----------------------------------------------------------------
% Arithmetic encoding of Markov Chains

% Make a count matrix (would be a probability matrix if rows divided by sum)
% NOTE: Sum of row = 1, P(i,j) = P(i->j)
P = [1 2 3 4; 5 6 7 8; 1 1 2 2; 1 3 4 2];
P0 = P./sum(P,2);

s = RandStream('mlfg6331_64');


n = 10000;
x = zeros(1,n);
x(1) = 1;

% Generate a sequence
for j = 1:(n-1)
    x(j+1) = datasample(s, 1:4, 1 ,'Weights',P0(x(j),:));
end

x_last = x(end);
y={};
for u=1:4
    if u == x_last
        idx = find(x(1:n-1)==u)+1;
    else
        idx = find(x==u)+1;
    end
    y{u} = x(idx);
    code{u} = arithenco(y{u}, P(u,:));
end

y_dec={};
for u=1:4
    y_dec{u} = arithdeco(code{u}, P(u,:), numel(y{u}));
end

% Decode the original sequence
cnt = ones(1,4);
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
[V,D,W] = eig(P0);
Q = -P0 .* log2(P0);
mu = W(:,1)';                            % L^2 unit vector
mu = mu./sum(mu);                       % Make it L^1 unit vector;
Hexact = sum(mu*Q)