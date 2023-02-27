function t = spcMovtrace(mov, r, c, binsz)
% Get trace from a movie

%% Dimensions
% Limits
dims = size(mov);

% Rows
r1 = max(r - binsz, 1);
r2 = min(r + binsz, dims(1));

% Columns
c1 = max(c - binsz, 1);
c2 = min(c + binsz, dims(2));

% Output
t3 = mov(r1:r2, c1:c2, :);
t = squeeze(sum(sum(t3,1),2));

end