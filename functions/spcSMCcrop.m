function smc = spcSMCcrop(smc, crop)
% Crop smc as if it's a stack
% crop is the top-left and bottom-right coordinates of the bounding box
% [Xtl, Ytl, Xbr, Ybr]. X is the column number and Y is the row number.

% SMC cell
if iscell(smc)
    smccell = smc;
    smc = smccell{1};
    ncell = length(smccell);
else
    ncell = 1;
end

% Get sizes
sz = smc(1).size;
sz_new = [crop(4)-crop(2)+1, crop(3)-crop(1)+1, sz(3), sz(4), 0];
n = size(smc,1);


%% Loop through and crop
for ii = 1 : ncell
    if ncell > 1
        smc = smccell{ii};
    end
    
    for i = 1 : n
        % Update dimention
        smc(i).size(1) = sz_new(1);
        smc(i).size(2) = sz_new(2);

        if isempty(smc(i).r)
            continue;
        end

        % Load
        r = smc(i).r;
        c = smc(i).c;
        v = full(smc(i).v) + 1;

        % Keep
        keeps = (r >= crop(2)) & (r <= crop(4)) & (c >= crop(1)) & (c <= crop(3));

        % Subtract
        r = r(keeps) - crop(2) + 1;
        c = c(keeps) - crop(1) + 1;
        v = v(keeps);

        % Putback
        smc(i).r = uint16(r);
        smc(i).c = uint16(c);
        smc(i).v = sparse(double(v-1));
    end
    
    if ncell > 1
        smccell{ii} = smc;
    end
end

% SMC cell put back
if ncell > 1
    smc = smccell;
end

end