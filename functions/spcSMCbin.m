function smc = spcSMCbin(smc, xybin, tbin)
% Bin smc files in xy and in t dimensions. 

if nargin < 3
    % Default no time binning
    tbin = 1;    
    if nargin  < 2
        % Default 2x xy bin
        xybin = 2;
    end
end
sz = smc(1).size;
sz_new = floor([sz(1)/xybin, sz(2)/xybin, sz(3)/tbin, sz(4), 0]);
n = size(smc,1);

if tbin == 1 && xybin == 1
    % Nothing to do
    return;
end

%% Binxy
if xybin > 1
    for i = 1 : n
        % Update dimensions
        smc(i).size(1) = sz_new(1);
        smc(i).size(2) = sz_new(2);
            
        if isempty(smc(i).r)
            continue;
        end

        % Binning
        r = ceil(smc(i).r / xybin);
        c = ceil(smc(i).c / xybin);
        li = r + (c-1) * sz_new(1);
        v = full(smc(i).v) + 1;
        
        % Find ones to merge
        [~, ind_in, ind_out] = unique(li, 'stable');
        
        % New values
        r_new = r(ind_in);
        c_new = c(ind_in);
        v_new = v(ind_in);
        
        % Find repeaters
        inds_old = (1 : length(r))';
        repeaters = inds_old(~ismember(inds_old,ind_in));
        
        for tind_old = repeaters'
            tind_new = ind_out(tind_old);
            v_new(tind_new) = v_new(tind_new) + v(tind_old);
        end
        
        % Debug
        if sum(v) ~= sum(v_new)
            fprintf('Error at T = %i, expected %i pixels and saw %i\n', i, sum(v), sum(v_new))
        end
        
        % Putback
        % Sparse matrix if it ends here
        if tbin == 1
            smc(i).v = sparse(v_new-1);
        else
            smc(i).v = v_new;
        end
        smc(i).r = uint16(r_new);
        smc(i).c = uint16(c_new);
        smc(i).size(5) = sum(v_new);
    end
end


%% Bin t
if tbin > 1
    keepvec = zeros(sz(3) * sz(4), 1);
    for ind2 = 1 : sz_new(4)
        for ind_new = 1 : sz_new(3)
            % Start and end row numbers
            istart = (ind2-1) * sz(3) + (ind_new - 1) * tbin + 1;
            iend = (ind2-1) * sz(3) + ind_new * tbin;
            
            % Update dimensions and tracker
            smc(istart).size(3) = sz_new(3);
            keepvec(istart) = 1;
            
            % New
            r = cat(1,smc(istart:iend).r);
            c = cat(1,smc(istart:iend).c);
            v = cat(1,smc(istart:iend).v);
            
            if xybin == 1
                % Unwrap if needed
                v = full(v) + 1;
            end
            
            % If empty nothing to do
            if isempty(r)
                continue;
            end
            
            % Linear indexing
            li = r + (c-1) * sz_new(1);
            
            % Find ones to merge
            [~, ind_in, ind_out] = unique(li, 'stable');
            
            % New values
            r_new = r(ind_in);
            c_new = c(ind_in);
            v_new = v(ind_in);
            
            % Find repeaters
            inds_old = (1 : length(r))';
            repeaters = inds_old(~ismember(inds_old,ind_in));

            for tind_old = repeaters'
                tind_new = ind_out(tind_old);
                v_new(tind_new) = v_new(tind_new) + v(tind_old);
            end
            
            % Debug
            if sum(v) ~= sum(v_new)
                fprintf('Error at T = %i, expected %i pixels and saw %i\n', i, sum(v), sum(v_new))
            end
            
            % Put back
            smc(istart).v = sparse(v_new-1);
            smc(istart).r = uint16(r_new);
            smc(istart).c = uint16(c_new);
            smc(istart).size(5) = sum(v_new);
        end
    end
    smc = smc(keepvec > 0);
end


end