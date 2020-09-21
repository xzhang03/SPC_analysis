%% Dataset
inputset = {'SZ333', 200302, 1; 'SZ333', 200305, 1; 'SZ333', 200311, 1;...
    'SZ333', 200312, 1; 'SZ333', 200314, 1; 'SZ333', 200315, 1;...
    'SZ333', 200317, 1; 'SZ334', 200302, 1; 'SZ334', 200306, 1;...
    'SZ334', 200311, 1; 'SZ334', 200312, 1; 'SZ334', 200314, 1;...
    'SZ334', 200315, 1; 'SZ335', 200303, 1; 'SZ335', 200307, 1;...
    'SZ335', 200311, 1; 'SZ335', 200312, 1; 'SZ335', 200314, 1;...
    'SZ335', 200315, 1};
inputset = inputset(6,:);
N = size(inputset,1);

%% Load
% Initialize
datastruct = struct('Expt', [],'Photons', [], 'Tm', [], 'NPPhotons', [], 'NPTm', []);
datastruct = repmat(datastruct, [N,1]);

for ind = 1 : N
    % Path
    spcpaths = spcPath(inputset{ind,1}, inputset{ind,2}, inputset{ind,3},...
        'server', 'nasquatch', 'user', 'stephen', 'slice', false, 'cdigit', 1);
    
    % Loading
    ROI_struct = load(fullfile(spcpaths.fp_out, spcpaths.xrun_mat));
    ROI_struct = ROI_struct.ROI_struct(inputset{ind,3});
    
    % Filling
    datastruct(ind).Expt = ind;
    datastruct(ind).Photons = ROI_struct.photon_vec';
    datastruct(ind).Tm = ROI_struct.tm_vec';  
    datastruct(ind).NPPhotons = ROI_struct.photon_np_vec';
    datastruct(ind).NPTm = ROI_struct.tm_np_vec';  
end

%% Combine
NP_photons = [datastruct(:).NPPhotons];
NP_tm = [datastruct(:).NPTm];
photons = [datastruct(:).Photons];
tm = [datastruct(:).Tm];
scatter(NP_photons, NP_tm)

%% Bin
% Trim
n = randi(length(NP_tm), [280 1]);
NP_photons2 = NP_photons(n);
NP_tm2 = NP_tm(n);
photons2 = photons(n);
tm2 = tm(n);

% Sort
[NP_photons2, NPorder] = sort(NP_photons2, 'ascend');
NP_tm2 = NP_tm2(NPorder);
[photons2, SMorder] = sort(photons2, 'ascend');
tm2 = tm2(SMorder);

% Reshape
NP_photons2 = reshape(NP_photons2, [20 14]);
NP_tm2 = reshape(NP_tm2, [20 14]);
photons2 = reshape(photons2, [20 14]);
tm2 = reshape(tm2, [20 14]);

% Plot
hold on
scatter(median(photons2), median(tm2,1))
scatter(median(NP_photons2), median(NP_tm2,1))
hold off