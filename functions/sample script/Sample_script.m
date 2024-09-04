% keep sigstruct
clear
mouse = 'SZ1133';
date = 240823;
runs = [2];
user = 'stephen';
cdigit = 2;


%% Check gross tm
spcPeeksdt(mouse, date, runs, 'user', user, 'cdigit', cdigit, 'crop',...
    [], 'bint', 2, 'autosave', true, 'makeplot', true, 'phthresh', [0 9000]);

%% Register from smcs
for runs = 1 : 2
    spcSMCReg(mouse, date, runs, 'user', user, 'force', false, 'binxy', 2, 'cdigit', cdigit,...
        'iterations', 3, 'previewlocalnorm', false, 'crop', [226 11 1027 510], 'outputbinxy', 2);
end

%% Demonsreg from smcs
for runs = 1
    spcSMCDemreg(mouse, date, runs, 'user', user, 'force', false, 'binxy', 2, 'cdigit', cdigit,...
        'iterations', 1, 'outputbinxy', 2, 'savesmc', true);
end

%% Prep for cellpose
spcCellposePrep(mouse, date, runs, 'user', user, 'sourcetype', 'smcdemreg', 'uselocalnorm', true,...
    'hp_norm_sigmas', [16 64], 'previewlocalnorm', true)
    
%% Cellpose
return;

%% Apply cellpose ROIs
spcApplyCellposeROI_smc(mouse, date, runs, 'user', 'stephen', 'cdigit', 2, 'GRIN', true, ...
    'StartFrameN', 1, 'binxy', 2, 'force', false, 'npsize', [30 2]);
