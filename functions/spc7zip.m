function spc7zip(mice, ziplevel, varargin)
% spc7zip uses 7zip to compress files and backup
% Zip level can be 0 (no compression), 1, 3, 5, 7, 9 (ultra). Leave empty for default 5
% There seems to be <5% difference in file sizes between 5 (normal) and 1 (fasted)
% or 9 (ultra). The fasted option seems to be 2.5x faster than normal

%% Parse inputs
if nargin < 3
    varargin = {};
    if nargin < 2
        ziplevel = [];
        if nargin < 1
            mice = {};
        end
    end
end

p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', 'stephen'); % user name for path

% 7zip variables
addOptional(p, 'allowupdate', true);
addOptional(p, 'zippath', 'C:\Program Files\7-Zip\7z.exe'); % 7Zip path
addOptional(p, 'targetpath', 'D:\User Folders\Stephen\Deep_backup\2p'); % Target path

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs
% No mice
% Target dir
if ~exist(fullfile(p.targetpath, 'Backup_tracker.mat'), 'file')
    donelist = dir(fullfile(p.targetpath, '*.7z'));
else
    donelist = load(fullfile(p.targetpath, 'Backup_tracker.mat'));
    donelist = donelist.donelist;
end
if isempty(mice)
    mice = uigetmice();
end
    
% Mice
if ~iscell(mice)
    mice = {mice};
end
nmice = length(mice);

% Check zip
if ~exist(p.zippath, 'file')
    disp('7z.exe not found. Please check path.')
    return;
end

%% Construct command
for imice = 1 : nmice
    % Case and type
    mouse = mice{imice};
    mouse = upper(mouse);

    % User (add yourself if needed)
    if isempty(p.user)
        switch mouse(1:2)
            case 'SZ'
                p.user = 'stephen';
            otherwise
                p.user = 'stephen';  
        end
    end

    % Input fp
    fp_in = sprintf('\\\\%s\\data\\2p\\%s\\%s', p.server, p.user, mouse);

    % Output fp
    fp_out = fullfile(p.targetpath, sprintf('%s.7z', mouse));
    if exist(fp_out, 'file')
        mode = 'a';
    else
        mode = 'u';
    end

    % Base command
    basecom = sprintf('"%s" %s "%s" "%s"', p.zippath, mode, fp_out, fp_in);

    %% Add switches
    % Compression level
    if ~isempty(ziplevel)
        basecom = sprintf('%s -mx%i', basecom, ziplevel);
    end

    %% Execute
    % Run
    if exist(fp_out, 'file')
        if p.allowupdate
            fprintf('%i/%i Updating %s...', imice, nmice, mouse);
        else
            fprintf('%i/%i Skipping %s...', imice, nmice, mouse);
            continue;
        end
    else
        fprintf('%i/%i Compressing %s...', imice, nmice, mouse);
    end
    tic;
    [status, result] = system(basecom);
    t = toc/60;
    fprintf(' Done (%0.1f min)\n', t);

    % In case something went wrong
    if status ~= 0
        disp('Something went wrong:')
        disp(result);
        return;
    end
    
    % Add to tracker
    donelist_new = dir(fp_out);
    irenew = strcmp({donelist(:).name}, donelist_new.name);
    if any(irenew)
        donelist(irenew) = donelist_new;
    else
        donelist(end+1) = donelist_new; %#ok<AGROW>
    end
    save(fullfile(p.targetpath, 'Backup_tracker.mat'), 'donelist', '-v7.3');
    
    % Parse message
    i1 = regexp(result, 'bytes (');
    i2 = regexp(result, ' GiB)');
    i3 = regexp(result, ' MiB)');
    
    % Sizes
    try
        s1 = str2double(result(i1(1)+7:i2(1)-1));
    catch
        try
            s1 = str2double(result(i1(1)+7:i3(1)-1)) / 1024;
        catch
            s1 = 1e-3;
        end
    end
    try
        s2 = str2double(result(i1(3)+7:i2(3)-1));
    catch
        try
            s2 = str2double(result(i1(3)+7:i3(3)-1)) / 1024;
        catch
            s2 = 1e-3;
        end
    end
    fprintf('%0.1f GB -> %0.1f GB (-%0.1f%%, %0.1f GB/min)\n', s1, s2, (s1-s2)/s1*100, s1/t);
    
end

%% Anonymouse functions
function mice = uigetmice()
    % base dir
    basedir = sprintf('\\\\%s\\data\\2p\\%s', p.server, p.user);
    micelist = dir(basedir);
    micelist = micelist(3:end);
    micelist = {micelist(:).name};
    
    % Add done
    micelist2 = micelist;    
    for i = 1 : length(donelist)
        idone = strcmpi(micelist, donelist(i).name(1:end-3));
        micelist2{idone} = sprintf('%s (done)', micelist{idone});
    end
    
    % Choose
    s = listdlg('ListString', micelist2, 'PromptString', 'Choose mice');
    mice = micelist(s);
end
end