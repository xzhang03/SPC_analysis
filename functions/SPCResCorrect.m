function SPCResCorrect(mouse, date, varargin)
% SPCResCorrect corrects the warping artifact caused by resonant mirrors.
% This is a wrap function to usee Hugo's xRegCorrect function.

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'user', ''); % user name for path
addOptional(p, 'slice', false); % Flag if data is slice
addOptional(p, 'sbxrun', []); % A single run number of the sbx file, used as reference (
addOptional(p, 'flimruns', []); % Run numbers of the flim files to warp

% Start with registered or unregistered data
addOptional(p, 'useregistered', true);

% Sbx parameters
addOptional(p, 'sbxframeinfo', [50 500]); % First and total number of frames to calculate reference image
addOptional(p, 'operator', @median); % What functions to use for calculating reference
addOptional(p, 'savesbxref', true); % Save sbx reference file

% Force saving
addOptional(p, 'force', false);

% Manual point sampling (shouldn't have to)
addOptional(p, 'getNewPoints', false);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% sbx reference
% Output paths
% run string
if p.sbxrun < 10
    runstr = sprintf('00%i', p.sbxrun);
elseif p.sbxrun < 100
    runstr = sprintf('0%i', p.sbxrun);
else
    runstr = num2str(p.sbxrun);
end

% Folders
datemouse = sprintf('%i_%s', date, mouse);
datemouserun = sprintf('%i_%s_run%s', date, mouse, runstr);
sbxrefpath = sprintf('\\\\%s\\data\\2p\\%s\\%s\\%s\\%s\\',...
p.server, p.user, mouse, datemouse, datemouserun);

% Filenames
sbxfn = sprintf('%i_%s_run%s.sbx', date, mouse, runstr);
sbxrefname = sprintf('%i_%s_run%s_ref.tif', date, mouse, runstr);



end