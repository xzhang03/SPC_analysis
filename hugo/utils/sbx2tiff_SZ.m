function sbx2tiff_SZ(fpath, pmt, xybin, tbin, frame_info)
% Inputs
if nargin < 5
    frame_info = []; % [Start frame, N frames]
    if nargin < 4
        tbin = 1;
        if nargin < 3
            xybin = 1;
            if nargin < 2
                pmt = 0;
                if nargin < 1
                    [fn, fp] = uigetfile('\\nasquatch\data\2p\stephen\*.sbx');
                    fpath = fullfile(fp, fn);
                    disp(fn);
                end
            end
        end
    end
end

% Use ui to get file
if isempty(fpath)
    [fn, fp] = uigetfile('\\nasquatch\data\2p\stephen\*.sbx');
    fpath = fullfile(fp, fn);
    disp(fn);
end


% Read
fprintf('Reading PMT %i... ', pmt);
if isempty(frame_info)
    x = sbxReadPMT(fpath, 0, -1, pmt);
else
    x = sbxReadPMT(fpath, frame_info(1) - 1, frame_info(2) - frame_info(1) + 1, pmt);
end
fprintf('Done. \n')

% Bin
fprintf('Binning... ');
if xybin > 1
    x = binxy(x, xybin);
end
if tbin > 1
    x = bint(x, tbin);
end
fprintf('Done. \n');

% Output
if isempty(frame_info)
    spath = sprintf('%s_PMT%i_binxy%i_bint%i.tif', fpath(1:strfind(fpath,'.')-1), pmt, xybin, tbin);
else
    spath = ...
        sprintf('%s_PMT%i_binxy%i_bint%i_Frames%i-%i.tif', ...
        fpath(1:strfind(fpath,'.')-1), pmt, xybin, tbin, ...
        frame_info(1), frame_info(2));
end
writetiff(x, spath);

end