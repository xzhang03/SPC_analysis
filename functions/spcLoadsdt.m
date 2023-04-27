function [mov, params] = spcLoadsdt(fpath, varargin)
% Load sdt file and align pixels. This version does everything from scratch
% [mov, params] = spcLoadsdt(fpath, varargin)
% Shoutout to Arnold Estrada, whose code is used here

%% Parse inputs
if nargin < 2
    varargin = {};
end
p = inputParser;

% General variables
addOptional(p, 'datatype', 'uint16');

% Frame variables
addOptional(p, 'tbins', 256);
addOptional(p, 'framesize', [1024 1024]);

% Rearrange
addOptional(p, 'transpose', true);
addOptional(p, 'fliplr', false);
addOptional(p, 'flipud', false);

% Crop
addOptional(p, 'crop', [512 1250]); % Leave empty to use full dimention (512 x 2048)

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Read headers
fp=fopen(fpath, 'rb', 'ieee-le');

% From Arnold Estrada
% https://github.com/aestrada71/Arnold_Estrada_Matlab_Code/blob/master/read_sdt-1.m
header.revision = fread(fp,1,'short');
header.info_offs = fread(fp,1,'long');
header.info_length = fread(fp,1,'short');
header.setup_offs = fread(fp,1,'long');
header.setup_length = fread(fp,1,'short');
header.data_block_offs = fread(fp,1,'long');
header.no_of_data_blocks = fread(fp,1,'short');
header.data_block_length = fread(fp,1,'long');
header.meas_desc_block_offs = fread(fp,1,'long');
header.no_of_meas_desc_blocks = fread(fp,1,'short');
header.meas_desc_block_length = fread(fp,1,'short');
header.header_valid = fread(fp,1,'ushort');
header.reserved1 = fread(fp,1,'ulong');
header.reserved2 = fread(fp,1,'ushort');
header.chksum = fread(fp,1,'ushort');

% check header
if header.header_valid ~= hex2dec('5555')
    fclose(fp);
    error('Invalid file header');
end

% Check frame size
typec1 = header.data_block_length / p.tbins / p.framesize(1) / p.framesize(2);
switch p.datatype
    case 'uint32'
        typec2 = 4;
    case 'uint16'
        typec2 = 2; 
    case 'uint8'
        typec2 = 1; 
end
if typec1 ~= typec2
    fclose(fp);
    error('Invalid frame sizes, must add up to %i pixels\n', header.data_block_length);
end

%% data block header
fseek(fp, header.data_block_offs, 'bof');
dataheader.block_no = fread(fp, 1, 'short');
dataheader.data_offs = fread(fp, 1, 'long');
dataheader.next_block_offs = fread(fp, 1, 'long');
dataheader.block_type = fread(fp, 1, 'ushort=>ushort');
dataheader.meas_desc_block_no = fread(fp, 1, 'short');
dataheader.lblock_no = fread(fp, 1, 'ulong');
dataheader.block_length = fread(fp, 1, 'ulong');

%% Read data block
fseek(fp, dataheader.data_offs, 'bof');
switch p.datatype
    case 'uint32'
        I = fread(fp, header.data_block_length/4, 'uint32'); 
    case 'uint16'
        I = fread(fp, header.data_block_length/2, 'uint16'); 
    case 'uint8'
        I = fread(fp, header.data_block_length, 'uint8'); 
end

%% Rearrange
% Get in shape
I2 = reshape(I, [p.tbins, p.framesize(1), p.framesize(2)]);

% Get the frames in the right order
for i = 1 : p.tbins
    % Get frame
    f = I2(i,:,:);
    f = squeeze(f);

    % Transpose and flip
    if p.transpose
        f = f';
    end
    if p.fliplr
        f = flip(f, 2);
    end
    if p.flipud
        f = flip(f, 1);
    end

    % Save frame
    if i == 1
        sz = size(f);
        mov = uint16(zeros(sz(1), sz(2), p.tbins));
    end
    mov(:,:,i) = f;
end

% Upwrap (f^*%ing hell who saves data this way man)
mov = cat(2, mov(1:2:end,:,:), mov(2:2:end, :, :));

% Crop
if ~isempty(p.crop)
    mov = mov(1:p.crop(1), 1:p.crop(2), :);
end

%% Save parameters
params = struct('header', header, 'dataheader', dataheader, 'p', p);

end