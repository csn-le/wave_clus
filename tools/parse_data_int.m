function parse_data_int(filename, max_memo_GB)
% max_memo_GB is an idea of the number of GB allocated for the data to be
% stored in RAM, so it is used to compute the number of segments in which
% the data should be split for processing
% Based in read_intan_data.m version 1.1, June 26, 2010
% (c) 2010, Intan Technologies, LLC
% For more information, see http://www.intantech.com
% For updates and latest version, see http://www.intantech.com/software.html


[unused, fname, ext] = fileparts(filename);
ext = lower(ext(2:end));

if strcmp(ext,'.int')
    error('Incorrect extension.');
end

with_memory=true;
try
	memory;
catch
	with_memory=false;
end
if with_memory
	[userview,systemview] = memory;
	memo_avaible = floor(systemview.PhysicalMemory.Available*0.80);
	if exist('max_memo_GB','var') && ~isempty(max_memo_GB)
        max_memo = max_memo_GB*(1024)^3;
		if max_memo > memo_avaible
			error('max_memo_GB > 80% of Physical Memory Available')
		end
	else
		max_memo = memo_avaible;
	end
else
	max_memo = max_memo_GB*(1024)^3;
end


fid = fopen(filename, 'r');

% Read first three header bytes encoding file version
for i=1:3
    header(i) = fread(fid, 1, 'uint8');
end

if (header(1) ~= 128)
    error('Improper data file format.');
end

if (header(2) ~= 1 || header(3) ~= 1)
    warning('Data file version may not be compatible with this m-file.');
end

% Now see which amplifier channels are saved in this file.

for i=1:64
    amp_on(i) = fread(fid, 1, 'uint8');
end

num_amps = sum(amp_on);
channels = find(amp_on);
s = dir(filename);
filesize = s.bytes;

lts = (filesize - 67)/(num_amps*4 + 1);

sr = 25000;
save('intan_meta_data','sr','lts','channels')

%-----------------------------------

fprintf(1, '\nData file contains %0.2f seconds of data from %d amplifier channels.\n', lts/sr, num_amps);
fprintf(1, 'Channels: %s \n',num2str(channels));

%--------------------------------------


% Go back to the beginning of the file...
frewind(fid);
% ...skip the header this time...
fread(fid, 3+64, 'uint8');

samples_per_channel = ceil(max_memo/(num_amps*4+1));
num_segments = ceil(lts/samples_per_channel);

outfile_handles = cell(1,num_amps);

for i = 1:num_amps
    outfile_handles{i} = fopen([fname '_' num2str(channels(i)) '.intch'],'w','l');  
end


for j=1:num_segments
    ini = (j-1)*samples_per_channel;
    fin = min(j*samples_per_channel,lts);
    data2 = fread(fid,(fin-ini)*(num_amps*4+1),'uint8=>uint8');
    data2((num_amps*4)+1:num_amps*4+1:end) = [];
    data2 = typecast(data2,'single');
    for ind = 1:num_amps
        fwrite(outfile_handles{ind},data2(ind:num_amps:end),'single');
    end
    fprintf('Segment %d out of %d processed.\n',j,num_segments);
end

fclose('all');
