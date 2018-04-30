new_folder = [];
max_memo_GB = [];
% This code requires the file read_Intan_RHD2000_file.m which is provided
% on the Intan Technologies web site.
% new_folder: output folder with the parsed channels and metadata. 
% If it doen't exist or is [], it will be 'wc_data'.
% max_memo_GB is an idea of the number of GB allocated for the data to be
% stored in RAM, so it is used to compute the number of segments in which
% the data should be split for processing
% Based in read_intan_data.m version 1.1, June 26, 2010
% (c) 2010, Intan Technologies, LLC
% For more information, see http://www.intantech.com
% For updates and latest version, see http://www.intantech.com/software.html

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


read_Intan_RHD2000_file;

sr = frequency_parameters.amplifier_sample_rate;
num_channels = length(amplifier_channels); % info sale de la header file

if ~exist('new_folder','var') || isempty(new_folder)
    new_folder = [pwd filesep 'wc_data'];
end

mkdir(new_folder)

% Read 'info.rhd'



% Read 'amplifier.dat'
fileinfo = dir([path filesep 'amplifier.dat']);
lts = fileinfo.bytes/(num_channels * 2); % int16 = 2 bytes
samples_per_channel = ceil(max_memo/(num_channels * 4)); %single (4 bytes)
num_segments = ceil(lts/samples_per_channel);
channels = 1:num_channels;
save([new_folder filesep 'intan_meta_data'],'sr','lts','channels')


fid = fopen([path filesep 'amplifier.dat'], 'r');

outfile_handles = cell(1,num_channels);

[fname] = 'ch';
for i = 1:num_channels
    outfile_handles{i} = fopen([new_folder filesep  fname '_' num2str(channels(i)) '.intch'],'w','l');  
end

for j=1:num_segments
    ini = (j-1)*samples_per_channel;
    fin = min(j*samples_per_channel,lts);
    data2 = single(fread(fid,(fin-ini)*(num_channels),'int16')) * 0.195;
    for ind = 1:num_channels
        fwrite(outfile_handles{ind},data2(ind:num_channels:end),'single');
    end
    fprintf('Segment %d out of %d processed.\n',j,num_segments);
end
fclose('all');
clearvars
%Auxiliar para gui_plot
%time = num_samples / sr; %total time amplifier
%j=-1;
%save('trial','j')