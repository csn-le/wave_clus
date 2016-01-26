function parse_data_pl2(filename)

% Requires Matlab Offline SDK (Plexon Inc) in the Matlab searchpath
% Based in code written by Azriel Ghadooshahy, Jan 2016

%Gather relevant info
pl2 = PL2GetFileIndex(filename);
[unused, fname, ext] = fileparts(filename);

%List of possible different continuous sources
%query source-channels for whether or not recording was enabled
source_names_cont = cellfun(@(x) x.SourceName,pl2.AnalogChannels,'uniformoutput',0); %it's possible that they haven't the same size
cont_record_chs = cellfun(@(x) x.RecordingEnabled,pl2.AnalogChannels,'uniformoutput',1);


%Load PL2 recorded continuous data
source_string = 'WB'; %could also be 'SPKC' or 'FP'
inds_in_source = ismember(source_names_cont,source_string);
inds_in_source = find(inds_in_source & cont_record_chs);
sourcenums = cellfun(@(x) x.Channel,pl2.AnalogChannels(inds_in_source),'uniformoutput',1);
%cdat = cell(1,numel(sourcenums));
channels = sourcenums;
fprintf('%d channels in the file\n', numel(sourcenums) )

for n = 1:numel(sourcenums) 
    num_in_source = sourcenums(n); 
    handler = PL2AdBySource(filename,source_string,num_in_source);
    lts = handler.FragCounts;
    sr = handler.ADFreq;
    fout = fopen([fname '_' num2str(num_in_source) '.pl2ch'],'w');
    fwrite(fout,handler.Values,'double');
    fclose(fout);
    fprintf('Channel %d out. \n', num_in_source);
end

save('pl2_meta_data','sr','lts','channels')
