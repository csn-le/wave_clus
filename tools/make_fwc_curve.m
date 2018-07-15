function make_fwc_curve(input,par_input)
%Handler to call wave_features to draw feature stats curve
%Inputs:
% matrix of spikes
% name of file with spikes inside

if nargin<2
	par_input = set_parameters();
end
fname = {};

if isnumeric(input) && min(size(input))>1
    fname = {[]};
elseif ischar(input)
    fname = {input};    
elseif isnumeric(input)
    dirnames = [dir('times_*.mat');dir('*_spikes.mat')];
    dirnames = {dirnames(~[dirnames.isdir]).name};
    
    for f =1:length(dirnames)
        aux = regexp(dirnames{f}, '\d+', 'match');
        n = str2num(aux{1});
        if ismember(n,input)
            input(input==n) = [];
            fname{end+1} = dirnames{f};    
        end
    end
end


for f = 1:length(fname)
    filename = fname{f}; 
    if isempty(filename)
        filename = datestr(now,'dd_mm_yyyy_HH_MM_SS');
        spikes = input;
        par = par_input;
    else
        load(filename,'spikes')
        
        aux = who('-file',filename);
        if ismember('par',aux) && ~isempty(strfind(filename,'times_'))
            load(filename,'par')
        else
            par = par_input;
        end
    end
        
        
    par.filename = filename;
    par.plot_feature_stats = true;
    wave_features(spikes,par);
end
end