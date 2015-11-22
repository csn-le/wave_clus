function formats = supported_wc_extensions()
    % Return a cell of arrays with the formats supported. Automatically updated
    formats ={};
    raw_readers_folder = fileparts(mfilename('fullpath'));
    files = dir(raw_readers_folder);
    for i = 1: length(files)
        fname = files(i).name;
        fend = regexp(fname,'_wc_reader.m');
        if ~isempty (fend)
             formats = [formats,{fname(1:fend-1)}];
        end
    end
end