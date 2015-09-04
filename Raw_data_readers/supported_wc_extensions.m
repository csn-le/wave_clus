function formats = supported_wc_extensions()
    % Return a cell of arrays with the formats supported. Automatically updated
    formats ={};
    raw_readers_folder = fileparts(mfilename('fullpath'));
    fnames = ls(raw_readers_folder);
    for i = 1: size(fnames,1)
        fname = fnames(i,:);
        fend = regexp(fname,'_wc_reader.m');
        if ~isempty (fend)
             formats = [formats,{fname(1:fend-1)}];
        end
    end
end