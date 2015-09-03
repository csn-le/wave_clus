classdef mat_reader < handle
	properties
        sr
        data
        max_segments
        opened_file
        sr_infile
        t0_segments
        segmentLength
        spikes_file
    end 
	methods 
        function obj = mat_reader(par, raw_filename)
            obj.sr = [];
            obj.spikes = [];
            obj.data = [];
            finfo = whos('-file', raw_filename);
            
            if ismember('sr',{finfo.name})  %if is possible, load sr from file; else from set_parameters
                load(raw_filename,'sr'); 
                obj.sr_infile = true;
                obj.sr = sr;
            else
                obj.sr_infile = false;
                obj.sr = par.sr; 
            end
            
            if ismember('spikes',{finfo.name}) && ismember('index',{finfo.name})
               obj.spikes_file = true;
            end
            if ismember('data',{finfo.name})
                load(raw_filename,'data');

                if strcmp(par.tmax,'all')
                    t0 = 0;
                    obj.data = data;
                else
                    min_index = floor(par.tmin * obj.sr);                   %min time to read (in micro-sec)
                    max_index = ceil(par.tmax * obj.sr); 
                    t0 = (min_index-1)/obj.sr*1000;
                    obj.data = data(min_index:max_index);
                end
                n = max(size(obj.data));

                %Segments the data in par.segments pieces
                obj.max_segments = ceil(n/ obj.sr / ...
                        (par.segments_length * 60));         %number of segments in which data is cutted


                obj.segmentLength = floor (n/obj.max_segments);

                obj.t0_segments = ones(1,obj.max_segments);
                obj.t0_segments(1) = t0;
                for i = 2:obj.max_segments
                    obj.t0_segments(i) = obj.t0_segments(i-1) + obj.segmentLength/obj.sr*1000;
                end
            end
                
                        
        end

        function [spikes, index] = load_spikes(obj)
             load(raw_filename, spikes); 
             load(raw_filename,index); 
        end
        
        function [sr,max_segments,with_raw,with_spikes] = get_info(obj)
            if isempty(obj.data)
                with_raw = false;
                max_segments = 0;
            else
                with_raw = true;
                max_segments = obj.max_segments;
            end
            
            with_spikes = obj.spikes_file;

            if obj.sr_infile
                sr = obj.sr;
            else
                sr = [];
            end
        end
        
        function index_ts = index2ts(obj,index,i)
            index_ts = (index-1)/obj.sr*1000 + obj.t0_segments(i);
        end
      
        function x = get_segment(obj,i)
            if i ~= obj.max_segments
                x = obj.data(obj.segmentLength*(i-1)+1:obj.segmentLength*i);
            else
                x = obj.data(obj.segmentLength*(i-1)+1:end);
            end
        end
    end
end
