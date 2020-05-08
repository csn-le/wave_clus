%Class for handle 'nc5' files. Wave_clus uses the firsts letters of the m file to now what files it can open.
%

%The classes used by readInData need the methods:
%   - A constructor, with inputs par (parameters) and 
%     raw_filename (name of the file to read).
%   - A get_info() method that return: [sr, max_segments, with_raw, with_spikes]
%       sr : sampling rate in Hz. A empty vector if the sr in the parameters
%            should be use.
%       max_segments: total number of segments of raw data to read. Numbers
%                   of times that set_segment(i) will be call.
%       with_raw: Boolean, true if the raw file has a continouos signal.
%       with_spikes: Boolean, true if the raw file has only the detected spikes.
%   - A get_segment(i) method, that should return the continuous segment
%       number i. It will be call in order from i=1 to i=max_segment (if it's necessary).
%   - A index2ts(index,i) method that return a vector of times in ms,
%       with the times of the 'index' elements of the segment i.

classdef nc5_wc_reader < handle
	properties      % Wave_clus will not access any internal property of the class
        sr
        max_segments
        opened_file
        segmentLength
        open_file
        t0_segments
    end 
	methods 
        function obj = nc5_wc_reader(par, raw_filename)
            stamps_name = [raw_filename(1:regexp(raw_filename,'_\d*.NC5')) '_TimeStamps.mat'];
            if exist(stamps_name,'file')
                load(stamps_name ,'lts', 'sr');
            else
                load('NSX_TimeStamps','lts', 'sr');
            end
            obj.sr = sr;
            initial_index = floor(par.tmin  * obj.sr);                   %min time to read (in sec)
            
            if strcmp(par.tmax,'all')
                tmax = lts/sr;
            else
                tmax = min(par.tmax,lts/sr);
            end
            obj.max_segments = ceil((tmax - par.tmin)/ ...
                    (par.segments_length *60));         %number of segments in which data is cutted
            obj.opened_file = fopen(raw_filename,'r','l');
			fseek(obj.opened_file,initial_index*2,'bof');
            
            obj.segmentLength = zeros(1,obj.max_segments);
            slen = floor(par.segments_length * 60 * sr);
            obj.segmentLength(1) = min(floor((tmax - par.tmin) * sr),slen);
            obj.t0_segments = zeros(1,obj.max_segments);
            obj.t0_segments(1) = initial_index / obj.sr *1e3;
            for i = 2:obj.max_segments
            	obj.t0_segments(i) = obj.t0_segments(i-1) + obj.segmentLength(i-1)/obj.sr*1000;
                obj.segmentLength(i) = min(floor((tmax - obj.t0_segments(i)/1000) * sr),slen);
            end

        end
        
        function [sr,max_segments,with_raw,with_spikes] = get_info(obj)
        	sr = obj.sr;
            max_segments = obj.max_segments;
            with_raw = true;
            with_spikes = false;
        end
        
        function index_ts = index2ts(obj,index,i)
            index_ts = (index)/obj.sr*1000 + obj.t0_segments(i);
        end
      
        function x = get_segment(obj,i)
            x=fread(obj.opened_file,obj.segmentLength(i),'int16=>double')'/4;
            if i == obj.max_segments
                fclose(obj.opened_file);
            end
        end
    end
end
