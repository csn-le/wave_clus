classdef nc5_wc_reader < handle
	properties
        sr
        max_segments
        %channel
        opened_file
        segmentLength
        TimeStamps
        dt
        open_file
    end 
	methods 
        function obj = ncs_wc_reader(par, raw_filename)
            load('NSX_TimeStamps','lts', 'sr');

            obj.sr = sr;
            obj.raw_filename = raw_filename;
            
            obj.open_file = fopen(raw_filename,'r','l');
            
            if strcmp(par.tmax,'all')
                initial_index = 0;
                obj.max_segments = ceil(lts/(par.segments_length  * sr * 60)); %number of segments in which data is cut
            else
                initial_index = floor(par.tmin * obj.sr);                   %min time to read (in micro-sec)
                tsmax = min(par.tsmax,lts/sr);
                obj.max_segments = ceil((tsmax - par.tsmin)/ ...
                    (par.segments_length *60));         %number of segments in which data is cutted
            end
            
            obj.opened_file = fopen(raw_filename,'r','l');
			fseek(obj.opened_file,initial_index*2,'bof');
            
            obj.segmentLength = floor (lts/obj.max_segments);
             
            obj.t0_segments = zeros(1,obj.max_segments);
            obj.t0_segments(1) = t0;
            for i = 2:obj.max_segments
            	obj.t0_segments(i) = obj.t0_segments(i-1) + obj.segmentLength/obj.sr*1000;
            end

        end
        
        function [sr,max_segments,with_raw,with_spikes] = get_info(obj)
        	sr = obj.sr;
            max_segments = obj.max_segments;
            with_raw = true;
            with_spikes = false;
        end
        
        function index_ts = index2ts(obj,index,i)
            index_ts = (index-1)/obj.sr*1000 + obj.t0_segments(i);
        end
      
        function x = get_segment(obj,i)
            x=fread(f1,obj.segmentLength,'int16=>double')/4;
            if i == obj.max_segments
                fclose(obj.opened_file);
            end
        end
    end
end
