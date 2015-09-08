classdef nse_wc_reader < handle
	properties
        raw_filename
    end 
	methods 
        function obj = ncs_wc_reader(par, raw_filename)
            obj.raw_filename = raw_filename;
        end
        
        function [sr,max_segments,with_raw,with_spikes] = get_info(obj)
        	sr = [];
            max_segments = 0;
            with_raw = false;
            with_spikes = true;
        end
        
        function [spikes, index] = load_spikes(obj)
            [index, Samples] = Nlx2MatSE(obj.raw_filename,1,0,0,0,1,0);
            spikes(:,:)= Samples(:,1,:);
            clear Samples
            spikes = spikes';
            index = index/1000;
        end

    end
end
