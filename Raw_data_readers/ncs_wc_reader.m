classdef ncs_wc_reader < handle
	properties
        sr
        max_segments
        raw_filename
        %channel
        opened_file
        recmin
        recmax
        num_scale_factor
        TimeStamps
        dt
    end 
	methods 
        function obj = ncs_wc_reader(par, raw_filename)
            obj.raw_filename = raw_filename;            
            obj.opened_file = fopen(obj.raw_filename, 'r', 'l');
            fseek(obj.opened_file,16384,'bof');                       %Skip Header, put pointer to the first record
            TimeStamps = fread(obj.opened_file,inf,'int64',(4+4+4+2*512));          %Read all TimeStamps, in us
            fseek(obj.opened_file, 16384+8+4+4+4, 'bof');                             %put pointer to the beginning of data
            
                                 
            % conversion to datapoints
            tsdiff = diff(TimeStamps);
            obj.dt = median(tsdiff);
            min_dt = min(diff(TimeStamps));
            if min_dt<=0, %check for corrupt TimeStamps that would violate the monotonic increase
                warning('corrupt TimeStamps - attempting correction')
                ind = find(tsdiff <= 0); %corrupt TimeStamps are usually too low
                TimeStamps(ind+1) = TimeStamps(ind) + obj.dt;
                clear tsdiff;
                min_dt = min(diff(TimeStamps));
                if min_dt <= 0
                    error('Automated correction of corrupt TimeStamps failed - try manually!');
                end
            end
            
            obj.sr = 512*1e6/obj.dt;            % sampling rate (in Hz).
            time0 = TimeStamps(1); 
            
            
            if strcmp(par.tmax,'all')
                index_tinitial = 0;
                index_tfinal = length(TimeStamps);
                tsmax = TimeStamps(end);
                tsmin = time0;
            else
                tsmin = time0 + par.tmin*1e6;                   %min time to read (in micro-sec)
                tsmax = time0 + par.tmax*1e6;                   %max time to read (in micro-sec)  
                index_tinitial = find(tsmin > TimeStamps);
                if isempty(index_tinitial) == 1;
                    index_tinitial = 0;
                else
                    index_tinitial = index_tinitial(end)-1;
                end    
                index_tfinal = find(tsmax < TimeStamps);
                if isempty(index_tfinal) ==1;
                    index_tfinal = length(TimeStamps);
                else
                    index_tfinal = index_tfinal(1);
                end 
            end
            fseek(obj.opened_file,16384+8+4+4+4+index_tinitial*(8+4+4+4+2*512),'bof');

            lts = index_tfinal - index_tinitial;
            %Segments the data in par.segments pieces
            obj.max_segments = ceil((tsmax - tsmin)/ ...
                    (par.segments_length * 1e6 * 60));         %number of segments in which data is cutted
            segmentLength = floor (lts/obj.max_segments);

            tsmin = 1 : segmentLength :lts;
            tsmin = tsmin(1:obj.max_segments);
            tsmax = tsmin - 1;
            tsmax = tsmax (2:end);
            tsmax = [tsmax, lts];
            obj.recmax = tsmax;
            obj.recmin = tsmin;
            obj.TimeStamps = TimeStamps;
            scale_factor = textread(obj.raw_filename,'%s',43);

            if(str2num(scale_factor{41})*1e6 > 0.5)
                obj.num_scale_factor = 1e6 * str2num(scale_factor{43}); %for the new CSC format
            else
                obj.num_scale_factor = 1e6 * str2num(scale_factor{41}); %for the old CSC format
            end
            
        end
        
        function [sr,max_segments,with_raw,with_spikes] = get_info(obj)
        	sr = obj.sr;
            max_segments = obj.max_segments;
            with_raw = true;
            with_spikes = false;
        end
        
        function index_ts = index2ts(obj,index,i)
            index_ts = (obj.TimeStamps(obj.recmin(i)+floor(index/512))'+mod(index,512)/512*obj.dt)/1000;
        end
      
        function x = get_segment(obj,i)
            
            Samples = fread(obj.opened_file,512*(obj.recmax(i)- ...
                obj.recmin(i)+1),'512*int16=>int16',8+4+4+4); %the recmax are in timestamps indexs
            
            x = double(Samples(:))' *  obj.num_scale_factor;
           
            if i == obj.max_segments
                fclose(obj.opened_file);
            end
        end
    end
end
