
function new_par = update_parameters(new_par, par, type,NaN_fill)
% This function overwrite or create fields from par in new_par. The fields
%   used will be only the ones with type 'type' present in 'par'.
% new_par a struct() could be empty.
% par a struct() with parameters.
% type a string: 'detect' for parameters used in detection an alignment of
%   spikes. 'clus' for parameters used in detection or clustering. 'relevant' for both.

if ~exist('NaN_fill','var')
    NaN_fill = false;
end
detection_params = {'channels','segments_length', 'sr','tmax','tmin','w_pre', ...
    'w_post','alignment_window', 'stdmin','stdmax', 'detect_fmin', ...
    'detect_fmax','sort_fmin','sort_fmax', 'ref_ms', 'detection', ...
    'int_factor','interpolation','alignment_window','sort_order','detect_order','detection_date'};



clus_params = {'min_inputs','max_inputs','scales','features','template_sdnum', 'template_k', ...
    'template_k_min','template_type','force_feature','match', ...
    'max_spk','permut','mintemp', 'maxtemp', 'tempstep','SWCycles',...
    'KNearNeighb', 'min_clus','max_clus','randomseed','force_auto','sorting_date','elbow_min','c_ov'};

batch_ploting_params = {'temp_plot','max_spikes_plot','print2file','plot_feature_stats'};


new_par_names = fieldnames(new_par);
load_par_names = fieldnames(par);
if strcmp(type,'detect') || strcmp(type,'relevant')
    for i= 1:length(detection_params)
        if ismember(detection_params(i),load_par_names)
            field = char(detection_params(i));
            new_par.(field ) = par.(field );
         elseif NaN_fill
            field = char(detection_params(i));
            new_par.(field) = NaN;
        end
    end
end

if strcmp(type,'batch_plot')
    for i= 1:length(batch_ploting_params)
        if ismember(batch_ploting_params(i),load_par_names)
            field = char(batch_ploting_params(i));
            new_par.(field ) = par.(field );
        end
    end
end

if strcmp(type,'clus') || strcmp(type,'relevant')
    for i= 1:length(clus_params)
        if ismember(clus_params(i),load_par_names)
            field = char(clus_params(i));
            new_par.(field ) = par.(field );
        end
    end
end

if strcmp(type,'none')
    for i= 1:length(new_par_names)
        if ~ (any((ismember(new_par_names(i),clus_params)) || any(ismember(new_par_names(i),detection_params))))
            field = char(new_par_names(i));
            if ismember(new_par_names(i),load_par_names)
                new_par.(field) = par.(field );
            else
                new_par.(field) = [];
            end
        end
    end

end


if isfield(par,'ref')
    new_par.ref_ms = par.ref /par.sr *1000;
end

end
