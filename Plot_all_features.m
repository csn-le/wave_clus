function Plot_all_features(handles)
% function Plot_all_features(handles)

USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
inspk = USER_DATA{7};
classes = USER_DATA{6};

colors = ['b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];

figure(10)
nclasses = max(classes);
inputs = size(inspk,2);
for i=1:inputs
    for j=i+1:inputs
        subplot(inputs,inputs,(i-1)*inputs+j)
        hold on
        for k=1:nclasses
            class_aux = find(classes==k);
            max_spikes = min(par.max_spikes,length(class_aux));
            plot(inspk(class_aux(1:max_spikes),i),inspk(class_aux(1:max_spikes),j),['.' colors(k)],'markersize',.5)
            axis off
        end
    end
end

    

