function Plot_all_features(handles)
% function Plot_all_features(handles)

USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
inspk = USER_DATA{7};
classes = USER_DATA{6};

colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 

maxc = size(colors,1);

figure(10)
nclasses = max(classes);
inputs = size(inspk,2);
for i=1:inputs
    for j=i+1:inputs
        %subplottight(inputs,inputs,(i-1)*inputs+j)
        subplottight(inputs-1,inputs-1,(i-1)*(inputs-1)+j-1)
        
        hold on
        for k=1:nclasses
            class_aux = find(classes==k);
            max_spikes = min(par.max_spikes_plot,length(class_aux));
            plot(inspk(class_aux(1:max_spikes),i),inspk(class_aux(1:max_spikes),j),'.','color', colors(mod(k-1,maxc)+1,:),'markersize',.5)
            axis off
        end
    end
end

    

function h = subplottight(n,m,i)
    [c,r] = ind2sub([m n], i);
    ax = subplot('Position', [(c-1)/m, 1-(r)/n, 1/m, 1/n]);
    if(nargout > 0)
      h = ax;
    end