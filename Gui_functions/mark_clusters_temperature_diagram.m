% MARK CLUSTERS IN TEMPERATURE DIAGRAM

function mark_clusters_temperature_diagram(handles,tree,clustering_results)

handles.par.min.clus = clustering_results(1,5);

% creates cluster-temperature vector to plot in the temperature diagram
nclasses = max(clustering_results(:,2));
class_plot = [];
for i=1:nclasses
    ind = find(clustering_results(:,2)==i);
    classgui_plot(i) = clustering_results(ind(1),2);
    class_plot(i) = clustering_results(ind(1),4);
    if class_plot(i) == 0 %null original cluster
		class_plot(i) =1; %plot like they were from cluster 1
    end
    temp_plot(i) = clustering_results(ind(1),3);  
end

num_temp = floor((handles.par.maxtemp ... 
-handles.par.mintemp)/handles.par.tempstep);     % total number of temperatures 

tree(num_temp+1,2) = handles.par.mintemp+(num_temp)*handles.par.tempstep; %added for handle selection of max temp

temperature = tree(clustering_results(1,1)+1,2);

colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 
maxc = size(colors,1);
auto_sort_info = getappdata(handles.temperature_plot,'auto_sort_info');

% draw temperature diagram and mark clusters 
cla(handles.temperature_plot);
switch handles.par.temp_plot
    case 'lin'
        % draw diagram
        hold(handles.temperature_plot, 'on');
        if ~isempty(auto_sort_info)
            [xp, yp] = find(auto_sort_info.peaks);
            ptemps = handles.par.mintemp+(xp)*handles.par.tempstep;
            psize = tree(sub2ind(size(tree), xp,yp+4));
            plot(handles.temperature_plot,ptemps,psize,'xk','MarkerSize',7,'LineWidth',0.9);
            area(handles.temperature_plot,handles.par.mintemp+handles.par.tempstep.*[auto_sort_info.elbow,num_temp],max(ylim(handles.temperature_plot)).*[1 1],'LineStyle','none','FaceColor',[0.9 0.9 0.9]);
        end
        plot(handles.temperature_plot, [handles.par.mintemp handles.par.maxtemp-handles.par.tempstep],[handles.par.min.clus2 handles.par.min.clus2],'k:',...
            handles.par.mintemp+(1:num_temp)*handles.par.tempstep, ...
            tree(1:num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
        % mark clusters
        for i=1:min(size(tree,2)-4,length(class_plot))
            tree_clus = tree(temp_plot(i),4+class_plot(i));
            tree_temp = tree(temp_plot(i)+1,2);
            plot(handles.temperature_plot, tree_temp,tree_clus,'.','color',colors(mod(classgui_plot(i)-1,maxc)+1,:),'MarkerSize',20);
        end
        set(get(gca,'ylabel'),'vertical','Baseline');
    case 'log'
        % draw diagram
        set(handles.temperature_plot,'yscale','log');
        hold(handles.temperature_plot, 'on');
        if ~isempty(auto_sort_info)
            [xp, yp] = find(auto_sort_info.peaks);
            ptemps = handles.par.mintemp+(xp)*handles.par.tempstep;
            psize = tree(sub2ind(size(tree), xp,yp+4));
            semilogy(handles.temperature_plot,ptemps,psize,'xk','MarkerSize',7,'LineWidth',0.9);
            area(handles.temperature_plot,handles.par.mintemp+handles.par.tempstep.*[auto_sort_info.elbow,num_temp],max(ylim(handles.temperature_plot)).*[1 1],'LineStyle','none','FaceColor',[0.9 0.9 0.9],'basevalue',1);
        end
        semilogy(handles.temperature_plot, [handles.par.mintemp handles.par.maxtemp-handles.par.tempstep], ...
            [handles.par.min.clus handles.par.min.clus],'k:',...
            handles.par.mintemp+(1:num_temp)*handles.par.tempstep, ...
            tree(1:num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
        % mark clusters
        for i=1:length(class_plot)
            if class_plot(i)+4>size(tree,2)
                continue
            end
            tree_clus = tree(temp_plot(i),4+class_plot(i));
            tree_temp = tree(temp_plot(i)+1,2);
            semilogy(handles.temperature_plot, tree_temp,tree_clus,'.','color',colors(mod(classgui_plot(i)-1,maxc)+1,:),'MarkerSize',20);
        end
        set(get(handles.temperature_plot,'ylabel'),'vertical','Baseline');
end

% xlim(handles.temperature_plot, [0 handles.par.maxtemp])
xlabel(handles.temperature_plot, 'Temperature','FontSize',8); 
ylabel(handles.temperature_plot, 'Clusters size','FontSize',8);