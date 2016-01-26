function plot_spikes(handles)
set(handles.file_name,'string','Plotting...'); drawnow;
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
spikes = USER_DATA{2};
classes = USER_DATA{6};
classes = classes(:)';
inspk = USER_DATA{7};
temp = USER_DATA{8};
ls = size(spikes,2);
minclus = handles.minclus;
clustering_results = USER_DATA{10};

% Closes aux figures
h_figs=get(0,'children');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2= findobj(h_figs,'tag','wave_clus_aux1');
h_fig3= findobj(h_figs,'tag','wave_clus_aux2');
h_fig4= findobj(h_figs,'tag','wave_clus_aux3');
h_fig5= findobj(h_figs,'tag','wave_clus_aux4');
h_fig6= findobj(h_figs,'tag','wave_clus_aux5');
close(h_fig1); close(h_fig2); close(h_fig3); close(h_fig4); close(h_fig5); close(h_fig6);
if ishandle(10)
	close(10)
end
% Extract spike features if needed
if get(handles.spike_shapes_button,'value') ==0
    if isempty(inspk) || (length(inspk)~=size(spikes,1))
        [inspk] = wave_features(spikes,handles);
        USER_DATA{7} = inspk;
    end
end

% Classes should be consecutive numbers
classes_names = sort(unique(classes));
classes_names = classes_names(classes_names>0);

% updates 'clustering_results_bk'
clustering_results_bk = clustering_results;

% Forcing
if handles.force==1
    for i = classes_names
        ind = find(clustering_results_bk(:,2)==i); % get index of GUI class
        oclass = clustering_results_bk(ind(1),4); % get original class
        otemp = clustering_results_bk(ind(1),3); % get original temperature
        ind2 = find(classes==i); % get index of forced class
        clustering_results(ind2,3) = otemp; % update original temperatures with forced class
        clustering_results(ind2,4) = oclass; % update original class with forced class
    end
end



for i = 1:min(length(classes_names),par.max_clus)
   c = classes_names(i);
   if c~= i
       classes(classes == c) = i;
   end
end

% Defines nclusters
cluster_sizes = zeros(1,par.max_clus);
ifixflag = zeros(1,par.max_clus);
for i=1:par.max_clus
    cluster_sizes(i) = nnz(classes==i);
end

if handles.setclus == 0 && handles.undo==0 && handles.merge==0 && handles.force==0  
    sizemin_clus = par.min_clus;
elseif handles.setclus == 1
	sizemin_clus = 1;
else
    sizemin_clus = par.min_clus;
end

clusn = find(cluster_sizes >= sizemin_clus);
nclusters = length(clusn);

% Get fixed clusters
fix_class2 = [];
nfix_class = [];

if ~isfield(handles,'new_manual')
	if get(handles.fix1_button,'value') ==1     
		nclusters = nclusters +1;
		fix_class = USER_DATA{20}';
		classes(classes==nclusters)=0;
		classes(fix_class)=nclusters;
		ifixflag(nclusters)=1;
		
		fix_class2 = [fix_class2 fix_class]; 
		nfix_class = [nfix_class 1];
		clusn = [clusn nclusters];
	end
	if get(handles.fix2_button,'value') ==1     
		nclusters = nclusters +1;
		fix_class = USER_DATA{21}';
		classes(classes==nclusters)=0;
		classes(fix_class)=nclusters;
		ifixflag(nclusters)=1;
		
		fix_class2 = [fix_class2 fix_class]; 
		nfix_class = [nfix_class 2];
		clusn = [clusn nclusters];
	end
	if get(handles.fix3_button,'value') ==1     
		nclusters = nclusters +1;
		fix_class = USER_DATA{22}';
		classes(classes==nclusters)=0;
		classes(fix_class)=nclusters;
		ifixflag(nclusters)=1;
		
		fix_class2 = [fix_class2 fix_class]; 
		nfix_class = [nfix_class 3];
		clusn = [clusn nclusters];
	end
	% Get fixed clusters from aux figures
	for i=4:par.max_clus
		eval(['fixx = par.fix' num2str(i) ';']);
		if fixx == 1 && ~isempty(USER_DATA{22+i-3})
			nclusters = nclusters +1;
			fix_class = USER_DATA{22+i-3}';
			classes(classes==nclusters) = 0;
			classes(fix_class) = nclusters;
			ifixflag(nclusters) = 1;
			
			fix_class2 = [fix_class2 fix_class];
			nfix_class = [nfix_class i];
			clusn = [clusn nclusters];
		end
	end
end

% Merge operations
if handles.merge == 1 && ~isempty(nfix_class)
    imerge = 1;% index for the original temperature that will represent all the fixed classes
    bigger_fix = 0;
    for i = 1:length(nfix_class)
        aux = nnz(clustering_results(:,2) == nfix_class(i));
        if aux > bigger_fix
            imerge = i;
            bigger_fix = aux;
        end
    end
    mtemp = clustering_results(imerge,3); % temperature that represents all the fixed classes
    classes(fix_class2) = nfix_class(1); % labels all the fixed classes with the new number
end


if handles.force==0  &&  handles.setclus==0
    forced = USER_DATA{13};
    USER_DATA{14} = forced;
    new_forced = false(size(forced));
    new_forced(fix_class2) = forced(fix_class2);
    clear forced
    USER_DATA{13} = new_forced;
end

% Defines classes
non_clustered = ones(1,size(spikes,1));
cont = 0;
for i = clusn
    class_temp = find(classes == i);
    if ((ifixflag(i)==1) && (~isempty(class_temp)))
        ifixflagc = 1;
    else
        ifixflagc = 0;
    end
    if ((length(class_temp) >= sizemin_clus) || (ifixflagc == 1))
        cont = cont+1;
        eval(['class' num2str(cont) '= class_temp;'])
        non_clustered(class_temp) = 0;
    end
end
nclusters = cont;
rejected = USER_DATA{15};
class0 = find(non_clustered & ~rejected);
clear non_clustered rejected

% Redefines classes
classes = zeros(size(spikes,1),1);
for i = 0:nclusters
    if ~ (isempty(class0) && i==0)
        eval(['classes(class' num2str(i) ') = ' num2str(i) ';']);
    end
end

% Saves new classes
USER_DATA{6} = classes;


% new temperature when merge
if handles.merge == 1 && ~isempty(nfix_class)
    clustering_results(fix_class2,3) = mtemp;
    clustering_results(fix_class2,4) = clustering_results(imerge,4);
end 
clustering_results(:,1) = temp; % GUI temperature
clustering_results(:,5) = minclus; % GUI minimum cluster
 
% Saves new classes and keep fixed classes in 'clustering_results'. 
% Keep the original temperature and cluster number in the fixed spikes.
% The temperature of the non-fixed spikes will be 
% the GUI temperature (temp) and cluster number will be 
% the GUI cluster number (classes)
if (~isempty(fix_class2)) && handles.merge==0 && handles.undo==0 && handles.force==0
    % selects the index of the non-fixed spikes
    % since those are the ones which are going to be updated
    ind_non_fix = 1:length(classes); 
    ind_non_fix(fix_class2) = []; 
    clustering_results(ind_non_fix,4) = classes(ind_non_fix); % classes of the non-fixed spikes in the original clusters column
    clustering_results(ind_non_fix,3) = temp; % temperature of the non-fixed spikes in the original temperature column
end

% update new classes
clustering_results(:,2) = classes;

% If there are no fix and rejected clusters and undo operations, 
% original classes are the same as current classes
% or 0 if they are rejected or manual selected
if isempty(fix_class2) && handles.undo==0 && handles.merge==0 && handles.force==0
    if isfield(handles,'new_manual')
		clustering_results(:,4) = clustering_results(:,4); % same as before
		clustering_results(:,3) = clustering_results(:,3);
		clustering_results(handles.new_manual|classes==0,4) = 0;
		clustering_results(handles.new_manual,3) = temp;
	else
		clustering_results(:,4) = clustering_results(:,2); % clusters
		clustering_results(:,3) = temp; % temperatures
    end
end
clear classes
% Updates clustering_results and clustering_results_bk in USER_DATA
USER_DATA{10} = clustering_results; 
USER_DATA{11} = clustering_results_bk; 
clear clustering_results_bk; 
for i=20:52
    USER_DATA{i} = [];
end

set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);
for i=4:par.max_clus
    eval(['par.fix' num2str(i) '=0;']);
end

set(handles.wave_clus_figure,'userdata',USER_DATA)

% Clear plots
for i=1:4
    cla(eval(['handles.spikes' num2str(i-1)]),'reset');
    cla(eval(['handles.isi' num2str(i-1)]),'reset');
end    
cla(handles.projections); 
hold(handles.projections,'on')
% Plot clusters
ylimit = [];
colors = ['k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];


forced = USER_DATA{13};
figs_num = 6;
opened_figs = cell(1,figs_num);
for i = 0:nclusters
    if ~ (isempty(class0) && i==0)
        %PLOTS SPIKES OR PROJECTIONS
        eval(['max_spikes=min(length(class' num2str(i) '), par.max_spikes_plot);']);
        eval(['sup_spikes=length(class' num2str(i) ');']);
        permut = randperm(sup_spikes);
        permut = permut(1:max_spikes);
        if get(handles.spike_shapes_button,'value') ==1 && get(handles.plot_all_button,'value') ==1
            eval(['line(1:ls,spikes(class' num2str(i) '(permut),:)'',''color'',''' colors(i+1) ''',''Parent'',handles.projections)']);
            xlim(handles.projections, [1 ls])
        elseif get(handles.spike_shapes_button,'value') ==1
            eval(['av   = mean(spikes(class' num2str(i) ',:));']);
            plot(handles.projections,1:ls,av,'color',colors(i+1),'linewidth',2);
            xlim(handles.projections,[1 ls])
        else
            eval(['plot(handles.projections,inspk(class' num2str(i) ',1),inspk(class' num2str(i) ',2),''.' colors(i+1) ''',''markersize'',.5);']);
            axis(handles.projections,'auto');
        end
    
        if i < 4
            clus_ax = eval(['handles.spikes' num2str(i)]); 
            hold(clus_ax,'on')
            eval(['av   = mean(spikes(class' num2str(i) '(:,permut),:));']); % JMG
            eval(['avup = av + par.to_plot_std * std(spikes(class' num2str(i) '(:,permut),:));']); % JMG 
            eval(['avdw = av - par.to_plot_std * std(spikes(class' num2str(i) '(:,permut),:));']); % JMG
            
            if get(handles.plot_all_button,'value') ==1
                eval(['line(1:ls,spikes(class' num2str(i) '(permut),:)'',''color'',''' colors(i+1) ''',''Parent'',clus_ax)']);
                if i==0
                    plot(clus_ax,1:ls,av,'c','linewidth',2)
                    plot(clus_ax,1:ls,avup,'c','linewidth',.5)
                    plot(clus_ax,1:ls,avdw,'c','linewidth',.5)
                else
                    plot(clus_ax,1:ls,av,'k','linewidth',2);
                    plot(clus_ax,1:ls,avup,1:ls,avdw,'color',[.4 .4 .4],'linewidth',.5)
                end
            else
                plot(clus_ax,1:ls,av,'color',colors(i+1),'linewidth',2)
                plot(clus_ax,1:ls,avup,1:ls,avdw,'color',[.65 .65 .65],'linewidth',.5)
            end
            xlim(clus_ax, [1 ls])
            eval(['aux=num2str(length(class' num2str(i) '));']);
            if i>0 
                ylimit = [ylimit; ylim(clus_ax)];
                title(clus_ax,['Cluster ' num2str(i) ':  # ' aux ' (' num2str(nnz(clustering_results(:,2)==i & ~forced(:))) ')'],'Fontweight','bold');
            else            
                title(clus_ax,['Cluster ' num2str(i) ':  # ' aux],'Fontweight','bold');
            end

        else
            par.axes_nr = i+1;
            par.ylimit = ylimit;
            eval(['par.class_to_plot = class' num2str(i) ';']);
            par.plot_all_button = get(handles.plot_all_button,'value');
            USER_DATA{1} = par;
            set(handles.wave_clus_figure,'userdata',USER_DATA)

            if i < 9
                opened_figs{1} = wave_clus_aux('Visible', 'off');
            elseif i < 14
                opened_figs{2} = wave_clus_aux1('Visible', 'off');
            elseif i < 19
                opened_figs{3} = wave_clus_aux2('Visible', 'off');
            elseif i < 24
                opened_figs{4} = wave_clus_aux3('Visible', 'off');
            elseif i < 29
                opened_figs{5} = wave_clus_aux4('Visible', 'off');
            else
                opened_figs{6} = wave_clus_aux5('Visible', 'off');
            %-------------------------------------------------------------------------
            end
        end
    end
end


draw_histograms(handles, 0:min(nclusters,3),USER_DATA);

%Resize axis
if ~isempty(ylimit)
    ymin = min(ylimit(:,1));
    ymax = max(ylimit(:,2));
    for i=1:4
        clus_ax = eval(['handles.spikes' num2str(i-1)]); 
        ylim(clus_ax,[ymin ymax]);
    end
end

for i =1:figs_num
    if ~isempty(opened_figs{i})  
    	set(opened_figs{i},'units','pixel','position',get(0,'screensize'))
	set(opened_figs{i},'Visible', 'on'); 
    end
end

if ~isempty(USER_DATA{5})
    mark_clusters_temperature_diagram(handles,USER_DATA{5},clustering_results)
end
set(handles.file_name,'string', par.file_name_to_show);
axes(handles.projections)
drawnow
