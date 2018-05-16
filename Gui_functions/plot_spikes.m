function plot_spikes(handles)
set(handles.file_name,'string','Plotting...'); 
drawnow;
if exist('groot','builtin')
    if isprop(handles.wave_clus_figure,'GraphicsSmoothing')
        set(handles.wave_clus_figure,'GraphicsSmoothing','off');
    end
    try
        set(groot,'defaultfiguregraphicssmoothing','off');
        set(groot,'DefaultAxesFontSize',8)
    end
end

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
h_figs = get(0,'children');
close(findobj(h_figs,'flat','tag','wave_clus_aux'));
for w = 1:5
    try
        close(findobj(h_figs,'flat','tag',['wave_clus_aux' num2str(w)]));
    catch
        h_figs = get(0,'children'); %to fix Invalid handle in old Matlabs
        close(findobj(h_figs,'flat','tag',['wave_clus_aux' num2str(w)]));
    end
end
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
USER_DATA{11} = clustering_results; 

% Forcing
if handles.force==1
    for i = classes_names
        ind = find(clustering_results(:,2)==i); % get index of GUI class
        oclass = clustering_results(ind(1),4); % get original class
        otemp = clustering_results(ind(1),3); % get original temperature
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
    sizemin_clus = minclus;
elseif handles.setclus == 1
	sizemin_clus = 1;
else
    sizemin_clus = minclus;
end

clusn = find(cluster_sizes >= sizemin_clus);
nclusters = length(clusn);

% Get fixed clusters
fix_class2 = [];
nfix_class = [];

if ~isfield(handles,'new_manual')
	for fi = 1:3
        if get(eval(['handles.fix' num2str(fi) '_button']),'value') ==1
            nclusters = nclusters +1;
            fix_class = USER_DATA{19+fi}';
            classes(classes==nclusters)=0;
            classes(fix_class)=nclusters;
            ifixflag(nclusters)=1;
            fix_class2 = [fix_class2 fix_class]; 
            nfix_class = [nfix_class fi];
            clusn = [clusn nclusters];    
        end
	end
	% Get fixed clusters from aux figures
	for i=4:min(par.max_clus,33)
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
    imerge = find(clustering_results(:,2) ==nfix_class(imerge),1);
    mtemp = clustering_results(imerge,3); % temperature that represents all the fixed classes
    classes(fix_class2) = nfix_class(1); % labels all the fixed classes with the new number
end

if ~isfield(handles,'unforce')
    if handles.force==0  &&  handles.setclus==0
        forced = USER_DATA{13};
        USER_DATA{14} = forced;
        new_forced = false(size(forced));
        new_forced(fix_class2) = forced(fix_class2);
        clear forced
        USER_DATA{13} = new_forced;
    elseif handles.force==0
        USER_DATA{14} = USER_DATA{13};
    end
end
% Defines classes
non_clustered = ones(1,size(spikes,1));
nclusters = 0;
for i = clusn
    class_temp = find(classes == i);
    if ((ifixflag(i)==1) && (~isempty(class_temp)))
        ifixflagc = 1;
    else
        ifixflagc = 0;
    end
    if ((length(class_temp) >= sizemin_clus) || (ifixflagc == 1))
        nclusters = nclusters+1;
        eval(['class' num2str(nclusters) '= class_temp;'])
        non_clustered(class_temp) = 0;
    end
end
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
    if isfield(handles,'new_spc_classes')
            clustering_results(ind_non_fix,4) = handles.new_spc_classes(ind_non_fix);
    end
    if handles.setclus == 0
        clustering_results(ind_non_fix,3) = temp; % temperature of the non-fixed spikes in the original temperature column
    end
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
	elseif handles.setclus == 0
        if isfield(handles,'new_spc_classes')
            clustering_results(:,4) = handles.new_spc_classes;
        else
            clustering_results(:,4) = clustering_results(:,2); % clusters
        end
		clustering_results(:,3) = temp; % temperatures
    end
end

clear classes
% Updates clustering_results in USER_DATA
USER_DATA{10} = clustering_results; 

for i=20:52
    USER_DATA{i} = [];
end

for i=4:par.max_clus
    eval(['par.fix' num2str(i) '=0;']);
end

set(handles.wave_clus_figure,'userdata',USER_DATA)
ax_v=[];
% Clear plots
for i=1:3
    set(eval(['handles.fix' num2str(i) '_button']),'value',0);
    cla(eval(['handles.spikes' num2str(i)]),'reset');
	ax_v(end+1)=eval(['handles.spikes' num2str(i)]);
    cla(eval(['handles.isi' num2str(i)]),'reset');
end    
% cla(handles.isi0);
cla(handles.spikes0,'reset');
cla(handles.projections);
hold(handles.projections,'on')

% Plot clusters
ylimit = [];
colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 
maxc = size(colors,1);
forced = USER_DATA{13};
figs_num = 6;
opened_figs = cell(1,figs_num);
for i = 0:nclusters
    tmpy = []; %as a flag to don't make the same vector twice
    if ~ (isempty(class0) && i==0)
        %PLOTS SPIKES OR PROJECTIONS
        class_i = eval(['class' num2str(i)]);
        sup_spikes = length(class_i);
        max_spikes = min(sup_spikes, par.max_spikes_plot);
        permut = randperm(sup_spikes,max_spikes);
        xlim(handles.projections,'manual');
        if get(handles.spike_shapes_button,'value') ==1 && (get(handles.plot_all_button,'value') ==1) && ~strcmp(par.all_classes_ax,'mean')
            % optimizing for speed:
            tmpy=spikes(class_i(permut),:);
            tmpn=size(tmpy,1);
            tmpx=repmat([1:ls NaN]',1,tmpn);
            tmpx=reshape(tmpx,numel(tmpx),1);
            tmpy=[tmpy'; repmat(NaN,1,tmpn)];
            tmpy=reshape(tmpy,numel(tmpy),1);
            line(tmpx,tmpy,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'Parent',handles.projections,'Visible','off');
			xlim(handles.projections, [1 ls])
        elseif get(handles.spike_shapes_button,'value') ==1
            av   = mean(spikes(class_i,:));
            plot(handles.projections,1:ls,av,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'linewidth',2);
            xlim(handles.projections,[1 ls])
        else
            plot(handles.projections,inspk(class_i,1),inspk(class_i,2),'.','Color',colors(mod(i-1,maxc)+1,:)*(i~=0),'markersize',.5);
            axis(handles.projections,'auto');
        end
    
        if i < 4
            clus_ax = eval(['handles.spikes' num2str(i)]); 
            xlim(clus_ax,'manual');
            xlim(clus_ax,[1 ls]);
            hold(clus_ax,'on')
            
            av   = mean(spikes(class_i,:));
            avup = av + par.to_plot_std * std(spikes(class_i,:));
            avdw = av - par.to_plot_std * std(spikes(class_i,:));
                      
            if get(handles.plot_all_button,'value') ==1
                % optimizing for speed:
                if isempty(tmpy)
                    tmpy=spikes(class_i(permut),:);
                    tmpn=size(tmpy,1);
                    tmpx=repmat([1:ls NaN]',1,tmpn);
                    tmpx=reshape(tmpx,numel(tmpx),1);
                    tmpy=[tmpy'; repmat(NaN,1,tmpn)];
                    tmpy=reshape(tmpy,numel(tmpy),1);
                end
                line(tmpx,tmpy,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'Parent',clus_ax);
				
				if i==0
                    line(1:ls,av,'color','c','linewidth',2,'Parent',clus_ax)
                    line(1:ls,avup,'color','c','linewidth',0.5,'Parent',clus_ax)
                    line(1:ls,avdw,'color','c','linewidth',0.5,'Parent',clus_ax)
                else
                    line(1:ls,av,'color','k','linewidth',2,'Parent',clus_ax)
                    line(1:ls,avdw,'color',[.4 .4 .4],'linewidth',0.5,'Parent',clus_ax)
                    line(1:ls,avup,'color',[.4 .4 .4],'linewidth',0.5,'Parent',clus_ax)
                end
            else
                plot(clus_ax,1:ls,av,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'linewidth',2)
                plot(clus_ax,1:ls,avup,1:ls,avdw,'color',[.65 .65 .65],'linewidth',.5)
            end
            
            eval(['aux=num2str(length(class' num2str(i) '));']);
            if i>0 
                ylim(clus_ax,'auto');
                ylimit = [ylimit; ylim(clus_ax)];
                title(clus_ax,['Cluster ' num2str(i) ':  # ' aux ' (' num2str(nnz(clustering_results(:,2)==i & ~forced(:))) ')'],'Fontweight','bold');
            else            
                title(clus_ax,['Cluster ' num2str(i) ':  # ' aux],'Fontweight','bold');
                xlim(clus_ax, [1 ls])
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
            elseif i < 34
                opened_figs{6} = wave_clus_aux5('Visible', 'off');
            %-------------------------------------------------------------------------
            end
        end
    end
end


draw_histograms(handles, 0:min(nclusters,3),USER_DATA);

if ~isempty(USER_DATA{5})
    mark_clusters_temperature_diagram(handles,USER_DATA{5},clustering_results)
end
set(handles.file_name,'string', par.file_name_to_show);

set(allchild(handles.projections),'Visible','on')

%Resize axis
if ~isempty(ylimit)
    ymin = min(ylimit(:,1));
    ymax = max(ylimit(:,2));
    ylim(handles.spikes0,[ymin ymax]);
    if get(handles.spike_shapes_button,'value') ==1
        ylim(handles.projections,[ymin ymax]);
    end
    linkaxes(ax_v,'xy'); %drawnow inside
    ylim(ax_v(1),[ymin ymax]);
else
    drawnow
end

for i =1:figs_num
    if ~isempty(opened_figs{i})  
        set(opened_figs{i},'units','normalized','outerposition',[0 0 1 1])
        set(opened_figs{i},'Visible', 'on'); 
    end
end

if exist('groot','builtin')
    set(groot,'defaultfiguregraphicssmoothing','remove')
    set(groot,'DefaultAxesFontSize','remove')
end
