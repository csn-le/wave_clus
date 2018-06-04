function wc_polytrode_figs(input,show_figs)
%Plot 
%if input is a string, load times from that file
%otherwise is a cell(USER_DATA) from wave_clus_gui

if ~exist('show_figs','var')
	show_figs = false;
end

if ~exist('input','var') || isempty(input)
    h_figs = get(0,'children');
    h_fig = findobj(h_figs,'tag','wave_clus_figure');
    
    if ~isempty(h_fig)
        input = get(h_fig,'userdata');
    else
        aux = dir('times_*.mat');
        if length(aux)==1
            input = aux.name;
        else
            error('empty input, without opened GUI and more that one times in pwd')
        end
    end
end
    
if isstr(input)
    load(input,'spikes','par','cluster_class')
    current_par = set_parameters();
    par = update_parameters(par,current_par,'batch_plot');
	classes = cluster_class(:,1); 
	[~,aux,~] = fileparts(input);
	filename = aux(7:end);
	par.to_plot_std = 1;
else
	par = input{1};
	spikes = input{2};
	classes = input{6};
	filename = par.nick_name;
end

if exist('groot','builtin')
    set(groot,'defaultfiguregraphicssmoothing','off');
end

nclusters = max(classes);
ls = size(spikes,2); % polytrode spike length

if par.channels == 1 || isnan(par.channels)
    lch = par.w_pre + par.w_post; % spike length
    nchannels = ls/lch;
else
    nchannels = par.channels;
    lch = ls/nchannels; % spike length
end
h_figs=get(0,'children');
h_fig1 = findobj(h_figs,'Name',['^' filename 'pol_plot_']);
close(h_fig1);

% PLOT POLYTRODE CLASSES
colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 
maxc = size(colors,1);

clas_ix = cell(nclusters+1,1);
for i = 1:nclusters
	clas_ix{i} = find(classes==i);
end
clas_ix{nclusters+1} = find(classes==0);

nfigs = ceil((nclusters+1)/5);
pp_figs = [];
for i = 1: nfigs
	pp_figs(i) = figure('numbertitle','off','name',[filename '_aux'],'Visible','off');
end

for j=1:nchannels
    ylimit = [];
    ch_ax = [];
    for i=1:nclusters+1
		nspikes = length(clas_ix{i});
        if nspikes
            av   = mean(spikes(clas_ix{i},(j-1)*lch + 1 : j*lch ));
            avup = av + par.to_plot_std * std(spikes(clas_ix{i},(j-1)*lch + 1 : j*lch));
			avdw = av - par.to_plot_std * std(spikes(clas_ix{i},(j-1)*lch + 1 : j*lch));
            max_spikes = min(nspikes,par.max_spikes_plot);
            permut = randperm(nspikes,max_spikes);

			
			nfig = ceil(i/5);
            set(0, 'CurrentFigure', pp_figs(nfig))

            if nchannels>4
                col_correction = (nfig-1)*5;
                subplot(nchannels,5,i+5*(j-1)-col_correction)
            else
                col_correction = (nfig-1)*5;
                subplot(4,5,i+5*(j-1)-col_correction)
            end
			hold on
            
            tmpy = spikes(clas_ix{i}(permut),(j-1) *lch + 1 : j*lch);
            tmpn = size(tmpy,1);
            tmpx = repmat([1:lch NaN]',1,tmpn);
            tmpx = reshape(tmpx,numel(tmpx),1);
            tmpy = [tmpy'; repmat(NaN,1,tmpn)];
            tmpy = reshape(tmpy,numel(tmpy),1);
                        
			if i==nclusters+1
                line(tmpx,tmpy,'color',[0 0 0]);
			else
				line(tmpx,tmpy,'color',colors(mod(i-1,maxc)+1,:));
			end
			plot(1:lch,av,'k','linewidth',2)
			plot(1:lch,avup,1:lch,avdw,'color' , [.4 .4 .4], 'linewidth' ,0.5)
			axis tight
			if i<nclusters+1, ylimit = [ylimit;get(gca,'YLim')]; end
			ch_ax = [ch_ax,gca];
			if i==1
				ylabel(['Channel '  num2str(j)],'Fontweight','bold');
			end
			if j==1
				if i==nclusters+1, k=0; else k=i; end
				title(['Cluster ' num2str(k) ':  # ' num2str(nspikes)],'Fontweight','bold');
			end
			
           
        end
        
  end
  ymin = min(ylimit(:,1)); ymax = max(ylimit(:,2)); 
  for k=1:length(ch_ax)
	set(ch_ax(k),'Ylim',[ymin ymax])
  end
end

 fig_name = @(x) ['Pol_plot_' filename '_fig' num2str(x) '.png'];

% SAVE FIGURE
 for h = 1:length(pp_figs)
	set( pp_figs(h),'papertype','usletter','paperorientation','portrait','paperunits','inches')
    set( pp_figs(h),'paperposition',[.25 .25 10.5 7.8])
    print( pp_figs(h),'-dpng',fig_name(h));
	 if ~show_figs
		close( pp_figs(h))
     else
         set( pp_figs(h),'Visible','on')
         set(pp_figs(h) ,'CloseRequestFcn',@(x,y) cellfun(@delete,num2cell(pp_figs)))
	 end
 end
 h = length(pp_figs)+1;
 while exist(fig_name(h), 'file')==2 
     delete(fig_name(h))
     h = h+1;
 end
 if exist('groot','builtin')
    set(groot,'defaultfiguregraphicssmoothing','remove')
end
