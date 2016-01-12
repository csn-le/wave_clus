function plot_spikes_aux(handles, plot_number)
if ~exist('plot_number','var') || (plot_number== 0)
	USER_DATA = get(handles.wave_clus_aux,'userdata');
	plot_number= 0;
else
	eval(['USER_DATA = get(handles.wave_clus_aux' num2str(plot_number) ',''userdata'');']);
end

par = USER_DATA{1};
spikes = USER_DATA{2};
spk_times = USER_DATA{3};
ls = size(spikes,2);
par.to_plot_std = 1;                % # of std from mean to plot
axes_nr = par.axes_nr;
ylimit = par.ylimit;
class_to_plot = par.class_to_plot;
max_spikes = min(par.max_spikes_plot, length(class_to_plot));
sup_spikes = length(class_to_plot);
forced = USER_DATA{13};
% Plot clusters
colors = ['k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];

sp_axes = eval(['handles.spikes' num2str(axes_nr-1)]); 
cla(sp_axes, 'reset');
hold(sp_axes, 'on');
av   = mean(spikes(class_to_plot,:));
avup = av + par.to_plot_std * std(spikes(class_to_plot,:));
avdw = av - par.to_plot_std * std(spikes(class_to_plot,:));
if par.plot_all_button ==1
    permut = randperm(sup_spikes);
    line(1:ls,spikes(class_to_plot(permut(1:max_spikes)),:)','color',colors(axes_nr),'Parent',sp_axes)
    plot(sp_axes, 1:ls,av,'k','linewidth',2);
    plot(sp_axes, 1:ls,avup,1:ls,avdw,'color',[.4 .4 .4],'linewidth',.5)
else
    plot(sp_axes, 1:ls,av,'color',colors(axes_nr),'linewidth',2)
    plot(sp_axes, 1:ls,avup,1:ls,avdw,'color',[.65 .65 .65],'linewidth',.5)
end
xlim(sp_axes, [1 ls]);
aux = length(class_to_plot);
nforced = nnz(forced(class_to_plot));
title(sp_axes, ['Cluster ' num2str(axes_nr-1) ':  # ' num2str(aux) ' (' num2str(aux-nforced) ')' ],'Fontweight','bold');

%Resize axis
ymin = min(ylimit(:,1));
ymax = max(ylimit(:,2));
ylim(sp_axes, [ymin ymax]);


isi_ax = eval(['handles.isi' num2str(axes_nr-1)]);
times = diff(spk_times(class_to_plot));
% Calculates # ISIs < 3ms  
multi_isi = nnz(times<3); 
% Builds and plots the histogram
eval(['[N,X]=hist(times,0:par.bin_step' num2str(axes_nr-1) ':par.nbins' num2str(axes_nr-1) ');']);
bar(isi_ax, X(1:end-1),N(1:end-1))
eval(['xlim(isi_ax, [0 par.nbins' num2str(axes_nr-1) ']);']);
%eval(['set(get(gca,''children''),''facecolor'',''' colors(axes_nr) ''',''edgecolor'',''' colors(axes_nr) ''',''linewidth'',0.01);']);  %  (FC) why is this commented?
title(isi_ax, [num2str(multi_isi) ' in < 3ms'])
xlabel(isi_ax, 'ISI (ms)');


eval(['set(handles.fix' num2str(4+plot_number*5) '_button,''value'',0);']);
eval(['set(handles.fix' num2str(5+plot_number*5) '_button,''value'',0);']);
eval(['set(handles.fix' num2str(6+plot_number*5) '_button,''value'',0);']);
eval(['set(handles.fix' num2str(7+plot_number*5) '_button,''value'',0);']);
eval(['set(handles.fix' num2str(8+plot_number*5) '_button,''value'',0);']);
