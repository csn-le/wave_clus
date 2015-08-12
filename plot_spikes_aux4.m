function Plot_Spikes_aux(handles)
USER_DATA = get(handles.wave_clus_aux4,'userdata');
par = USER_DATA{1};
spikes = USER_DATA{2};
spk_times = USER_DATA{3};
inspk = USER_DATA{7};
ls = size(spikes,2);
par.to_plot_std = 1;                % # of std from mean to plot
axes_nr = par.axes_nr;
ylimit = par.ylimit;
class_to_plot = par.class_to_plot;
max_spikes = min(par.max_spikes, length(class_to_plot));
sup_spikes = length(class_to_plot);

% Plot clusters
colors = ['k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];
eval(['axes(handles.spikes' num2str(axes_nr-1) ');']); 
cla reset
hold on
av   = mean(spikes(class_to_plot,:));
avup = av + par.to_plot_std * std(spikes(class_to_plot,:));
avdw = av - par.to_plot_std * std(spikes(class_to_plot,:));
if par.plot_all_button ==1
    permut=randperm(sup_spikes);
    plot(spikes(class_to_plot(permut(1:max_spikes)),:)','color',colors(axes_nr));
    plot(1:ls,av,'k','linewidth',2);
    plot(1:ls,avup,1:ls,avdw,'color',[.4 .4 .4],'linewidth',.5)
else
    plot(1:ls,av,'color',colors(axes_nr),'linewidth',2)
    plot(1:ls,avup,1:ls,avdw,'color',[.65 .65 .65],'linewidth',.5)
end
xlim([1 ls])
aux = num2str(length(class_to_plot));
eval(['title([''Cluster ' num2str(axes_nr-1) ':  # ' aux '''],''Fontweight'',''bold'')']);
eval(['axes(handles.isi' num2str(axes_nr-1) ');']); 
times = diff(spk_times(class_to_plot));
% Calculates # ISIs < 3ms  
bin_step_temp = 1;
eval(['[N,X]=hist(times,0:bin_step_temp:par.nbins' num2str(axes_nr-1) ');']);
multi_isi= sum(N(1:3)); 
% Builds and plots the histogram
eval(['[N,X]=hist(times,0:par.bin_step' num2str(axes_nr-1) ':par.nbins' num2str(axes_nr-1) ');']);
bar(X(1:end-1),N(1:end-1))
eval(['xlim([0 par.nbins' num2str(axes_nr-1) ']);']);
eval(['set(get(gca,''children''),''facecolor'',''' colors(axes_nr) ''',''edgecolor'',''' colors(axes_nr) ''',''linewidth'',0.01);']);    
title([num2str(multi_isi) ' in < 3ms'])
%title([num2str(sum(N(1:3))) ' in < 3ms'])
xlabel('ISI (ms)');

%Resize axis
ymin = min(ylimit(:,1));
ymax = max(ylimit(:,2));
eval(['axes(handles.spikes' num2str(axes_nr-1) '); ylim([ymin ymax])'])

set(handles.fix24_button,'value',0);
set(handles.fix25_button,'value',0);
set(handles.fix26_button,'value',0);
set(handles.fix27_button,'value',0);
set(handles.fix28_button,'value',0);
