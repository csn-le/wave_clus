function plot_spikes_aux(handles, USER_DATA,plot_number)

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
colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 
maxc = size(colors,1);


sp_axes = eval(['handles.spikes' num2str(axes_nr-1)]); 
cla(sp_axes, 'reset');
hold(sp_axes, 'on');
av   = mean(spikes(class_to_plot,:));
avup = av + par.to_plot_std * std(spikes(class_to_plot,:));
avdw = av - par.to_plot_std * std(spikes(class_to_plot,:));
ylim(sp_axes, 'manual');
xlim(sp_axes, 'manual');
if par.plot_all_button ==1
    permut = randperm(sup_spikes);
    
    tmpy=spikes(class_to_plot(permut(1:max_spikes)),:);
    tmpn=size(tmpy,1);
    tmpx=repmat([1:ls NaN]',1,tmpn);
    tmpx=reshape(tmpx,numel(tmpx),1);
    tmpy=[tmpy'; repmat(NaN,1,tmpn)];
    tmpy=reshape(tmpy,numel(tmpy),1);
    line(tmpx,tmpy,'color',colors(mod(axes_nr-2,maxc)+1,:),'Parent',sp_axes);    

    plot(sp_axes, 1:ls,av,'k','linewidth',2);
    plot(sp_axes, 1:ls,avup,1:ls,avdw,'color',[.4 .4 .4],'linewidth',.5)
else
    plot(sp_axes, 1:ls,av,'color',colors(mod(axes_nr-2,maxc)+1,:),'linewidth',2)
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
try
    xlim(isi_ax,'manual');
    eval(['[N,X]=hist(times,0:par.bin_step' num2str(axes_nr-1) ':par.nbins' num2str(axes_nr-1) ');']);
    bar(isi_ax, X(1:end-1),N(1:end-1))
    eval(['xlim(isi_ax, [0 par.nbins' num2str(axes_nr-1) ']);']);
    %eval(['set(get(gca,''children''),''facecolor'',''' colors(axes_nr) ''',''edgecolor'',''' colors(axes_nr) ''',''linewidth'',0.01);']);  %  (FC) why is this commented?
    title(isi_ax, [num2str(multi_isi) ' in < 3ms'])
    xlabel(isi_ax, 'ISI (ms)');
catch
    warning(['Error in the ISI plot of the Cluster ' num2str(axes_nr-1)] )
end



eval(['set(handles.fix' num2str(4+plot_number*5) '_button,''value'',0);']);
eval(['set(handles.fix' num2str(5+plot_number*5) '_button,''value'',0);']);
eval(['set(handles.fix' num2str(6+plot_number*5) '_button,''value'',0);']);
eval(['set(handles.fix' num2str(7+plot_number*5) '_button,''value'',0);']);
eval(['set(handles.fix' num2str(8+plot_number*5) '_button,''value'',0);']);
