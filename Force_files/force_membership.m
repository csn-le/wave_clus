function force_membership_new(channel)

handles.par.template_sdnum = 3;             % max radius of cluster in std devs.
handles.par.template_type = 'center';        % nn, center, ml, mahal
clear spikes;
eval(['load times_CSC' num2str(channel) ';']);
classes = cluster_class(:,1);
spk_times = cluster_class(:,2);

f_in  = spikes(find(classes~=0),:);
f_out = spikes(find(classes==0),:);
class_in = classes(find(classes~=0),:);
class_out = force_membership_wc(f_in, class_in, f_out, handles);
classes(find(classes==0)) = class_out;
cluster_class(:,1)=classes;

% Makes Plot
subs_spk = [12 7:11];
subs_isi = [18 13:17];
colors = ['k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];
figure
set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]) 
subplot(3,6,3); title([pwd '   Channel  ' num2str(channel)],'fontsize',14); axis off;
ylimit = [];
for i=1:max(classes)+1
    if ~ (isempty(find(classes==0)) & i==1)
        class_aux = find(classes==i-1);
        subplot(3,6,subs_spk(i))
        hold
        plot(spikes(class_aux,:)','color',colors(i))
        plot(mean(spikes(class_aux,:),1),'k','linewidth',2)
        aux=num2str(length(class_aux));
        eval(['title([''Cluster ' num2str(i-1) ':  # ' aux '''],''Fontweight'',''bold'')']);
        xlim([1 size(spikes,2)]);
        ylimit = [ylimit;ylim];
        
        subplot(3,6,subs_isi(i))
        times=diff(spk_times(class_aux));
        [N,X]=hist(times,0:1:100);
        bar(X(1:end-1),N(1:end-1))
        xlim([0 100]);
        set(get(gca,'children'),'facecolor',colors(i),'linewidth',0.01); 
        title([num2str(sum(N(1:3))) ' in < 3ms'])
        xlabel('ISI (ms)');
    end
end

% Rescale spike's axis 
if ~isempty(ylimit)
    ymin = min(ylimit(:,1));
    ymax = max(ylimit(:,2));
end
for i=1:max(classes)+1
    if ~ (isempty(find(classes==0)) & i==1)
        subplot(3,6,subs_spk(i))
        ylim([ymin ymax]);
    end
end

outfile=['times_CSCf' num2str(channel)];
save(outfile, 'cluster_class')
eval(['print -djpeg fig2print_ch' num2str(channel) 'f;']);

