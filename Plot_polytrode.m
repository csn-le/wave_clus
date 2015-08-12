function Plot_polytrode_button(handles)

USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
spikes = USER_DATA{2};
classes = USER_DATA{6};
nclusters = max(classes);
par.to_plot_std = 1;
ls = size(spikes,2); % polytrode spike length
lch = par.w_pre + par.w_post; % spike length
nchannels = ls/lch;
filename = par.filename;


h_figs=get(0,'children');
h_fig1 = findobj(h_figs,'Name','polytrode');
h_fig2= findobj(h_figs,'Name','polytrode_aux');
h_fig3= findobj(h_figs,'Name','polytrode_aux1');
h_fig4= findobj(h_figs,'name','polytrode_aux2');
close(h_fig1); close(h_fig2); close(h_fig3); close(h_fig4);

% PLOT POLYTRODE CLASSES
colors = ['b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];

for i = 1:nclusters
    eval(['class' num2str(i) '= find( classes==' num2str(i) ');']);
end
z = nclusters + 1;
eval(['class' num2str(z) '= find( classes==' num2str(0) ');']);
 
for j=1:nchannels
    ylimit = [];
    h_figs = [];
    for i=1:nclusters+1
        eval(['nspikes = length(class' num2str(i) ');'])
        if nspikes
            eval(['av   = mean(spikes(class' num2str(i) ',' num2str(j-1) '*lch + 1 : j*lch ));']);
            eval(['avup = av + par.to_plot_std * std(spikes(class' num2str(i) ',' num2str(j-1) '*lch + 1 : j*lch));']);
            eval(['avdw = av - par.to_plot_std * std(spikes(class' num2str(i) ',' num2str(j-1) '*lch + 1 : j*lch));']);
            eval(['max_spikes=min(length(class' num2str(i) '),par.max_spikes);']);
            eval(['sup_spikes=length(class' num2str(i) ');']);
            permut=randperm(sup_spikes);

            if i < 6
                h_fig1 = 100;
                figure(h_fig1)
                t = strcat(char(filename(1:9)) );
                set(gcf,'numbertitle','off','name',t)
                if nchannels>4
                    subplot(nchannels,5,i+5*(j-1))
                else
                    subplot(4,5,i+5*(j-1))
                end
                hold on
                if i==nclusters+1,colors(i)='k'; end;
                eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',''' colors(i) ''')']);
                plot(1:lch,av,'k','linewidth',2)
                plot(1:lch,avup,1:lch,avdw,'color' , [.4 .4 .4], 'linewidth' ,0.5)
    %             xlim([1 lch])
                axis tight
                if i<nclusters+1, ylimit = [ylimit;get(gca,'YLim')]; end
                h_figs = [h_figs,gca];
                if i==1
                    eval(['ylabel([''Channel '  num2str(j) '''],''Fontweight'',''bold'')']);
                end
                if j==1
                    if i==nclusters+1, k=0; else k=i; end
                    eval(['aux=num2str(length(class' num2str(i) '));']);
                    eval(['title([''Cluster ' num2str(k) ':  # ' aux '''],''Fontweight'',''bold'')']);
                end

            elseif i < 11
                h_fig2 = 101;
                figure(h_fig2)
                t = strcat(char(filename(1:9)), '_aux' );
                set(gcf,'numbertitle','off','name',t)
                if nchannels>4
                    subplot(nchannels,5,i-5+5*(j-1))
                else
                    subplot(4,5,i-5+5*(j-1))
                end
                hold on
                if i==nclusters+1,colors(i)='k'; end
                eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',''' colors(i) ''')']);
                plot(1:lch,av,'k','linewidth',2)
                plot(1:lch,avup,1:lch,avdw,'color' , [.4 .4 .4], 'linewidth' ,0.5)
                axis tight
                if i<nclusters+1, ylimit = [ylimit;get(gca,'YLim')]; end
                h_figs = [h_figs,gca];
                if i==6
                    eval(['ylabel([''Channel '  num2str(j) '''],''Fontweight'',''bold'')']);
                end
                if j==1
                    if i==nclusters+1, k=0; else k=i; end
                    eval(['aux=num2str(length(class' num2str(i) '));']);
                    eval(['title([''Cluster ' num2str(k) ':  # ' aux '''],''Fontweight'',''bold'')']);
                end

            elseif i < 16
                h_fig3 = 102;
                figure(h_fig3)
                t = strcat(char(filename(1:9)), '_aux1' );
                set(gcf,'numbertitle','off','name',t)
                if nchannels>4
                    subplot(nchannels,5,i-10+5*(j-1))
                else
                    subplot(4,5,i-10+5*(j-1))
                end
                hold on
                if i==nclusters+1,colors(i)='k'; end
                eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',''' colors(i) ''')']);
                plot(1:lch,av,'k','linewidth',2)
                plot(1:lch,avup,1:lch,avdw,'color' , [.4 .4 .4], 'linewidth' ,0.5)
                axis tight
                if i<nclusters+1, ylimit = [ylimit;get(gca,'YLim')]; end
                h_figs = [h_figs,gca];
                if i==11
                    eval(['ylabel([''Channel '  num2str(j) '''],''Fontweight'',''bold'')']);
                end
                if j==1
                    if i==nclusters+1, k=0; else k=i; end
                    eval(['aux=num2str(length(class' num2str(i) '));']);
                    eval(['title([''Cluster ' num2str(k) ':  # ' aux '''],''Fontweight'',''bold'')']);
                end

            elseif i < 21
                h_fig4 = 103;
                figure(h_fig4)
                t = strcat(char(filename(1:9)), '_aux2' );
                set(gcf,'numbertitle','off','name',t)
                if nchannels>4
                    subplot(nchannels,5,i-15+5*(j-1))
                else
                    subplot(4,5,i-15+5*(j-1))
                end
                hold on
                if i==nclusters+1,colors(i)='k'; end
                eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',''' colors(i) ''')']);
                plot(1:lch,av,'k','linewidth',2)
                plot(1:lch,avup,1:lch,avdw,'color' , [.4 .4 .4], 'linewidth' ,0.5)
                axis tight
                if i<nclusters+1, ylimit = [ylimit;get(gca,'YLim')]; end
                h_figs = [h_figs,gca];
                if i==11
                    eval(['ylabel([''Channel '  num2str(j) '''],''Fontweight'',''bold'')']);
                end
                if j==1
                    if i==nclusters+1, k=0; else k=i; end
                    eval(['aux=num2str(length(class' num2str(i) '));']);
                    eval(['title([''Cluster ' num2str(k) ':  # ' aux '''],''Fontweight'',''bold'')']);
                end
            end
        end;
        
  end
  ymin = min(ylimit(:,1)); ymax = max(ylimit(:,2)); 
  for k=1:length(h_figs), set(h_figs(k),'Ylim',[ymin ymax]),end
end

% SAVE FIGURE
if ~isempty(h_fig1)
    figure(100); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
    set(gcf,'paperposition',[.25 .25 10.5 7.8])
    eval(['print(h_fig1,''-djpeg'',''fig2print_ch_' filename(1:end-4) ''')' ]);
end
if ~isempty(h_fig2)
    figure(101); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
    set(gcf,'paperposition',[.25 .25 10.5 7.8])
    eval(['print(h_fig2,''-djpeg'',''fig2print_ch_' filename(1:end-4) 'a' ''')' ]);
end
if ~isempty(h_fig3)
    figure(102); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
    set(gcf,'paperposition',[.25 .25 10.5 7.8])
    eval(['print(h_fig3,''-djpeg'',''fig2print_ch_' filename(1:end-4) 'b' ''')' ]);
end
if ~isempty(h_fig4)
    figure(103); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
    set(gcf,'paperposition',[.25 .25 10.5 7.8])
    eval(['print(h_fig4,''-djpeg'',''fig2print_ch_' filename(1:end-4) 'd' ''')' ]);
end


 
