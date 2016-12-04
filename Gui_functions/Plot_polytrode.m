function Plot_polytrode_button(handles)

USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
spikes = USER_DATA{2};
classes = USER_DATA{6};
nclusters = max(classes);
ls = size(spikes,2); % polytrode spike length
lch = par.w_pre + par.w_post; % spike length
nchannels = ls/lch;
filename = par.nick_name;


h_figs=get(0,'children');
h_fig1 = findobj(h_figs,'Name','polytrode');
h_fig2= findobj(h_figs,'Name','polytrode_aux');
h_fig3= findobj(h_figs,'Name','polytrode_aux1');
h_fig4= findobj(h_figs,'name','polytrode_aux2');
close(h_fig1); close(h_fig2); close(h_fig3); close(h_fig4);

% PLOT POLYTRODE CLASSES
colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 
maxc = size(colors,1);

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
            eval(['max_spikes=min(length(class' num2str(i) '),par.max_spikes_plot);']);
            eval(['sup_spikes=length(class' num2str(i) ');']);
            permut=randperm(sup_spikes);

            if i < 6
                h_fig1 = 100;
                figure(h_fig1)
                set(gcf,'numbertitle','off','name',filename)
                if nchannels>4
                    subplot(nchannels,5,i+5*(j-1))
                else
                    subplot(4,5,i+5*(j-1))
                end
                hold on
                if i==nclusters+1
                    eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',[0 0 0])']);
                else
                    eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',[' num2str(colors(mod(i-1,maxc)+1,:)) '])']);
                end
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
                set(gcf,'numbertitle','off','name',[filename '_aux'])
                if nchannels>4
                    subplot(nchannels,5,i-5+5*(j-1))
                else
                    subplot(4,5,i-5+5*(j-1))
                end
                hold on
                if i==nclusters+1
                    eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',[0 0 0])']);
                else
                    eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',[' num2str(colors(mod(i-1,maxc)+1,:)) '])']);
                end
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
                set(gcf,'numbertitle','off','name',[filename '_aux1'])
                if nchannels>4
                    subplot(nchannels,5,i-10+5*(j-1))
                else
                    subplot(4,5,i-10+5*(j-1))
                end
                hold on
                if i==nclusters+1
                    eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',[0 0 0])']);
                else
                    eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',[' num2str(colors(mod(i-1,maxc)+1,:)) '])']);
                end
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
                set(gcf,'numbertitle','off','name',[filename '_aux2'])
                if nchannels>4
                    subplot(nchannels,5,i-15+5*(j-1))
                else
                    subplot(4,5,i-15+5*(j-1))
                end
                hold on
                if i==nclusters+1
                    eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',[0 0 0])']);
                else
                    eval(['plot(spikes(class' num2str(i) '(permut(1:max_spikes)),' num2str(j-1) '*lch + 1 : j*lch)'',''color'',[' num2str(colors(mod(i-1,maxc)+1,:)) '])']);
                end
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
    eval(['print(h_fig1,''-djpeg'',''fig2print_ch_' filename ''')' ]);
end
if ~isempty(h_fig2)
    figure(101); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
    set(gcf,'paperposition',[.25 .25 10.5 7.8])
    eval(['print(h_fig2,''-djpeg'',''fig2print_ch_' filename 'a' ''')' ]);
end
if ~isempty(h_fig3)
    figure(102); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
    set(gcf,'paperposition',[.25 .25 10.5 7.8])
    eval(['print(h_fig3,''-djpeg'',''fig2print_ch_' filename 'b' ''')' ]);
end
if ~isempty(h_fig4)
    figure(103); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
    set(gcf,'paperposition',[.25 .25 10.5 7.8])
    eval(['print(h_fig4,''-djpeg'',''fig2print_ch_' filename 'd' ''')' ]);
end


 
