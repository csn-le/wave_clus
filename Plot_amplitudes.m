function Plot_polytrode_amplitudes(handles)

USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
spikes = USER_DATA{2};

colors = ['b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];

h_fig = 100;
figure(h_fig)
nclasses = max(classes);
nchannels = par.channels;
inputs = nchannels;
ls = size(spikes,2); % pseudospike length
lch = ls/nchannels; % spike channel length
amps = zeros(ls,nchannels);
filename = par.nick_name;

% AMPLITUDES
for i=1:nchannels
    aux = spikes(1:length(classes) ,(i-1)*lch + 1 : i*lch );
    amps(1:length(classes),i) = max(aux,[],2);
end  
         
% PLOTS         
for i=1:inputs
    for j=i+1:inputs
        subplot(inputs,inputs,(i-1)*inputs+j)
        hold on
        for k=1:nclasses
            class_aux = find(classes==k);
            max_spikes = min(par.max_spikes_plot,length(class_aux));
            plot(amps(class_aux(1:max_spikes),i),amps(class_aux(1:max_spikes),j),['.' colors(k)],'markersize',.5)
            axis off
        end
    end
end

set(gcf,'numbertitle','off','name',filename,'menubar','none')

% SAVE FIGURE
set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
set(gcf,'paperposition',[.25 .25 10.5 7.8])
eval(['print(h_fig,''-djpeg'',''fig2print_amp_' filename ''')' ]);


