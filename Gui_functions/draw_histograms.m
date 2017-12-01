function draw_histograms(handles, c2plot,USER_DATA)
%This function plot the histogram of the cluster number 'c2plot'
%handles should have the axes of that histogram. USER_DATA is the USER_DATA
%of wave_clus.


spk_times = USER_DATA{3};
classes = USER_DATA{6};
par = USER_DATA{1};

for i = c2plot 
    if classes == 0
        rejected = USER_DATA{15};
        times = diff(spk_times(classes(:)==i & ~rejected(:)));
        clear rejected 
    else
        times = diff(spk_times(classes==i));
    end
    % Calculates # ISIs < 3ms  
    multi_isi = nnz(times < 3); 
    % Builds and plots the histogram
    isi_ax = eval(['handles.isi' num2str(i)]);
    xlim(isi_ax,'manual');
    eval(['[N,X]=hist(times,0:par.bin_step' num2str(i) ':par.nbins' num2str(i) ');']);
    bar(isi_ax,X(1:end-1),N(1:end-1))
    eval(['xlim(isi_ax,[0 par.nbins' num2str(i) ']);']);
    title(isi_ax,[num2str(multi_isi) ' in < 3ms'])
    xlabel(isi_ax,'ISI (ms)');
end
