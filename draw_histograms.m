function draw_histograms(handles, c2plot)

USER_DATA = get(handles.wave_clus_figure,'userdata');
spk_times = USER_DATA{3};
classes = USER_DATA{6};
par = USER_DATA{1};

for i = c2plot 
    eval(['axes(handles.isi' num2str(i) ');']); 
    times=diff(spk_times(classes==i));
    % Calculates # ISIs < 3ms  
    multi_isi = nnz(times < 3); 
    % Builds and plots the histogram
    eval(['[N,X]=hist(times,0:par.bin_step' num2str(i) ':par.nbins' num2str(i) ');']);
    bar(X(1:end-1),N(1:end-1))
    eval(['xlim([0 par.nbins' num2str(i) ']);']);
    title([num2str(multi_isi) ' in < 3ms'])
    xlabel('ISI (ms)');
end