% PROGRAM Do_clustering.
% Does clustering on all files in Files.txt
% Runs after Get_spikes.

function Do_clustering_sim(simulations)

print2file =1;                              %for saving printouts.
%print2file =0;                              %for printing printouts.

handles.par.w_pre=20;                       %number of pre-event data points stored
handles.par.w_post=44;                      %number of post-event data points stored
handles.par.detection = 'pos';              %type of threshold
handles.par.stdmin = 5.00;                  %minimum threshold
handles.par.stdmax = 50;                    %maximum threshold
handles.par.interpolation = 'y';            %interpolation for alignment
handles.par.int_factor = 2;                 %interpolation factor
handles.par.detect_fmin = 300;              %high pass filter for detection (default 300)
handles.par.detect_fmax = 3000;             %low pass filter for detection (default 3000)
handles.par.sort_fmin = 300;                %high pass filter for sorting (default 300)
handles.par.sort_fmax = 3000;               %low pass filter for sorting (default 3000)

handles.par.max_spk = 20000;                % max. # of spikes before starting templ. match.
handles.par.template_type = 'center';       % nn, center, ml, mahal
handles.par.template_sdnum = 3;             % max radius of cluster in std devs.
handles.par.permut = 'y';                   % for selection of random 'par.max_spk' spikes before starting templ. match. 
% handles.par.permut = 'n';                 % for selection of the first 'par.max_spk' spikes before starting templ. match.

handles.par.features = 'wav';               %choice of spike features
handles.par.inputs = 10;                    %number of inputs to the clustering
handles.par.scales = 4;                     %scales for wavelet decomposition
if strcmp(handles.par.features,'pca');      %number of inputs to the clustering for pca
    handles.par.inputs=3; 
end

handles.par.mintemp = 0;                    %minimum temperature
handles.par.maxtemp = 0.201;                %maximum temperature
handles.par.tempstep = 0.01;                %temperature step
handles.par.num_temp = floor(...
(handles.par.maxtemp - ...
handles.par.mintemp)/handles.par.tempstep); %total number of temperatures 
handles.par.stab = 0.8;                     %stability condition for selecting the temperature
handles.par.SWCycles = 100;                 %number of montecarlo iterations
handles.par.KNearNeighb = 11;               %number of nearest neighbors
handles.par.randomseed = 0;                 % if 0, random seed is taken as the clock value
%handles.par.randomseed = 147;              % If not 0, random seed   
handles.par.fname_in = 'tmp_data';          % temporary filename used as input for SPC

handles.par.min_clus_abs = 20;              %minimum cluster size (absolute value)
handles.par.min_clus_rel = 0.005;           %minimum cluster size (relative to the total nr. of spikes)
%handles.par.temp_plot = 'lin';               %temperature plot in linear scale
handles.par.temp_plot = 'log';              %temperature plot in log scale
handles.par.force_auto = 'y';               %automatically force membership if temp>3.
handles.par.max_spikes = 5000;              %maximum number of spikes to plot.

handles.par.sr = 24000;                     %sampling frequency, in Hz.

figure
set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]) 


for k=1:length(simulations)
    tic
    file_to_cluster = ['C_sim_script_' num2str(simulations(k))];

    %Load continuous data (for ploting)
    % for spikes only data this is not possible
    % therefore we use a try-catch case for this
    continuous_data_av = 1;
    flag0 = 0;
    try 
        eval(['load ' char(file_to_cluster) '.mat;']);
        if length(data)>60*handles.par.sr
            x = data(1:60*handles.par.sr)';       %will plot just 60 sec.
        else
            x = data;
        end
        clear data
        eval(['load ' char(file_to_cluster) '_spikes.mat;']);
        
        %Filters and gets threshold
        [b,a] = ellip(2,0.1,40,[handles.par.sort_fmin handles.par.sort_fmax]*2/(handles.par.sr));
        xf = filtfilt(b,a,x);
        thr = handles.par.stdmin * median(abs(xf))/0.6745;
        thrmax = handles.par.stdmax * median(abs(xf))/0.6745;
    catch 
        continuous_data_av = 0; 
        % LOAD SC DATA
        eval(['load ' char(file_to_cluster) '_aligned ;'],'flag0=1;');
        if flag0==1
            eval(['load ' char(file_to_cluster) ]); 
        end;
    end; 

    % LOAD SPIKES
    handles.par.fname = ['data_' char(file_to_cluster)];   %filename for interaction with SPC
    nspk = size(spikes,1);
    naux = min(handles.par.max_spk,size(spikes,1));
    handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*naux);

    % CALCULATES INPUTS TO THE CLUSTERING ALGORITHM. 
    inspk = wave_features(spikes,handles);              %takes wavelet coefficients.
    
    % SELECTION OF SPIKES FOR SPC 
    if handles.par.permut == 'n'
        % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
        if size(spikes,1)> handles.par.max_spk;
            % take first 'handles.par.max_spk' spikes as an input for SPC
            inspk_aux = inspk(1:naux,:);
        else
            inspk_aux = inspk;
        end

        %INTERACTION WITH SPC
        save(handles.par.fname_in,'inspk_aux','-ascii');
        [clu, tree] = run_cluster(handles);
        [temp] = find_temp(tree,handles);

        %DEFINE CLUSTERS
        class1=find(clu(temp,3:end)==0);
        class2=find(clu(temp,3:end)==1);
        class3=find(clu(temp,3:end)==2);
        class4=find(clu(temp,3:end)==3);
        class5=find(clu(temp,3:end)==4);
        class0=setdiff(1:size(spikes,1), sort([class1 class2 class3 class4 class5]));
        whos class*
        
    else
        % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
        if size(spikes,1)> handles.par.max_spk;
            % random selection of spikes for SPC 
            ipermut = randperm(length(inspk));
            ipermut(naux+1:end) = [];
            inspk_aux = inspk(ipermut,:);
        else
            ipermut = randperm(length(inspk));
            inspk_aux = inspk(ipermut,:);
        end

        %INTERACTION WITH SPC
        save(handles.par.fname_in,'inspk_aux','-ascii');
        [clu, tree] = run_cluster(handles);
        [temp] = find_temp(tree,handles);

        %DEFINE CLUSTERS
        class1=ipermut(find(clu(temp,3:end)==0));
        class2=ipermut(find(clu(temp,3:end)==1));
        class3=ipermut(find(clu(temp,3:end)==2));
        class4=ipermut(find(clu(temp,3:end)==3));
        class5=ipermut(find(clu(temp,3:end)==4));
        class0=setdiff(1:size(spikes,1), sort([class1 class2 class3 class4 class5]));
        whos class*
    end
   
    % IF TEMPLATE MATCHING WAS DONE, THEN FORCE
    if (size(spikes,1)> handles.par.max_spk | ...
            (handles.par.force_auto == 'y'));
        classes = zeros(size(spikes,1),1);
        if length(class1)>=handles.par.min_clus; classes(class1) = 1; end
        if length(class2)>=handles.par.min_clus; classes(class2) = 2; end
        if length(class3)>=handles.par.min_clus; classes(class3) = 3; end
        if length(class4)>=handles.par.min_clus; classes(class4) = 4; end
        if length(class5)>=handles.par.min_clus; classes(class5) = 5; end
        f_in  = spikes(classes~=0,:);
        f_out = spikes(classes==0,:);
        class_in = classes(find(classes~=0),:);
        class_out = force_membership_wc(f_in, class_in, f_out, handles);
        classes(classes==0) = class_out;
        class0=find(classes==0);        
        class1=find(classes==1);        
        class2=find(classes==2);        
        class3=find(classes==3);        
        class4=find(classes==4);        
        class5=find(classes==5);        
    end    
        
    %PLOTS
    clf
    clus_pop = [];
    ylimit = [];
    subplot(3,5,11)
    temperature=handles.par.mintemp+temp*handles.par.tempstep;
    switch handles.par.temp_plot
            case 'lin'
                plot([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep], ...
            [handles.par.min_clus handles.par.min_clus],'k:',...
            handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
            tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
            case 'log'
                semilogy([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep], ...
                [handles.par.min_clus handles.par.min_clus],'k:',...
                handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
                tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
        end
    subplot(3,5,6)
    hold on
    cluster=zeros(nspk,2);
    cluster(:,2)= index';
    num_clusters = length(find([length(class1) length(class2) length(class3)...
            length(class4) length(class5) length(class0)] >= handles.par.min_clus));
    clus_pop = [clus_pop length(class0)];
    if length(class0) > handles.par.min_clus; 
        subplot(3,5,6); 
            max_spikes=min(length(class0),handles.par.max_spikes);
            plot(spikes(class0(1:max_spikes),:)','k'); 
            xlim([1 size(spikes,2)]);
        subplot(3,5,10); 
            hold on
            plot(spikes(class0(1:max_spikes),:)','k');  
            plot(mean(spikes(class0,:),1),'c','linewidth',2)
            xlim([1 size(spikes,2)]); 
            title('Cluster 0','Fontweight','bold')
        subplot(3,5,15)
            xa=diff(index(class0));
            [n,c]=hist(xa,0:1:100);
            bar(c(1:end-1),n(1:end-1))
            xlim([0 100])
            xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
            title([num2str(length(class0)) ' spikes']);
    end
    if length(class1) > handles.par.min_clus; 
        clus_pop = [clus_pop length(class1)];
        subplot(3,5,6); 
            max_spikes=min(length(class1),handles.par.max_spikes);
            plot(spikes(class1(1:max_spikes),:)','b'); 
            xlim([1 size(spikes,2)]);
        subplot(3,5,7); 
            hold
            plot(spikes(class1(1:max_spikes),:)','b'); 
            plot(mean(spikes(class1,:),1),'k','linewidth',2)
            xlim([1 size(spikes,2)]); 
            title('Cluster 1','Fontweight','bold')
            ylimit = [ylimit;ylim];
        subplot(3,5,12)
        xa=diff(index(class1));
        [n,c]=hist(xa,0:1:100);
        bar(c(1:end-1),n(1:end-1))
        xlim([0 100])
        set(get(gca,'children'),'facecolor','b','linewidth',0.01)    
        xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
        title([num2str(length(class1)) ' spikes']);
        cluster(class1(:),1)=1;
    end
    if length(class2) > handles.par.min_clus;
        clus_pop = [clus_pop length(class2)];
        subplot(3,5,6); 
            max_spikes=min(length(class2),handles.par.max_spikes);
            plot(spikes(class2(1:max_spikes),:)','r');  
            xlim([1 size(spikes,2)]);
        subplot(3,5,8); 
            hold
            plot(spikes(class2(1:max_spikes),:)','r');  
            plot(mean(spikes(class2,:),1),'k','linewidth',2)
            xlim([1 size(spikes,2)]); 
            title('Cluster 2','Fontweight','bold')
            ylimit = [ylimit;ylim];
        subplot(3,5,13)
            xa=diff(index(class2));
            [n,c]=hist(xa,0:1:100);
            bar(c(1:end-1),n(1:end-1))
            xlim([0 100])
            set(get(gca,'children'),'facecolor','r','linewidth',0.01)    
            xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
            cluster(class2(:),1)=2;
            title([num2str(length(class2)) ' spikes']);
    end
    if length(class3) > handles.par.min_clus;
        clus_pop = [clus_pop length(class3)];
        subplot(3,5,6); 
            max_spikes=min(length(class3),handles.par.max_spikes);
            plot(spikes(class3(1:max_spikes),:)','g');  
            xlim([1 size(spikes,2)]);
        subplot(3,5,9); 
            hold
            plot(spikes(class3(1:max_spikes),:)','g');  
            plot(mean(spikes(class3,:),1),'k','linewidth',2)
            xlim([1 size(spikes,2)]); 
            title('Cluster 3','Fontweight','bold')
            ylimit = [ylimit;ylim];
        subplot(3,5,14)
            xa=diff(index(class3));
            [n,c]=hist(xa,0:1:100);
            bar(c(1:end-1),n(1:end-1))
            xlim([0 100])
            set(get(gca,'children'),'facecolor','g','linewidth',0.01)    
            xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
            cluster(class3(:),1)=3;
            title([num2str(length(class3)) ' spikes']);
    end
    if length(class4) > handles.par.min_clus;
        clus_pop = [clus_pop length(class4)];
        subplot(3,5,6); 
            max_spikes=min(length(class4),handles.par.max_spikes);
            plot(spikes(class4(1:max_spikes),:)','c');  
            xlim([1 size(spikes,2)]);
        cluster(class4(:),1)=4;
    end
    if length(class5) > handles.par.min_clus; 
        clus_pop = [clus_pop length(class5)];
        subplot(3,5,6); 
            max_spikes=min(length(class5),handles.par.max_spikes);
            plot(spikes(class5(1:max_spikes),:)','m');  
            xlim([1 size(spikes,2)]);
        cluster(class5(:),1)=5;
    end

% Rescale spike's axis 
    if ~isempty(ylimit)
        ymin = min(ylimit(:,1));
        ymax = max(ylimit(:,2));
    else
        ymin = -200;
        ymax = 200;
    end
    if length(class1) > handles.par.min_clus; subplot(3,5,7); ylim([ymin ymax]); end
    if length(class2) > handles.par.min_clus; subplot(3,5,8); ylim([ymin ymax]); end
    if length(class3) > handles.par.min_clus; subplot(3,5,9); ylim([ymin ymax]); end
    if length(class0) > handles.par.min_clus; subplot(3,5,10); ylim([ymin ymax]); end
        
    subplot(3,1,1)
    box off; hold on
    %% these lines are for plotting continuous data; CPG 3/4/2008
    if continuous_data_av == 1
        plot((1:length(xf))/handles.par.sr,xf(1:length(xf)))
        line([0 length(xf)/handles.par.sr],[thr thr],'color','r')
        ylim([-thrmax/2 thrmax])
    end;
    % end of continuous data plotting; CPG 3/4/2008
    title([pwd '/' char(file_to_cluster)],'Interpreter','none','Fontsize',14)
    features_name = handles.par.features;
    toc
    if print2file==0;
        print
    else       
        set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print -djpeg fig2print_' char(file_to_cluster)]);
    end
        
    %SAVE FILES
    par = handles.par;
    cluster_class = cluster;
    outfile=['times_' char(file_to_cluster)];
    
    if handles.par.permut == 'n'
        save(outfile, 'cluster_class', 'par', 'spikes', 'inspk')
    else
        save(outfile, 'cluster_class', 'par', 'spikes', 'inspk', 'ipermut')
    end
    
    numclus=length(clus_pop)-1;
    outfileclus='cluster_results.txt';
    fout=fopen(outfileclus,'at+');
    fprintf(fout,'%s\t %s\t %g\t %d %g\t', char(file_to_cluster), features_name, temperature, numclus, handles.par.stdmin);
    for ii=1:numclus
        fprintf(fout,'%d\t',clus_pop(ii));
    end
    fprintf(fout,'%d\n',clus_pop(end));
    fclose(fout);
end
clear all

