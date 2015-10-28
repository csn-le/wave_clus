function Do_clustering(input, parallel, par_input)

% PROGRAM Do_clustering.
% Does clustering on all files in Files.txt
% Runs after Get_spikes.

% function Do_clustering(input, par_input)
% Saves spikes, spike times (in ms), coefficients used (inspk), used 
%parameters, random spikes selected for clustering (ipermut), spikes forced
%in a class (forced) results (cluster_class)

%input must be: 
%               A .txt file with the names of the spikes files to use.
%               A matlab cell with the names of the spikes files to use.
%               A vector, in this case the function will proccessall the
%                   '_spikes.mat' files with that numbers in the folder.
%                   (ipunt=2 don't implies 20 or viceversa)
%               'all', in this case the functions will process all the
%                '_spikes.mat' files in the folder.
%par_input must be a struct with some of the detecction parameters. All the
%parameters included will overwrite the parameters load from set_parameters()
if exist('parallel','var') && parallel == true
    if exist('matlabpool','file')
        matlabpool('open')
    else
        parpool
    end
end

if isnumeric(input) || any(strcmp(input,'all'))
    filenames = {};
    dirnames = dir();
    dirnames = {dirnames.name};
    
    for i = 1:length(dirnames)
        fname = dirnames{i};
        
        if length(fname) < 12 
            continue
        end
        if ~ strcmp(fname(end-10:end),'_spikes.mat')
            continue
        end
        if strcmp(input,'all')
            filenames = [filenames {fname}];
        else
            aux = regexp(fname, '\d+', 'match');
            if ismember(str2num(aux{1}),input)
                filenames = [filenames {fname}];   
            end
        end
    end
    
elseif ischar(input) && length(input) > 4
    if  strcmp (input(end-3:end),'.txt')
        filenames =  textread(input,'%s');
    else
        filenames = {input};
    end

elseif iscellstr(input)
    filenames = input;
else
    ME = MException('MyComponent:noValidInput', 'Invalid input arguments');
    throw(ME)
end


min_spikes4SPC = 16;

par_file = set_parameters();

for fnum = 1:length(filenames)
    filename = filenames{fnum};
    par = struct;
    par = update_parameters(par,par_file,'clus');
    
    if isfield(par,'channels')
        par.inputs = par.inputs * par.channels;
    end
    par.filename = filename;
    par.reset_results = true;
    
    data_handler = readInData(par);
    par = data_handler.par;
    
    if isfield(par,'channels')
        par.inputs = par.inputs * par.channels;
    end
    
    par.fname_in = 'tmp_data_wc';                       % temporary filename used as input for SPC
    par.fname = ['data_' data_handler.nick_name];
    par.nick_name = data_handler.nick_name;
    par.fnamespc = par.fname;                  		%filename if "save clusters" button is pressed

    if exist('par_input','var')
        par = update_parameters(par,par_input,'clus');
    end
    
    if data_handler.with_spikes            			%data have some time of _spikes files
        [spikes, index] = data_handler.load_spikes(); 
    else
        warning('MyComponent:noValidInput', 'File: %s doesn''t include spikes', filename);
        throw(ME)
        continue
    end
        
    % LOAD SPIKES
    nspk = size(spikes,1);
    naux = min(par.max_spk,size(spikes,1));
    par.min_clus = max(par.min_clus,par.min_clus_rel*naux);
    
	
    if nspk < min_spikes4SPC     
        warning('MyComponent:noValidInput', 'Not enough spikes in the file');
        continue
    end
    
    % CALCULATES INPUTS TO THE CLUSTERING ALGORITHM. 
    inspk = wave_features(spikes,par);     %takes wavelet coefficients.

    % SELECTION OF SPIKES FOR SPC 
    if par.permut == 'n'
        % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
        if size(spikes,1)> par.max_spk;
            % take first 'par.max_spk' spikes as an input for SPC
            inspk_aux = inspk(1:naux,:);
        else
            inspk_aux = inspk;
        end

        %INTERACTION WITH SPC
        save(par.fname_in,'inspk_aux','-ascii');
        [clu, tree] = run_cluster(par);
        [temp] = find_temp(tree,par);

        %DEFINE CLUSTERS
        class1 = find(clu(temp,3:end)==0);
        class2 = find(clu(temp,3:end)==1);
        class3 = find(clu(temp,3:end)==2);
        class4 = find(clu(temp,3:end)==3);
        class5 = find(clu(temp,3:end)==4);

    else
        % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
        if size(spikes,1)> par.max_spk;
            % random selection of spikes for SPC 
            ipermut = randperm(length(inspk));
            ipermut(naux+1:end) = [];
            inspk_aux = inspk(ipermut,:);
        else
            ipermut = randperm(size(inspk,1));
            inspk_aux = inspk(ipermut,:);
        end

        %INTERACTION WITH SPC
        save(par.fname_in,'inspk_aux','-ascii');
        [clu, tree] = run_cluster(par);
        [temp] = find_temp(tree,par);

        %DEFINE CLUSTERS
        class1=ipermut(clu(temp,3:end)==0);
        class2=ipermut(clu(temp,3:end)==1);
        class3=ipermut(clu(temp,3:end)==2);
        class4=ipermut(clu(temp,3:end)==3);
        class5=ipermut(clu(temp,3:end)==4);
        
    end
    
    class0 = setdiff(1:size(spikes,1), sort([class1 class2 class3 class4 class5]));
    whos class*
    
    
    % IF TEMPLATE MATCHING WAS DONE, THEN FORCE
    if (size(spikes,1)> par.max_spk || ...
            (par.force_auto))
        classes = zeros(size(spikes,1),1);
        if length(class1)>=par.min_clus; classes(class1) = 1; end
        if length(class2)>=par.min_clus; classes(class2) = 2; end
        if length(class3)>=par.min_clus; classes(class3) = 3; end
        if length(class4)>=par.min_clus; classes(class4) = 4; end
        if length(class5)>=par.min_clus; classes(class5) = 5; end
        f_in  = spikes(classes~=0,:);
        f_out = spikes(classes==0,:);
        class_in = classes(classes~=0,:);
        class_out = force_membership_wc(f_in, class_in, f_out, par);
        forced = classes==0;
        classes(classes==0) = class_out;
        forced(classes==0) =0;
        class0 = find(classes==0);
        class1 = find(classes==1);
        class2 = find(classes==2);
        class3 = find(classes==3);
        class4 = find(classes==4);
        class5 = find(classes==5);
        
    else
        forced = zeros(1, size(spikes,1));
    end
    current_par = par;
    par = struct;
    par = update_parameters(par, current_par, 'relevant');
    par.min_clus_rel = current_par.min_clus_rel;
    cluster = zeros(nspk,2);
    cluster(:,2)= index';
    for i = 1:5
        eval(['cluster(class' num2str(i) '(:),1)=' num2str(i) ';']);
    end
    cluster_class = cluster;
 
    
    temp_used = temp;
    for i = 1:5
       if eval(['length(class' num2str(i) ')']) > par.min_clus
            Temp(i) = temp_used;
       end
    end
    
    save(['times_' data_handler.nick_name], 'cluster_class','spikes', 'index', 'par','inspk','forced','Temp')
    if exist('ipermut','var')
        save(['times_' data_handler.nick_name],'ipermut','-append')
    end
end


for fnum = 1:length(filenames)
    filename = filenames{fnum};
    par = struct;
    par.filename = filename;

    par.cont_segment = true;  %maybe true and save the sample in spikes

    data_handler = readInData(par);
    par = data_handler.update_par(par);
    if ~data_handler.with_wc_spikes       			%data should have spikes
        continue
    end
    par.nick_name = data_handler.nick_name;

    figure('Visible','Off')
    set(gcf, 'PaperUnits', 'inches', 'PaperType', 'A4', 'PaperPositionMode', 'auto','units','normalized','outerposition',[0 0 1 1]) 
    subplot(3,1,1)
    if par.cont_segment && data_handler.with_psegment
        box off; hold on
        %these lines are for plotting continuous data 
        [xd_sub, sr_sub] = data_handler.get_signal_sample();
        lx = length(xd_sub);
        plot((1:lx)/sr_sub,xd_sub)
        noise_std_detect = median(abs(xd_sub))/0.6745;
        xlim([0 lx/sr_sub])
        thr = par.stdmin * noise_std_detect; 
        thrmax = thr*par.stdmax/par.stdmin;
        if strcmp(par.detection,'pos')
            line([0 length(xd_sub)/sr_sub],[thr thr],'color','r')
            ylim([-thrmax/2 thrmax])
        elseif strcmp(par.detection,'neg')
            line([0 length(xd_sub)/sr_sub],[-thr -thr],'color','r')
            ylim([-thrmax thrmax/2])
        else
            line([0 length(xd_sub)/sr_sub],[thr thr],'color','r')
            line([0 length(xd_sub)/sr_sub],[-thr -thr],'color','r')
            ylim([-thrmax thrmax])
        end
    end
    title([pwd '/' filename],'Interpreter','none','Fontsize',14)

    if ~data_handler.with_spc
       continue
    end
        
    % LOAD SPIKES        
    [clu, tree, spikes, index, inspk, ipermut, classes, forced] = data_handler.load_results();
    nspk = size(spikes,1);
    [temp] = find_temp(tree,par);
    class0 = find(classes==0);
    class1 = find(classes==1);
    class2 = find(classes==2);
    class3 = find(classes==3);
    class4 = find(classes==4);
    class5 = find(classes==5);  
    

    %PLOTS
    clus_pop = [];
    ylimit = [];
    subplot(3,5,11)
    temperature = par.mintemp+temp*par.tempstep;
    color = get(gca,'ColorOrder');
    hold on 
    switch par.temp_plot
            case 'lin'
                plot([par.mintemp par.maxtemp-par.tempstep], ...
                [par.min_clus par.min_clus],'k:',...
                par.mintemp+(1:par.num_temp)*par.tempstep, ...
                tree(1:par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
            
                for i=1:5
                    if eval(['length(class' num2str(i) ')']) > par.min_clus
                        tree_clus = tree(temp,4+i);
                        tree_temp = tree(temp+1,2);
                        plot(tree_temp,tree_clus,'.','color',color(i,:),'MarkerSize',20);
                    end
                end
            case 'log'
                semilogy([par.mintemp par.maxtemp-par.tempstep], ...
                [par.min_clus par.min_clus],'k:',...
                par.mintemp+(1:par.num_temp)*par.tempstep, ...
                tree(1:par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
                
                for i=1:5
                    if eval(['length(class' num2str(i) ')']) > par.min_clus
                        tree_clus = tree(temp,4+i);
                        tree_temp = tree(temp+1,2);
                        semilogy(tree_temp,tree_clus,'.','color',color(i,:),'MarkerSize',20);
                    end
                end
    end
    subplot(3,5,6)
    hold on
    
    clus_pop = [clus_pop length(class0)];
    if length(class0) > par.min_clus; 
        subplot(3,5,6); 
            max_spikes=min(length(class0),par.max_spikes_plot);
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
    if length(class1) > par.min_clus; 
        clus_pop = [clus_pop length(class1)];
        subplot(3,5,6); 
            max_spikes=min(length(class1),par.max_spikes_plot);
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
        
    end
    if length(class2) > par.min_clus;
        clus_pop = [clus_pop length(class2)];
        subplot(3,5,6); 
            max_spikes=min(length(class2),par.max_spikes_plot);
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
            title([num2str(length(class2)) ' spikes']);
    end
    if length(class3) > par.min_clus;
        clus_pop = [clus_pop length(class3)];
        subplot(3,5,6); 
            max_spikes=min(length(class3),par.max_spikes_plot);
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
            title([num2str(length(class3)) ' spikes']);
    end
    if length(class4) > par.min_clus;
        clus_pop = [clus_pop length(class4)];
        subplot(3,5,6); 
            max_spikes=min(length(class4),par.max_spikes_plot);
            plot(spikes(class4(1:max_spikes),:)','c');  
            xlim([1 size(spikes,2)]);
    end
    if length(class5) > par.min_clus; 
        clus_pop = [clus_pop length(class5)];
        subplot(3,5,6); 
            max_spikes=min(length(class5),par.max_spikes_plot);
            plot(spikes(class5(1:max_spikes),:)','m');  
            xlim([1 size(spikes,2)]);
    end

    % Rescale spike's axis 
    if ~isempty(ylimit)
        ymin = min(ylimit(:,1));
        ymax = max(ylimit(:,2));
    else
        ymin = -200;
        ymax = 200;
    end
    if length(class1) > par.min_clus; subplot(3,5,7); ylim([ymin ymax]); end
    if length(class2) > par.min_clus; subplot(3,5,8); ylim([ymin ymax]); end
    if length(class3) > par.min_clus; subplot(3,5,9); ylim([ymin ymax]); end
    if length(class0) > par.min_clus; subplot(3,5,10); ylim([ymin ymax]); end


    features_name = par.features;

    numclus=length(clus_pop)-1;
    outfileclus='cluster_results.txt';
    fout=fopen(outfileclus,'at+');
    if isfield(par,'stdmin')
        stdmin = par.stdmin;
    else
        stdmin = NaN;
    end
    fprintf(fout,'%s\t %s\t %g\t %d\t %g\t', char(filename), features_name, temperature, numclus, stdmin);
    for ii=1:numclus
        fprintf(fout,'%d\t',clus_pop(ii));
    end
    fprintf(fout,'%d\n',clus_pop(end));
    fclose(fout);

    
    
    if par.print2file;
        print(gcf,'-dpng',['fig2print_' filename(1:end-5) '.png'],'-r200');
    else
        print
    end 
    
end


if exist('parallel','var') && parallel == true
    if exist('matlabpool','file')
        matlabpool('close')
    else
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
 
end





	





    
	
