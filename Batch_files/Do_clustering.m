function Do_clustering(input)

% PROGRAM Do_clustering.
% Does clustering on all files in Files.txt
% Runs after Get_spikes.


if isnumeric(input) || strcmp(input,'all')
    filenames = {};
    se = supported_wc_extensions();
    dirnames = dir();
    dirnames = {dirnames.name};
    
    for i = 1:length(dirnames)
        fname = dirnames{i};
        [unused, f, ext] = fileparts(fname);
        ext = lower(ext(2:end));
        if any(strcmp(ext,se)) 
            if strcmp(ext,'mat')
                sprintf('Skipped file ''%s''. The ''.mat'' files should be added by name.\n',fname);
                continue
            end
            if strcmp(input,'all')
                filenames = [filenames {fname}];
            else
                aux = regexp(f, '\d+', 'match');
                if ismember(str2num(aux{1}),input)
                    filenames = [filenames {fname}];   
                end
            end
        end
    end
    
elseif ischar(input) && length(input) > 4
    if  strcmp (input(end-3,end),'.txt')
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



for i = 1: size(filenames,1)
    
    par = set_parameters();
    par.filename = filename;
    par.reset_results = true;

    par.sample_segment = true;  %maybe true and save the sample in spikes

    data_handler = readInData(par);
    par = data_handler.par;


    if data_handler.with_spikes            %data have some time of _spikes files
        [spikes, index] = data_handler.load_spikes(); 
        if ~data_handler.with_wc_spikes
            [spikes] = spike_alignment(spikes,par);
        end
    else    
        set(handles.file_name,'string','Detecting spikes ...'); drawnow
        index = [];
        spikes = [];
        for n = 1:data_handler.max_segments
            x = data_handler.get_segment();
                %<----  Add here extra processing of the signal (x)
            [new_spikes, temp_aux_th, new_index]  = amp_detect(x, handles);
            index = [index data_handler.index2ts(new_index)]; %new_index to ms
            spikes = [spikes; new_spikes];
        end
    end

       

        
        
    
    
    
    
    
    
    
    current_par = par;
    par = struct;
    par = update_parameters(par, current_par, 'relevant');

    %<----  Add here auxiliar parameters

    save(['times_' data_handler.nick_name], 'spikes', 'index', 'par')        
end
end









% LOAD SPIKES
nspk = size(spikes,1);
naux = min(handles.par.max_spk,size(spikes,1));
par.min_clus = max(par.min_clus_abs,par.min_clus_rel*naux);
    
	
if nspk < 16

end    
	
	
	
	
figure
set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]) 


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
	[clu, tree] = run_cluster(handles.par);
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
		ipermut = randperm(size(inspk,1));
		inspk_aux = inspk(ipermut,:);
	end

	%INTERACTION WITH SPC
	save(handles.par.fname_in,'inspk_aux','-ascii');
	[clu, tree] = run_cluster(handles.par);
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
		max_spikes=min(length(class0),handles.par.max_spikes_plot);
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
		max_spikes=min(length(class1),handles.par.max_spikes_plot);
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
		max_spikes=min(length(class2),handles.par.max_spikes_plot);
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
		max_spikes=min(length(class3),handles.par.max_spikes_plot);
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
		max_spikes=min(length(class4),handles.par.max_spikes_plot);
		plot(spikes(class4(1:max_spikes),:)','c');  
		xlim([1 size(spikes,2)]);
	cluster(class4(:),1)=4;
end
if length(class5) > handles.par.min_clus; 
	clus_pop = [clus_pop length(class5)];
	subplot(3,5,6); 
		max_spikes=min(length(class5),handles.par.max_spikes_plot);
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
	
%SAVE FILES
	par = handles.par;
	cluster_class = cluster;
	outfile=['times_' char(file_to_cluster)];

	if handles.par.permut == 'n'
		save(outfile, 'cluster_class', 'par', 'spikes', 'inspk')
	else
		save(outfile, 'cluster_class', 'par', 'spikes', 'inspk', 'ipermut')
	end
	features_name = handles.par.features;

	numclus=length(clus_pop)-1;
	outfileclus='cluster_results.txt';
	fout=fopen(outfileclus,'at+');
	fprintf(fout,'%s\t %s\t %g\t %d\t %g\t', char(file_to_cluster), features_name, temperature, numclus, handles.par.stdmin);
	for ii=1:numclus
		fprintf(fout,'%d\t',clus_pop(ii));
	end
	fprintf(fout,'%d\n',clus_pop(end));
	fclose(fout);
	clear inspk; clear cluster_class
	
	subplot(3,1,1)      
end
box off; hold on
%% these lines are for plotting continuous data 
if continuous_data_av == 1
	plot((1:length(xf))/handles.par.sr,xf(1:length(xf)))
	if strcmp(handles.par.detection,'pos')
		line([0 length(xf)/handles.par.sr],[thr thr],'color','r')
		ylim([-thrmax/2 thrmax])
	elseif strcmp(handles.par.detection,'neg')
		line([0 length(xf)/handles.par.sr],[-thr -thr],'color','r')
		ylim([-thrmax thrmax/2])
	else
		line([0 length(xf)/handles.par.sr],[thr thr],'color','r')
		line([0 length(xf)/handles.par.sr],[-thr -thr],'color','r')
		ylim([-thrmax thrmax])
	end
end;
title([pwd '/' char(file_to_cluster)],'Interpreter','none','Fontsize',14)
%     title([pwd '   Channel  ' num2str(channel)],'Interpreter','none','Fontsize',14)
if print2file==0;
	print
else
	set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
	set(gcf,'paperposition',[.25 .25 10.5 7.8])
%         eval(['print -djpeg fig2printNSX_ch' num2str(channel)]);
	eval(['print -dpng fig2print_' char(file_to_cluster)]);
end 

clear spikes; 

%     
%     subplot(3,1,1)
% 
%     box off; hold on
%     %% these lines are for plotting continuous data 
%     if continuous_data_av == 1
%         plot((1:length(xf))/handles.par.sr,xf(1:length(xf)))
%         if strcmp(handles.par.detection,'pos')
%             line([0 length(xf)/handles.par.sr],[thr thr],'color','r')
%             ylim([-thrmax/2 thrmax])
%         elseif strcmp(handles.par.detection,'neg')
%             line([0 length(xf)/handles.par.sr],[-thr -thr],'color','r')
%             ylim([-thrmax thrmax/2])
%         else
%             line([0 length(xf)/handles.par.sr],[thr thr],'color','r')
%             line([0 length(xf)/handles.par.sr],[-thr -thr],'color','r')
%             ylim([-thrmax thrmax])
%         end
%     end;
%     % end of continuous data plotting
%     title([pwd '/' char(file_to_cluster)],'Interpreter','none','Fontsize',14)
%     features_name = handles.par.features;
%     toc
%     if print2file==0;
%         print
%     else       
%         set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
%         set(gcf,'paperposition',[.25 .25 10.5 7.8])
%         eval(['print -djpeg fig2print_' char(file_to_cluster)]);
%     end
%         
%     %SAVE FILES
%     par = handles.par;
%     cluster_class = cluster;
%     outfile=['times_' char(file_to_cluster)];
%     
%     if handles.par.permut == 'n'
%         save(outfile, 'cluster_class', 'par', 'spikes', 'inspk')
%     else
%         save(outfile, 'cluster_class', 'par', 'spikes', 'inspk', 'ipermut')
%     end
%     
%     numclus=length(clus_pop)-1;
%     outfileclus='cluster_results.txt';
%     fout=fopen(outfileclus,'at+');
%     fprintf(fout,'%s\t %s\t %g\t %d %g\t', char(file_to_cluster), features_name, temperature, numclus, handles.par.stdmin);
%     for ii=1:numclus
%         fprintf(fout,'%d\t',clus_pop(ii));
%     end
%     fprintf(fout,'%d\n',clus_pop(end));
%     fclose(fout);
% clear spikes; 
