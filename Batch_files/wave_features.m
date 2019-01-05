function [inspk] = wave_features(spikes,par)
%Calculates the spike features

scales = par.scales;
feature = par.features;
nspk = size(spikes,1);
ls = size(spikes,2);
if strcmp(par.max_inputs,'all')
    par.max_inputs = ls;
elseif par.max_inputs < 1
    par.max_inputs = ceil(par.max_inputs*ls);
end

% CALCULATES FEATURES
switch feature
    case 'wav'
		cc=zeros(nspk,ls);
		try
			spikes_l = reshape(spikes',numel(spikes),1);
			if exist('wavedec')
				[c_l,l_wc] = wavedec(spikes_l,scales,'haar');
			else
				[c_l,l_wc] = fix_wavedec(spikes_l,scales);
			end
			wv_c = [0;l_wc(1:end-1)];
			nc = wv_c/nspk;
			wccum = cumsum(wv_c);
			nccum = cumsum(nc);
			for cf = 2:length(nc)
				cc(:,nccum(cf-1)+1:nccum(cf)) = reshape(c_l(wccum(cf-1)+1:wccum(cf)),nc(cf),nspk)';
			end
		catch
		    if exist('wavedec')                             % Looks for Wavelets Toolbox
				for i=1:nspk                                % Wavelet decomposition
					[c,l] = wavedec(spikes(i,:),scales,'haar');
					cc(i,1:ls) = c(1:ls);
				end
			else
				for i=1:nspk                                % Replaces Wavelets Toolbox, if not available
					[c,l] = fix_wavedec(spikes(i,:),scales);
					cc(i,1:ls) = c(1:ls);
				end
			end
		
		end
		
        for i=1:ls
            thr_dist = std(cc(:,i)) * 3;
            thr_dist_min = mean(cc(:,i)) - thr_dist;
            thr_dist_max = mean(cc(:,i)) + thr_dist;
            aux = cc(find(cc(:,i)>thr_dist_min & cc(:,i)<thr_dist_max),i);

            if length(aux) > 10
				          ks(i) = test_ks(aux);
            else
				          ks(i) = 0;
            end
        end

        [A,ind] = sort(ks);
        A = A(length(A)-par.max_inputs+1:end);
        ncoeff = length(A);
        maxA = max(A);
        nd = 10;
        d = (A(nd:end)-A(1:end-nd+1))/maxA*ncoeff/nd;
        all_above1 = find(d>=1);
	if numel(all_above1) >=2
		%temp_bla = smooth(diff(all_above1),3);
		aux2 = diff(all_above1);
		temp_bla = conv(aux2(:),[1 1 1]/3);
		temp_bla = temp_bla(2:end-1);
		temp_bla(1) = aux2(1);
		temp_bla(end) = aux2(end);

		thr_knee_diff = all_above1(find(temp_bla(2:end)==1,1))+(nd/2); %ask to be above 1 for 3 consecutive coefficients
		inputs = par.max_inputs-thr_knee_diff+1;
	else
		inputs = par.min_inputs;
	end


        if  isfield(par,'plot_feature_stats') && par.plot_feature_stats
            [path,name,ext] = fileparts(par.filename);
            if isempty(path)
                path='.';
            end
            fig = figure('visible','off');
            stairs(sort(ks))
            hold on
            ylabel('ks_stat','interpreter','none')
            xlabel('#features')
        			if ~isempty(inputs)
          				line([numel(ks)-inputs+1 numel(ks)-inputs+1],ylim,'color','r')
        			end
            line([numel(ks)-par.max_inputs numel(ks)-par.max_inputs],ylim,'LineStyle','--','color','k')
            title(sprintf('%s \n number of spikes = %d.  inputs_selected = %d.',name,nspk,inputs),'interpreter','none');
            print(fig,'-dpng',[path filesep 'feature_select_' name '.png'])
            close(fig)
        end

        if inputs > par.max_inputs
            inputs = par.max_inputs;
        elseif isempty(inputs) || inputs < par.min_inputs
            inputs = par.min_inputs;
        end

        coeff(1:inputs)=ind(ls:-1:ls-inputs+1);
    case 'pca'
        if exist('pca','file')
        	[C,S] = pca(spikes);
        else
            [C,S] = princomp(spikes);
        end
        cc = S;
        coeff = 1:size(S,2);
        warning('PCA uses 10 features')
        inputs = 10;
    case 'waveform'
        cc = spikes;
        coeff = 1:ls;
        inputs = ls;
end

%CREATES INPUT MATRIX FOR SPC
inspk=zeros(nspk,inputs);
for i=1:nspk
    for j=1:inputs
        inspk(i,j)=cc(i,coeff(j));
    end
end
