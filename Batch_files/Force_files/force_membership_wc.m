function class_out = force_membership_wc(f_in, class_in, f_out, par)
% class = function force_membership_wc(f_in, class_in, f_out, par)
% Given classified points, try to classify new points via template matching
%
% f_in:          features of classified points  (# input spikes x n_features)
% class_in:      classification of those points
% f_out:         features of points to be classified (nspk x n_features)
% par        environment variables, of which the following are
%                required: 
%                    o par.template_sdnum - max radius of cluster,
%                                                   in std devs.
%                    o par.template_k     - # of nearest neighbors
%                    o par.template_k_min - min # of nn for vote
%                    o par.template_type  - nn, center, ml, mahal

nspk = size(f_out,1);
class_out = zeros(1,size(f_out,1));
switch par.template_type
    case 'nn'
        sdnum = par.template_sdnum;
        k     = par.template_k;
        k_min = par.template_k_min;
        sd    = sqrt(sum(var(f_in,1)))*ones(1,size(f_in,1));
        for i=1:nspk,
            nn = nearest_neighbor(f_out(i,:),f_in,sdnum*sd,Inf*ones(size(f_in)),Inf,k);
            if( nn )
                winner = mode(class_in(nn));
                if nnz(class_in(nn)==winner)<k_min
                    class_out(i) = 0;
                else
                    class_out(i) = winner;
                end
            else
                class_out(i) = 0;
            end
        end
      
    case 'center'
        [centers, sd, pd] = build_templates(class_in,f_in); % we are going to ignore pd
        sdnum = par.template_sdnum;
        for i=1:nspk,
            class_out(i) = nearest_neighbor(f_out(i,:),centers,sdnum*sd);        
        end
        
    case 'ml'
        [mu inv_sigma] = fit_gaussian(f_in,class_in);
        for i=1:nspk,
            class_out(i) = ML_gaussian(f_out(i,:),mu,inv_sigma);
        end
    case 'mahal'
        classes = unique(class_in);
        mdistance = zeros(length(classes), nspk);
        maxdist   = zeros(1, length(classes));
        for ci = 1:length(classes)
           i = classes(ci);
           mdistance(i,:) = mahal(f_out, f_in(class_in ==i, :));
           maxdist(i) = sqrt(mean(mahal(f_in(class_in ==i, :), f_in(class_in ==i, :))));
        end
        sdnum = par.template_sdnum;
        for i = 1:nspk
             [d winner] = min(mdistance(:,i));
             if sqrt(d) < sdnum*maxdist(winner)
                 class_out(i) = classes(winner);
             end
        end
        
        
    otherwise
        sprintf('force_membership(): <%s> is not a known template type.\n',par.template_type);
        
end