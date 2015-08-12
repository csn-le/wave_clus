function class_out = force_membership_wc(f_in, class_in, f_out, handles);
% class = function force_membership(f_in, class_in, f_out, handles)
% Given classified points, try to classify new points via template matching
%
% f_in:          features of classified points  (# input spikes x n_features)
% class_in:      classification of those points
% f_out:         features of points to be classified (nspk x n_features)
% handles        environment variables, of which the following are
%                required: 
%                    o handles.par.template_sdnum - max radius of cluster,
%                                                   in std devs.
%                    o handles.par.template_k     - # of nearest neighbors
%                    o handles.par.template_k_min - min # of nn for vote
%                    o handles.par.template_type  - nn, center, ml, mahal

nspk = size(f_out,1);
class_out = zeros(1,size(f_out,1));
switch handles.par.template_type
    case 'nn'
        sdnum = handles.par.template_sdnum;
        k     = handles.par.template_k;
        k_min = handles.par.template_k_min;
        sd    = sqrt(sum(var(f_in,1)))*ones(1,size(f_in,1));
        for i=1:nspk,
            nn = nearest_neighbor(f_out(i,:),f_in,sdnum*sd,Inf*ones(size(f_in)),Inf,k);
            if( nn )
                winner = vote(class_in(nn),k_min);
                class_out(i) = winner;
            else
                class_out(i) = 0;
            end
        end
      
    case 'center'
        [centers, sd, pd] = build_templates(class_in,f_in); % we are going to ignore pd
        sdnum = handles.par.template_sdnum;
        for i=1:nspk,
            class_out(i) = nearest_neighbor(f_out(i,:),centers,sdnum*sd);        
        end
        
    case 'ml'
        [mu sigma] = fit_gaussian(f_in,class_in);
        for i=1:nspk,
            class_out(i) = ML_gaussian(f_out(i,:),mu,sigma);
        end
    case 'mahal'
        [mu sigma] = fit_gaussian(f_in,class_in);
        for i=1:nspk,
            class_out(i) = nearest_mahal(f_out(i,:),mu,sigma);
        end
        
    otherwise
        sprintf('force_membership(): <%s> is not a known template type.\n',handles.par.template_type);
        
end