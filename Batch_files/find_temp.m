function [temp] = find_temp(tree, par)
% Selects the temperature.

num_temp = floor((par.maxtemp-par.mintemp)/par.tempstep);     % total number of temperatures 
min_clus = par.min_clus;

aux = diff(tree(:,5));   % Changes in the first cluster size
aux1 = diff(tree(:,6));   % Changes in the second cluster size
aux2 = diff(tree(:,7));   % Changes in the third cluster size
aux3 = diff(tree(:,8));   % Changes in the third cluster size

temp = 1;         % Initial value

for t = 1:num_temp-1;
    % Looks for changes in the cluster size of any cluster larger than min_clus.
    if ( aux(t) > min_clus || aux1(t) > min_clus || aux2(t) > min_clus || aux3(t) >min_clus )    
        temp = t+1;         
    end
end

%In case the second cluster is too small, then raise the temperature a little bit 
if (temp == 1 && tree(temp,6) < min_clus)
    temp = 2;
end    
   