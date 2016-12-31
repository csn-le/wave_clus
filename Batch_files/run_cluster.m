function [clu, tree] = run_cluster(par, multi_files)
dim = par.inputs;
fname = par.fnamespc;
fname_in = par.fname_in;

% DELETE PREVIOUS FILES
if exist([fname '.dg_01.lab'],'file')
    delete([fname '.dg_01.lab']);
    delete([fname '.dg_01']);
end

dat = load(fname_in);
n = length(dat);
fid = fopen(sprintf('%s.run',fname),'wt');
fprintf(fid,'NumberOfPoints: %s\n',num2str(n));
fprintf(fid,'DataFile: %s\n',fname_in);
fprintf(fid,'OutFile: %s\n',fname);
fprintf(fid,'Dimensions: %s\n',num2str(dim));
fprintf(fid,'MinTemp: %s\n',num2str(par.mintemp));
fprintf(fid,'MaxTemp: %s\n',num2str(par.maxtemp));
fprintf(fid,'TempStep: %s\n',num2str(par.tempstep));
fprintf(fid,'SWCycles: %s\n',num2str(par.SWCycles));
fprintf(fid,'KNearestNeighbours: %s\n',num2str(par.KNearNeighb));
fprintf(fid,'MSTree|\n');
fprintf(fid,'DirectedGrowth|\n');
fprintf(fid,'SaveSuscept|\n');
fprintf(fid,'WriteLables|\n');
fprintf(fid,'WriteCorFile~\n');
if par.randomseed ~= 0
    fprintf(fid,'ForceRandomSeed: %s\n',num2str(par.randomseed));
end    
fclose(fid);

system_type = computer;
switch system_type
    case {'PCWIN'}    
%         if exist([pwd '\cluster.exe'])==0
%             directory = which('cluster.exe');
%             copyfile(directory,pwd);
%         end
        [status,result] = dos(sprintf('"%s" %s.run',which('cluster.exe'),fname));
        %[status,result] = dos(sprintf('cluster.exe %s.run',fname));
    case {'PCWIN64'}    
%         if exist([pwd '\cluster_64.exe'])==0
%             directory = which('cluster_64.exe');
%             copyfile(directory,pwd);
%         end
        [status,result] = dos(sprintf('"%s" %s.run',which('cluster_64.exe'),fname));
        %[status,result] = dos(sprintf('cluster_64.exe %s.run',fname));
    case {'MAC'}
        if exist([pwd '/cluster_mac.exe'])==0
            directory = which('cluster_mac.exe');
            copyfile(directory,pwd);
        end
        fileattrib([pwd '/cluster_mac.exe'],'+x')
        run_mac = sprintf('./cluster_mac.exe %s.run',fname);
	    [status,result] = unix(run_mac);
   case {'MACI','MACI64'}
        if exist([pwd '/cluster_maci.exe'])==0
            directory = which('cluster_maci.exe');
            copyfile(directory,pwd);
        end
        fileattrib([pwd '/cluster_maci.exe'],'+x')
        run_maci = sprintf('./cluster_maci.exe %s.run',fname);
	    [status,result] = unix(run_maci);
    case {'GLNX86'}      
        
        run_linux = sprintf('''%s'' %s.run',which('cluster_linux.exe'),fname);
        fileattrib(which('cluster_linux.exe'),'+x')
        
        [stat,mess]=fileattrib(which('cluster_linux.exe'));
        
        if mess.UserExecute==0
            if exist([pwd '/cluster_linux.exe'],'file') == 0
                directory = which('cluster_linux.exe');
                copyfile(directory,pwd);
            end
            run_linux = sprintf('./cluster_linux.exe %s.run',fname);
        end
	    [status,result] = unix(run_linux);
        
    case {'GLNXA64', 'GLNXI64'}
       
        run_linux = sprintf('''%s'' %s.run',which('cluster_linux64.exe'),fname);
        fileattrib(which('cluster_linux64.exe'),'+x')
        
        [stat,mess]=fileattrib(which('cluster_linux64.exe'));
        
        if mess.UserExecute==0
            if exist([pwd '/cluster_linux64.exe'],'file') == 0
                directory = which('cluster_linux64.exe');
                copyfile(directory,pwd);
            end
            run_linux = sprintf('./cluster_linux64.exe %s.run',fname);
        end
        [status,result] = unix(run_linux);
    otherwise 
    	ME = MException('MyComponent:NotSupportedArq', '%s type of computer not supported.',com_type);
    	throw(ME)
end

if status ~= 0
    disp(result)
end



if exist('multi_files','var') && multi_files==true
	log_name = [par.filename 'spc_log.txt'];
	f = fopen(log_name,'w');
	fprintf(f,['----------\nSPC result of file: ' par.filename '\n']);
	fprintf(f,result);
	fclose(f);
else
	log_name = 'spc_log.txt';
	f = fopen(log_name,'w');
	fprintf(f,result);
	fclose(f);
end

clu = load([fname '.dg_01.lab']);
tree = load([fname '.dg_01']); 

try
delete(sprintf('%s.run',fname));    
delete([fname '*.mag']);
delete([fname '*.edges']);
delete([fname '*.param']);
end

if exist([fname '.knn'],'file')
    delete([fname '.knn']);
end

delete(fname_in); 
