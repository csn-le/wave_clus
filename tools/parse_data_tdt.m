function parse_data_tdt(Block, max_memo_GB)
%This code required OpenDeveloper (Tucker-Davis Technologies)

folder = fileparts(Block);

if isempty(folder)
    Tank = pwd;
else
    Tank = folder;
end 

with_memory=true;
try
	memory;
catch
	with_memory=false;
end
if with_memory
	[userview,systemview] = memory;
	memo_avaible = floor(systemview.PhysicalMemory.Available*0.80);
	if exist('max_memo_GB','var') && ~isempty(max_memo_GB)
        max_memo = max_memo_GB*(1024)^3;
		if max_memo > memo_avaible
			error('max_memo_GB > 80% of Physical Memory Available')
		end
	else
		max_memo = memo_avaible;
	end
else
	max_memo = max_memo_GB*(1024)^3;
end


TTX = actxcontrol('TTank.X', [0 0 0 0]);
TTX.ConnectServer('Local','Me');

if TTX.ConnectServer('Local','Me') == 0
    error('Error connecting to server');
end
if TTX.OpenTank(Tank,'R') == 0
    error('Error opening tank');
end
if TTX.SelectBlock(Block) == 0
    error('Error opening block');
end

TTX.SetGlobalV('WavesMemLimit',max_memo);
set(gcf,'visible','off')


N = TTX.ReadEventsV(10000, 'Wave', 0, 0, 0, 0, 'NODATA');
channels = unique(TTX.ParseEvInfoV(0, N, 4));
fprintf('%d channels in the block\n', max(channels))
sr = unique(TTX.ParseEvInfoV(0,N,9));
%timestamps = TTX.ParseEvInfoV(0,N,6);
 
for ch = channels

    N = TTX.ReadEventsV(1000000, 'Wave', ch, 0, 0, 0, 'ALL');
    if N == 1000000
        N = TTX.ReadEventsV(N*10000, 'Wave', ch, 0, 0, 0, 'ALL');
    end
    
    fout = fopen([Block '_' num2str(ch) '.tdtch'],'w','l');
    data = TTX.ParseEvV(0,1);
    rec_readed = 1;
    rec_to_read = N - 1;
    smp_per_rec = length(data);
    B_per_rec = whos('data');
    B_per_rec = B_per_rec.bytes;
    loops = ceil(B_per_rec/max_memo);
    rec_per_loop = ceil(rec_to_read/loops);
    fwrite(fout,data,'single');
    for seg = 1:loops
        Nreading = min(rec_per_loop,rec_to_read);
		%TTX.ReadEventsV(N, 'Wave', ch, 0, 0, 0, 'ALL');
        data = TTX.ParseEvV(rec_readed, Nreading);
        data = reshape(data,1,[]);
        fwrite(fout,data,'single');
        rec_readed = rec_readed + Nreading;
        rec_to_read = rec_to_read - Nreading;
        clear data
    end
	fprintf('Channel %d out. \n', ch)
    fclose(fout);

    
end

TTX.CloseTank
TTX.ReleaseServer
close('all')

lts = N * smp_per_rec;
save('tdt_meta_data','sr','lts','channels')



