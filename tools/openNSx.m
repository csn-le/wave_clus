function varargout = openNSx(varargin)

% openNSx
% 
% Opens and reads an NSx file then returns all file information in a NSx
% structure. Works with File Spec 2.1, 2.2 and 2.3.
% Use OUTPUT = openNSx(fname, 'read', 'report', 'e:xx:xx', 'c:xx:xx', 't:xx:xx', 'mode', 'precision', 'skipfactor').
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   fname:        Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file. 
%                 DEFAULT: Will open Open File UI.
%
%   'read':       Will read the data in addition to the header information
%                 if user passes this argument.
%                 DEFAULT: will only read the header information.
%
%   'report':     Will show a summary report if user passes this argument.
%                 DEFAULT: will not show report.
%
%   'e:XX:YY':    User can specify which electrodes need to be read. The
%                 number of electrodes can be greater than or equal to 1
%                 and less than or equal to 128. The electrodes can be
%                 selected either by specifying a range (e.g. 20:45) or by
%                 indicating individual electrodes (e.g. 3,6,7,90) or both.
%                 Note that, when individual channels are to be read, all
%                 channels in between will also be read. The prorgam will
%                 then remove the unwanted channels. This may result in a
%                 large memory footprint. If memory issues arrise please
%                 consider placing openNSx in a for loop and reading
%                 individual channels.
%                 This field needs to be preceded by the prefix 'e:'. See
%                 example for more details. If this option is selected the
%                 user will be promped for a CMP mapfile (see: KTUEAMapFile)
%                 provided by Blackrock Microsystems. This feature required
%                 KTUEAMapFile to be present in path.
%                 DEFAULT: will read all existing electrodes.
%
%   'c:XX:YY':    User can specify which channels need to be read. The
%                 number of channels can be greater than or equal to 1
%                 and less than or equal to 255. The channels can be
%                 selected either by specifying a range (e.g. 20:45) or by
%                 indicating individual channels (e.g. 3,6,7,90) or both.
%                 Note that, when individual channels are to be read, all
%                 channels in between will also be read. The prorgam will
%                 then remove the unwanted channels. This may result in a
%                 large memory footprint. If memory issues arrise please
%                 consider placing openNSx in a for loop and reading
%                 individual channels.
%                 This field needs to be preceded by the prefix 'c:'. See
%                 example for more details.
%                 DEFAULT: will read all existing analog channels.
%
%   't:XX:YY':    User can specify the beginning and end of the data
%                 segment to be read. If the start time is greater than the
%                 length of data the program will exit with an error
%                 message. If the end time is greater than the length of
%                 data the end packet will be selected for end of data. The
%                 user can specify the start and end values by comma 
%                 (e.g. [20,50]) or by a colon (e.g. [20:50]). To use this
%                 argument the user must specify the [electrodes] or the
%                 interval will be used for [electrodes] automatically.
%                 This field needs to be preceded by the prefix 't:'. 
%                 Note that if 'mode' is 'sample' the start duration cannot
%                 be less than 1. The duration is inclusive.
%                 See example for more details.
%                 DEFAULT: will read the entire file.
%
%   'mode':       The user can specify the mode of duration in [duration],
%                 such as 'sec', 'min', 'hour', or 'sample'. If 'sec' is
%                 specified the numbers in [duration] will correspond to
%                 the number of seconds. The same is true for 'min', 'hour'
%                 and 'sample'.
%                 DEFAULT: reads 'sample'.
%
%   'precision':  This will specify the precision for NSx file. If set to
%                 'double' the NSx data will be read as 'double' and if set
%                 to 'short', the NSx data will be read as 'int16' data
%                 type. While reading the file as 'short' may have a much
%                 smaller memory footprint and a faster read time, some 
%                 post data analysis such as multiplying the signal by a 
%                 factor that will make the data larger than (-32,768 to 
%                 32,767 -- refer to MATLAB documentation for more 
%                 information) may result in unexpected behavior. 
%                 Always use caution when using short. If you are not sure
%                 of what to use then do not specify this option.
%                 DEFAULT: will read data in 'int16'.
%
%   'skipfactor': This option will allow the user to read a decimated
%                 version of the data. The skipfactor will determine how
%                 many samples to skip. For example, if skipfactor is 2
%                 then every other sample is read. If skipfactor is 5 then
%                 every fifth sample is read. This is useful to briefly
%                 looking at the data in a large datafile when reading the
%                 entire dataset would overflow memory.
%                 DEFAULT: is set to 1, so every sample will be read.
%
%   'ver':        If this argument is passed to the function it will return
%                 the version number of the function without reading any
%                 data files.
%
%   OUTPUT:       Contains the NSx structure.
%
%   Example 1: 
%   openNSx('report','read','c:\data\sample.ns5', 'e:15:30', 't:3:10','min', 'p:short', 's:5');
%
%   or equivalently
%   openNSx('report','read','c:\data\sample.ns5', 'electrodes', 15:30, 'duration', 3:10, 'min', 'precision', 'short', 'skipfactor', 5);
%
%   In the example above, the file c:\data\sample.ns5 will be used. A
%   report of the file contents will be shown. The data will be read from
%   electrodes 15 through 50 in the 3-10 minute time interval. A decimated 
%   version of the datafile will be read, where only every 5th sample is read.
%   If any of the arguments above are omitted the default values will be used.
%
%   Example 2:
%   openNSx('read','c:15:30');
%
%   In the example above, the file user will be prompted for the file. The
%   file will be read using 'int16' precision as default. All time points 
%   of Only channels 15 through 30 will be read. If any of the arguments 
%   above are omitted the default values will be used.
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   Version 5.1.6.0
%

%% Defining the NSx data structure and sub-branches.
NSx          = struct('MetaTags',[],'Data',[], 'RawData', []);
NSx.MetaTags = struct('FileTypeID',[],'SamplingLabel',[],'ChannelCount',[],'SamplingFreq',[], 'TimeRes', [], ...
                      'ChannelID',[],'DateTime',[],'DateTimeRaw',[], 'Comment', [], 'FileSpec', [], ...
                      'Timestamp', [], 'DataPoints', [], 'DataDurationSec', [], 'openNSxver', [], 'Filename', [], 'FilePath', [], ...
                      'FileExt', []);

NSx.MetaTags.openNSxver = '5.1.6.0';

% Defining constants
ExtHeaderLength = 66;
elecReading     = 0;
maxNSPChannels  = 128;
NSx.RawData.PausedFile = 0;

%% Validating the input arguments. Exit with error message if error occurs.
next = '';
for i=1:length(varargin)
    inputArgument = varargin{i};
    if strcmpi(inputArgument, 'ver')
        varargout{1} = NSx.MetaTags.openNSxver;
        return;
    elseif strcmpi(inputArgument, 'channels')
        next = 'channels';
    elseif strcmpi(inputArgument, 'skipfactor')
        next = 'skipfactor';
    elseif strcmpi(inputArgument, 'electrodes')
        next = 'electrodes';
    elseif strcmpi(inputArgument, 'duration')
        next = 'duration';
    elseif strcmpi(inputArgument, 'precision')
        next = 'precision';
    elseif strcmpi(inputArgument, 'report')
        Report = inputArgument;
    elseif strcmpi(inputArgument, 'noread')
        ReadData = inputArgument;
    elseif strcmpi(inputArgument, 'read')
        ReadData = inputArgument;
    elseif (strncmp(inputArgument, 't:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'duration')
        if strncmp(inputArgument, 't:', 2)
            inputArgument(1:2) = [];
        end
        modifiedTime = 1;
        inputArgument = str2num(inputArgument);
        StartPacket = inputArgument(1);
        EndPacket = inputArgument(end);
        next = '';
    elseif (strncmp(inputArgument, 'e:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'electrodes')
        if exist('KTUEAMapFile', 'file') == 2
            Mapfile = KTUEAMapFile;
            Elec = str2num(inputArgument(3:end)); %#ok<ST2NM>
            if min(Elec)<1 || max(Elec)>128
                disp('The electrode number cannot be less than 1 or greater than 128.');
                if nargout; varargout{1} = []; end
                return;
            end
            for chanIDX = 1:length(Elec)
                userRequestedChannels(chanIDX) = Mapfile.Electrode2Channel(Elec(chanIDX));
            end
            elecReading = 1;
        else
            disp('To read data by ''electrodes'' the function KTUEAMapFile needs to be in path.');
            clear variables;
            return;
        end
        next = '';
    elseif (strncmp(inputArgument, 's:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'skipFactor')
        if strncmp(inputArgument, 's:', 2)
            skipFactor = str2num(inputArgument(3:end)); %#ok<ST2NM>
        else
            skipFactor = inputArgument;
        end
        next = '';
    elseif (strncmp(inputArgument, 'c:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'channels')
        if strncmp(inputArgument, 'c:', 2)
            userRequestedChanRow = str2num(inputArgument(3:end)); %#ok<ST2NM>
        else
            userRequestedChanRow = inputArgument;
        end
        next = '';
    elseif (strncmp(varargin{i}, 'p:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'precision')
        if strncmp(varargin{i}, 'p:', 2)
            precisionTypeRaw = varargin{i}(3:end);
        else
            precisionTypeRaw = varargin{i};
        end
        switch precisionTypeRaw
			case 'int16'
				precisionType = '*int16=>int16';
            case 'short'
                precisionType = '*short=>short';
            case 'double'
                precisionType = '*int16';
            otherwise
                disp('Read type is not valid. Refer to ''help'' for more information.');
                if nargout; varargout{1} = []; end
                return;
        end
        clear precisionTypeRaw;
        next = '';
    elseif strfind(' hour min sec sample ', [' ' inputArgument ' ']) ~= 0
        TimeScale = inputArgument;
    else
        temp = inputArgument;
        if length(temp)>3 && strcmpi(temp(end-3),'.')
            fname = inputArgument;
            if exist(fname, 'file') ~= 2
                disp('The file does not exist.');
                if nargout; 
                    varargout{1} = []; 
                end
                return;
            end
        else
            disp(['Invalid argument ''' inputArgument ''' .']);
            if nargout; varargout{1} = []; end
            return;
        end
    end
end
clear next;

%% Popup the Open File UI. Also, process the file name, path, and extension
%  for later use, and validate the entry.
if ~exist('fname', 'var')
    if ~ismac
        [fname, path] = getFile('*.ns*', 'Choose an NSx file...');
    else
        [fname, path] = getFile('*.*', 'Choose an NSx file...');
    end
    if fname == 0
        disp('No file was selected.');
        if nargout
            clear variables;
        end
        return;
    end
    fext = fname(end-3:end);
else
    if isempty(fileparts(fname))
        fname = which(fname);
    end
    [path,fname, fext] = fileparts(fname);
    fname = [fname fext];
    path  = [path '/'];
end
if fname==0
    if nargout; varargout{1} = []; end
    return; 
end

tic;

%% Give all input arguments a default value. All input argumens are
%  optional.
if ~exist('Report', 'var');        Report = 'noreport';    end
if ~exist('ReadData', 'var');      ReadData = 'noread';    end
if ~exist('StartPacket', 'var');   StartPacket = 0;        end
if ~exist('TimeScale', 'var');     TimeScale = 'sample';   end
if ~exist('precisionType', 'var'); precisionType = '*short=>short'; end
if ~exist('skipFactor', 'var');    skipFactor = 1;         end
if ~exist('modifiedTime', 'var');  modifiedTime = 0;       end

if strcmpi(ReadData, 'noread')
    disp('NOTE: Reading the header information only. To read the data use with parameter ''read'': openNSx(''read'')');
end

if strcmp(Report, 'report')
    disp(['openNSx ' NSx.MetaTags.openNSxver]);
end

%% Reading Basic Header from file into NSx structure.
FID                       = fopen([path fname], 'r', 'ieee-le');
NSx.MetaTags.Filename     = fname;
NSx.MetaTags.FilePath     = path(1:end-1);
NSx.MetaTags.FileExt      = fext;
NSx.MetaTags.FileTypeID   = fread(FID, [1,8]   , '*char');
if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALSG')
	NSx.MetaTags.FileSpec      = '2.1';
    NSx.MetaTags.SamplingLabel = fread(FID, [1,16]  , '*char');
    NSx.MetaTags.TimeRes       = 30000;
    NSx.MetaTags.SamplingFreq  = NSx.MetaTags.TimeRes / fread(FID, 1 , 'uint32=>double');
    ChannelCount               = double(fread(FID, 1       , 'uint32=>double'));
    NSx.MetaTags.ChannelCount  = ChannelCount;
    NSx.MetaTags.ChannelID     = fread(FID, [ChannelCount 1], '*uint32');
    try
    	t                          = dir([path fname]);
    	NSx.MetaTags.DateTime      = t.date;
    end
elseif strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD')
    BasicHeader                = fread(FID, 306, '*uint8');
    NSx.MetaTags.FileSpec      = [num2str(double(BasicHeader(1))) '.' num2str(double(BasicHeader(2)))];
    HeaderBytes                = double(typecast(BasicHeader(3:6), 'uint32'));
    NSx.MetaTags.SamplingLabel = char(BasicHeader(7:22))';
    NSx.MetaTags.Comment       = char(BasicHeader(23:278))';
    NSx.MetaTags.TimeRes       = double(typecast(BasicHeader(283:286), 'uint32'));
    NSx.MetaTags.SamplingFreq  = NSx.MetaTags.TimeRes / double(typecast(BasicHeader(279:282), 'uint32'));
    t                          = double(typecast(BasicHeader(287:302), 'uint16'));
    ChannelCount               = double(typecast(BasicHeader(303:306), 'uint32'));
    NSx.MetaTags.ChannelCount  = ChannelCount;
    readSize                   = double(ChannelCount * ExtHeaderLength);
    ExtendedHeader             = fread(FID, readSize, '*uint8');
    
    %% Removing extra garbage characters from the Comment field.
    NSx.MetaTags.Comment(find(NSx.MetaTags.Comment==0,1):end) = 0;    
    
	%% Populating extended header information
	for headerIDX = 1:ChannelCount
		offset = double((headerIDX-1)*ExtHeaderLength);
		NSx.ElectrodesInfo(headerIDX).Type = char(ExtendedHeader((1:2)+offset))';
		if (~strcmpi(NSx.ElectrodesInfo(headerIDX).Type, 'CC'))
			disp('extended header not supported');
			fclose(FID);
			if nargout; varargout{1} = []; end
			return;			
		end
		NSx.ElectrodesInfo(headerIDX).ElectrodeID = typecast(ExtendedHeader((3:4)+offset), 'uint16');
		NSx.ElectrodesInfo(headerIDX).Label = char(ExtendedHeader((5:20)+offset))';
		NSx.ElectrodesInfo(headerIDX).ConnectorBank = char(ExtendedHeader(21+offset) + ('A' - 1));
		NSx.ElectrodesInfo(headerIDX).ConnectorPin   = ExtendedHeader(22+offset);
		NSx.ElectrodesInfo(headerIDX).MinDigiValue   = typecast(ExtendedHeader((23:24)+offset), 'int16');
		NSx.ElectrodesInfo(headerIDX).MaxDigiValue   = typecast(ExtendedHeader((25:26)+offset), 'int16');
		NSx.ElectrodesInfo(headerIDX).MinAnalogValue = typecast(ExtendedHeader((27:28)+offset), 'int16');
		NSx.ElectrodesInfo(headerIDX).MaxAnalogValue = typecast(ExtendedHeader((29:30)+offset), 'int16');
		NSx.ElectrodesInfo(headerIDX).AnalogUnits    = char(ExtendedHeader((31:46)+offset))';
		NSx.ElectrodesInfo(headerIDX).HighFreqCorner = typecast(ExtendedHeader((47:50)+offset), 'uint32');
		NSx.ElectrodesInfo(headerIDX).HighFreqOrder  = typecast(ExtendedHeader((51:54)+offset), 'uint32');
		NSx.ElectrodesInfo(headerIDX).HighFilterType = typecast(ExtendedHeader((55:56)+offset), 'uint16');
		NSx.ElectrodesInfo(headerIDX).LowFreqCorner  = typecast(ExtendedHeader((57:60)+offset), 'uint32');
		NSx.ElectrodesInfo(headerIDX).LowFreqOrder   = typecast(ExtendedHeader((61:64)+offset), 'uint32');
		NSx.ElectrodesInfo(headerIDX).LowFilterType  = typecast(ExtendedHeader((65:66)+offset), 'uint16');
	end
	clear ExtendedHeader;
	%% Parsing and validating FileSpec and DateTime variables
	NSx.MetaTags.DateTimeRaw = t.';
	NSx.MetaTags.DateTime = [num2str(t(2)) '/'  num2str(t(4)) '/' num2str(t(1))...
		' ' datestr(t(3), 'dddd') ' ' num2str(t(5)) ':'  ...
		num2str(t(6)) ':'  num2str(t(7)) '.' num2str(t(8))] ;
	clear t;
else
    disp('This version of openNSx can only read File Specs 2.1, 2.2 and 2.3');
    disp(['The selected file spec is ' NSx.MetaTags.FileSpec '.']);
    fclose(FID);
    clear variables;
    return;
end

% Determining the length of file and storing the value of fEOF
f.EOexH = double(ftell(FID));
fseek(FID, 0, 'eof');
f.EOF = double(ftell(FID));

% Read Raw Header for saveNSx
fseek(FID, 0, 'bof');
NSx.RawData.Headers = fread(FID, f.EOexH, 'uint8');
NSx.RawData.DataHeader = fread(FID, 9, 'uint8');
fseek(FID, f.EOexH, 'bof');

% Reading all data headers and calculating all the file pointers for data
% and headers
if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALSG')
    % Determining DataPoints
    f.BOData = f.EOexH;
    f.EOData = f.EOF;
    NSx.MetaTags.DataPoints = (f.EOF-f.EOexH)/(ChannelCount*2);
elseif strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD')    
    segmentCount = 0;
    while double(ftell(FID)) < f.EOF
        if (fread(FID, 1, 'uint8') ~= 1)
            % Fixing another bug in Central 6.01.00.00 TOC where DataPoints is
            % not written back into the Data Header
            %% BIG NEEDS TO BE FIXED
            NSx.MetaTags.DataPoints = double(f.EOF - f.BOData)/(ChannelCount*2);
            break;
        end
        segmentCount = segmentCount + 1;
        NSx.MetaTags.Timestamp(segmentCount)  = fread(FID, 1, 'uint32');
        NSx.MetaTags.DataPoints(segmentCount) = fread(FID, 1, 'uint32');
        f.BOData(segmentCount) = double(ftell(FID));
        fseek(FID, NSx.MetaTags.DataPoints(segmentCount) * ChannelCount * 2, 'cof');
        f.EOData(segmentCount) = double(ftell(FID));
        % Fixing the bug in 6.01.00.00 TOC where DataPoints is not
        % updated and is left as 0
        NSx.MetaTags.DataPoints(segmentCount) = (f.EOData(segmentCount)-f.BOData(segmentCount))/(ChannelCount*2);
    end
end

% Fixing a bug in 6.03.00.00 TOC where an extra data packet (length 9) was
% written for no reason. Removing the information read for the extra
% invalid packet
if length(NSx.MetaTags.DataPoints) > 1 && all(NSx.MetaTags.Timestamp(1:2) == [0,0])
    NSx.MetaTags.DataPoints(1) = [];
    NSx.MetaTags.Timestamp(1) = [];
    f.BOData(1) = [];
    f.EOData(1) = [];
    segmentCount = 1;
end

% Determining if the file has a pause in it
if length(NSx.MetaTags.DataPoints) > 1
    NSx.RawData.PausedFile = 1;
    if modifiedTime == 1
        disp('This data file contains pauses.');
        disp('openNSx cannot read files with pauses using the ''t:XX'' parameter.');
        fclose(FID); clear variables; return;
    end
end

%% Copying ChannelID to MetaTags for filespec 2.2 and 2.3 for compatibility with filespec 2.1
if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD')
    NSx.MetaTags.ChannelID = [NSx.ElectrodesInfo.ElectrodeID]';
end

%% Determining the number of channels to read and validating the input
if ~elecReading
    if ~exist('userRequestedChanRow', 'var')
        userRequestedChannels = NSx.MetaTags.ChannelID;
    else
        if any(userRequestedChanRow > ChannelCount)
            disp(['Channel file only contains ' num2str(ChannelCount) ' channels.']);
            fclose(FID); clear variables; return;
        else
            userRequestedChannels = NSx.MetaTags.ChannelID(userRequestedChanRow);
        end
    end
end
for idx = 1:length(userRequestedChannels)
    if ~any(ismember(NSx.MetaTags.ChannelID, userRequestedChannels(idx)))
        disp(['Electrode ' num2str(Mapfile.Channel2Electrode(userRequestedChannels(idx))) ' does not exist in this file.']);
        fclose(FID); clear variables; return;
    end
    userRequestedChanRow(idx) = find(NSx.MetaTags.ChannelID == userRequestedChannels(idx),1);
end

%% Adjusts StartPacket and EndPacket based on what time setting (sec, min,
%  hour, or packets) the user has indicated in the input argument.
if ~NSx.RawData.PausedFile
    if ~exist('EndPacket', 'var')
        EndPacket = NSx.MetaTags.DataPoints;
    end
    switch TimeScale
        case 'sec'
            StartPacket = double(StartPacket) * NSx.MetaTags.SamplingFreq;
            EndPacket = EndPacket * NSx.MetaTags.SamplingFreq;
        case 'min'
            StartPacket = StartPacket * NSx.MetaTags.SamplingFreq * 60;
            EndPacket = EndPacket * NSx.MetaTags.SamplingFreq * 60;
        case 'hour'
            StartPacket = StartPacket * NSx.MetaTags.SamplingFreq * 3600;
            EndPacket = EndPacket * NSx.MetaTags.SamplingFreq * 3600;
    end
	EndPacket = EndPacket - 1;
    %% Validate StartPacket and EndPacket to make sure they do not exceed the
    %  length of packets in the file. If EndPacket is over then the last packet
    %  will be set for EndPacket. If StartPacket is over then will exist with an
    %  error message.
    if EndPacket >= NSx.MetaTags.DataPoints
        if StartPacket >= NSx.MetaTags.DataPoints
            disp('The starting packet is greater than the total data duration.');
            disp('The file was not read.');
            fclose(FID);
            if nargout; varargout{1} = []; end
            return;
        end
        disp('The time interval specified is longer than the data duration.');
        disp('Last data point will be used instead.');
        disp('Press enter to continue...');
        pause;
        EndPacket = NSx.MetaTags.DataPoints - 1;
    end
    DataLength = EndPacket - StartPacket + 1;
end
% from now StartPacket and EndPacket are in terms of Samples and are zero-based
clear TimeScale

% Reading the data if flag 'read' is used
if strcmp(ReadData, 'read')

    % Determine what channels to read
    numChansToRead = double(length(min(userRequestedChanRow):max(userRequestedChanRow)));
    if NSx.RawData.PausedFile
        for dataIDX = 1:segmentCount
            fseek(FID, f.BOData(dataIDX), 'bof');
            % Skip the file to the first channel to read
            fseek(FID, (find(NSx.MetaTags.ChannelID == min(userRequestedChannels))-1) * 2, 'cof');        
            NSx.Data{dataIDX} = fread(FID, [numChansToRead NSx.MetaTags.DataPoints(dataIDX)], [num2str(numChansToRead) precisionType], double((ChannelCount-numChansToRead)*2 + ChannelCount*(skipFactor-1)*2));
        end    
    else
        fseek(FID, f.BOData, 'bof');
        % Skip the file to the beginning of the time requsted, if not 0
        fseek(FID, StartPacket * 2 * ChannelCount, 'cof');
        
        % Skip the file to the first channel to read
        fseek(FID, (find(NSx.MetaTags.ChannelID == min(userRequestedChannels))-1) * 2, 'cof');        
        NSx.Data = fread(FID, [numChansToRead DataLength], [num2str(numChansToRead) precisionType], double((ChannelCount-numChansToRead)*2 + ChannelCount*(skipFactor-1)*2));
    end
end

%% Fixing a bug in 6.03 TOC where an extra 0-length packet is introduced
if NSx.RawData.PausedFile && strcmp(ReadData, 'read')
    if isempty(NSx.Data{1})
        NSx.Data = cell2mat(NSx.Data(2));
    end
end

%% Fixing a bug in 6.03 where data packets with 0 lengh may be added
if any(NSx.MetaTags.DataPoints == 0) && strcmp(ReadData, 'read')
    segmentsThatAreZero = find(NSx.MetaTags.DataPoints == 0);
    NSx.MetaTags.DataPoints(segmentsThatAreZero) = [];
    NSx.MetaTags.Timestamp(segmentsThatAreZero) = [];
    NSx.Data(segmentsThatAreZero) = [];
end

%% Removing extra channels that were read, but weren't supposed to be read
channelThatWereRead = min(userRequestedChanRow):max(userRequestedChanRow);
if ~isempty(setdiff(channelThatWereRead,userRequestedChanRow))
	deleteChannels = setdiff(channelThatWereRead, userRequestedChanRow) - min(userRequestedChanRow) + 1;
    if NSx.RawData.PausedFile
        for segIDX = 1:size(NSx.Data,2)
            NSx.Data{segIDX}(deleteChannels,:) = [];
        end
    else
        NSx.Data(deleteChannels,:) = [];
    end
end

%% Adjusting the ChannelID variable to match the read electrodes
channelIDToDelete = setdiff(1:ChannelCount, userRequestedChanRow);
NSx.MetaTags.ChannelID(channelIDToDelete) = [];

%% Calculating the DataPoints in seconds and adding it to MetaData
NSx.MetaTags.DataDurationSec = double(NSx.MetaTags.DataPoints)/NSx.MetaTags.SamplingFreq;

%% Displaying a report of basic file information and the Basic Header.
if strcmp(Report, 'report')
    disp( '*** FILE INFO **************************');
    disp(['File Path          = '  path]);
    disp(['File Name          = '  fname   ]);
	disp(['File Version       = '  NSx.MetaTags.FileSpec   ]);
    disp(['Duration (seconds) = '  num2str(NSx.MetaTags.DataDurationSec)]);
    disp(['Total Data Points  = '  num2str(NSx.MetaTags.DataPoints)]);
    disp(' ');
    disp( '*** BASIC HEADER ***********************');
    disp(['File Type ID       = '          NSx.MetaTags.FileTypeID      ]);
    disp(['Sample Frequency   = '  num2str(double(NSx.MetaTags.SamplingFreq))         ]);
    disp(['Electrodes Read    = '  num2str(double(NSx.MetaTags.ChannelCount))   ]);
    disp(['Data Point Read    = '  num2str(size(NSx.Data,2))]);
end

%% If user does not specify an output argument it will automatically create
%  a structure.
outputName = ['NS' fext(4)];
if (nargout == 0),
    assignin('caller', outputName, NSx);
else
    varargout{1} = NSx;
end

if strcmp(Report, 'report')
	disp(['The load time for ' outputName ' file was ' num2str(toc, '%0.1f') ' seconds.']);
end
fclose(FID);

end