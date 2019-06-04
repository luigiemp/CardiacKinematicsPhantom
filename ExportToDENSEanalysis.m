% % Luigi Perotti and Ilya Verzhbinsky - 2019 %
% pulled from https://github.com/luigiemp/CardiacKinematicsPhantom

% This code takes the phantom image from DENSEphantom_Driver.m and converts
% it into a .dns file that can be read into the DENSEanalysis software,
% where displacement data can be re-processed from the phantom DENSE
% images. 

% The DENSEanalysis software can be found at:
% https://github.com/denseanalysis/denseanalysis


close all
clc
clear
%% Inputs
slice = 12; %Location of imaging slice in mm

%Where all exported image matfiles are stored
MatPath =  ''; 

%Matfile with entire Phantom workspace from DENSEphantom_Driver.m
PhantomData = load(''); 

OutputPath = PhantomData.OutputPath;                
OutputName = ''; % Output filename for DENSEanalysis .dns file

xres = PhantomData.NvoxelX;
yres = PhantomData.NvoxelY;
nFrames = PhantomData.TimeSteps;
PixelSpacing    = [PhantomData.hx PhantomData.hy]; 
k_e            = PhantomData.ke_x;  

SNR = [PhantomData.SNR, Inf];
reps = PhantomData.reps;
%% Load phase data (in cycles) from Phantom

direc = dir(MatPath);
dirIndex = [direc.isdir];  % Find the index for directories
dirList = {direc(dirIndex).name}';

del = zeros(1,length(dirList));
for j = 1:length(dirList)
    tmp = dirList{j};
    if strcmp(tmp(1),'.')
        del(j) = 1;
    end
end
del = find(del);
dirList(del) = [];
dirList = sort_nat(dirList);

nSNR = length(dirList);

Xwrap = zeros(xres,yres,nFrames,nSNR,reps);
Ywrap = zeros(xres,yres,nFrames,nSNR,reps);
Zwrap = zeros(xres,yres,nFrames,nSNR,reps);
Mag = zeros(xres,yres,nFrames,reps);

count = 1;
for s = 1:nSNR
    SNRfile = dirList{s};
    direc = dir(fullfile(MatPath,SNRfile));
    dirIndex = [direc.isdir];  % Find the index for directories
    repList = {direc(dirIndex).name}';

    del = zeros(1,length(repList));
    for j = 1:length(repList)
        tmp = repList{j};
        if strcmp(tmp(1),'.')
            del(j) = 1;
        end
    end
    del = find(del);
    repList(del) = [];
    repList = sort_nat(repList);
    
    for r = 1:reps
        
        repfile = repList{r};
        direc = dir(fullfile(MatPath,SNRfile,repfile));
        dirIndex = [direc.isdir];  % Find the index for directories
        fileList = {direc(~dirIndex).name}';
        
        del = zeros(1,length(fileList));
        for j = 1:length(fileList)
            tmp = fileList{j};
            if strcmp(tmp(1),'.')
                del(j) = 1;
            end
        end
        del = find(del);
        fileList(del) = [];
        fileList = sort_nat(fileList);
        
        
        
        for j = 1:length(fileList)
            tmp = fileList{j};
            dat = load(fullfile(MatPath,SNRfile,repfile,tmp)); 

            if strcmp(tmp(1),'m')
                t = str2double(tmp(5:6));
                Mag(:,:,t,s,r) = dat.magnitude_wN; 
            elseif strcmp(tmp(6),'X')
                t = str2double(tmp(8:9));
                Ywrap(:,:,t,s,r) = dat.phaseX_wN; %%% NOTE - x and y need to be swapped when loaded
            elseif strcmp(tmp(6),'Y')
                t = str2double(tmp(8:9));
                Xwrap(:,:,t,s,r) = dat.phaseY_wN; %%% NOTE - x and y need to be swapped when loaded
            elseif strcmp(tmp(6),'Z')
                t = str2double(tmp(8:9));
                Zwrap(:,:,t,s,r) = dat.phaseZ_wN;
            else
                error('Invalid File Name')
            end
        end
        clc
        fprintf(['loading image files ... ',num2str(count/(reps*nSNR) * 100), '%%\n'])
        count = count + 1;
    end

end

    
% convert cycle data to uint16
Mag = uint16(Mag .* 2^12);
Xwrap = uint16((Xwrap + 0.5) .* 2^12);
Ywrap = uint16((Ywrap + 0.5) .* 2^12); 
Zwrap = uint16((Zwrap + 0.5) .* 2^12);

   

%% load in template .dns file 
for SNRid = 1:nSNR
    template = fullfile('template.dns'); %% Template DENSEanalysis worspace.

    tmp = load(template, '-mat');

    names = fieldnames(tmp.seq);

    seq = tmp.seq;

    count = 0;
    for s = SNRid  
        for r = 1:reps
            for n = 1:4
                num = count*4 + n;
                for f = 1:length(names)
                    seq(num).(names{f}) = tmp.seq(n).(names{f});
                end
            end
            count = count + 1;
        end
    end


    pos_orig = seq(1).ImagePositionPatient;

    empty = repmat({''}, [nFrames, 1]);

    count = 0;
    for s = SNRid
        for r = 1:reps
            for n = 1:4
                num = count*4 + n;

                seq(num).Width       = size(Mag,2);
                seq(num).Height      = size(Mag,1);
                seq(num).Filename    = empty;
                seq(num).FileModDate = empty;
                seq(num).FileSize    = empty;
                seq(num).MediaStorageSOPInstanceUID = empty;
                seq(num).InstanceCreationTime = empty;

                seq(num).SOPInstanceUID = empty;
                seq(num).AcquisitionTime = empty;
                seq(num).ContentTime = empty;
                seq(num).TriggerTime = num2cell((1:nFrames).');
                seq(num).NominalInterval = num2cell(1000 * ones(nFrames,1));

                seq(num).InstanceNumber = num2cell((1:nFrames).');

                seq(num).LargestPixelValue = 2^12;
                seq(num).WindowCenter = 2048;
                seq(num).WindowWidth = 2048;

                seq(num).CardiacNumberOfImages = nFrames;
                seq(num).NumberInSequence = nFrames;

                %shift slice positition
                ori = seq(n).ImageOrientationPatient;
                a = ori(1:3); b = ori(4:6);

                v = cross(a,b);
                v = v * slice;

                seq(num).ImagePositionPatient = pos_orig + v;
            end
            count = count + 1;
        end
    end

    %% Fill in Mag, x, y and z metadata

    count = 0;
    for s = SNRid
        for r = 1:reps
            num = count*4;

            % Mag data
            seq(num + 1).DENSEid = 'mag.overall'; 
            seq(num + 1).DENSEindex = num2cell((1:nFrames).');
            seq(num + 1).PixelSpacing = PixelSpacing;
            seq(num + 1).DENSEdata = struct('Number',     nFrames,...
                                            'Partition',  [1 1],...
                                            'Scale',      [],...
                                            'EncFreq',    [],...
                                            'SwapFlag',   0,...
                                            'NegFlag',    [0 0 0]);

            imagecomments = empty;
            for i = 1:nFrames
                fmt = ['DENSE overall mag - Rep:0/1 Slc:0/1 Par:0/1 Phs:%d/%d ',...
                        'RCswap:0 RCSflip:0/0/0'];
                imagecomments{i} = sprintf(fmt, i-1, nFrames);
            end
            seq(num + 1).ImageComments = imagecomments;


            % x phase
            seq(num + 2).DENSEid = 'pha.x';
            seq(num + 2).DENSEindex = num2cell((1:nFrames).');
            seq(num + 2).PixelSpacing = PixelSpacing;
            seq(num + 2).DENSEdata = struct('Number',     nFrames,...
                                            'Partition',  [1 1],...
                                             'Scale',      1,...
                                             'EncFreq',    k_e,...
                                             'SwapFlag',   0,...
                                             'NegFlag',    [0 0 0]);

            imagecomments = empty;
            for i = 1:nFrames
                fmt = ['DENSE x-enc pha - Scale:1.000000 EncFreq:%0.2f Rep:0/1 ',...
                       'Slc:0/1 Par:0/1 Phs:%d/%d RCswap:0 RCSflip:0/0/0'];
                imagecomments{i} = sprintf(fmt, k_e, i-1, nFrames);
            end
            seq(num + 2).ImageComments = imagecomments;

            % y phase
            seq(num + 3).DENSEid = 'pha.y';
            seq(num + 3).DENSEindex = num2cell((1:nFrames).');
            seq(num + 3).PixelSpacing = PixelSpacing;
            seq(num + 3).DENSEdata = struct('Number',     nFrames,...
                                            'Partition',  [1 1],...
                                            'Scale',      1,...
                                            'EncFreq',    k_e,...
                                            'SwapFlag',   0,...
                                            'NegFlag',    [0 0 0]);

            imagecomments = empty;
            for i = 1:nFrames
                fmt = ['DENSE y-enc pha - Scale:1.000000 EncFreq:%0.2f Rep:0/1 ',...
                       'Slc:0/1 Par:0/1 Phs:%d/%d RCswap:0 RCSflip:0/0/0'];
                imagecomments{i} = sprintf(fmt, k_e, i-1, nFrames);
            end
            seq(num + 3).ImageComments = imagecomments;

            % z phase
            seq(num + 4).DENSEid = 'pha.z';
            seq(num + 4).DENSEindex = num2cell((1:nFrames).');
            seq(num + 4).PixelSpacing = PixelSpacing;
            seq(num + 4).DENSEdata = struct('Number',     nFrames,...
                                            'Partition',  [1 1],...
                                            'Scale',      1,...
                                            'EncFreq',    k_e,...
                                            'SwapFlag',   0,...
                                            'NegFlag',    [0 0 0]);

            imagecomments = empty;
            for i = 1:nFrames
                fmt = ['DENSE z-enc pha - Scale:1.000000 EncFreq:%0.2f Rep:0/1 ',...
                       'Slc:0/1 Par:0/1 Phs:%d/%d RCswap:0 RCSflip:0/0/0'];
                imagecomments{i} = sprintf(fmt, k_e, i-1, nFrames);
            end
            seq(num + 4).ImageComments = imagecomments;

            dns.img{num + 1} = Mag(:,:,:,s,r); 
            dns.img{num + 2} = Xwrap(:,:,:,s,r); 
            dns.img{num + 3} = Ywrap(:,:,:,s,r);
            dns.img{num + 4} = Zwrap(:,:,:,s,r);

            dnsname = ['SNR_',num2str(SNR(s)),'_Rep_',num2str(r)];

            dns.dns(count+1) = struct('Name',      dnsname,...
                             'UID',          dicomuid,...
                             'Type',         'xyz',...
                             'MagIndex',     [num+1 num+1 num+1],...
                             'PhaIndex',     [num+2 num+3 num+4],...
                             'Number',       nFrames,...
                             'PixelSpacing', PixelSpacing,...
                             'Scale',        [1 1 1],...
                             'EncFreq',      [k_e k_e k_e],...
                             'SwapFlag',     false,...
                             'NegFlag',      [false,false,false]);

         count = count + 1;
        end




    end
    

    dns.roi = struct([]);
    dns.img = dns.img';
    dns.dns = dns.dns';

    dns.seq = seq;
    %% Load roi

    DNS = dns;


    save(fullfile(OutputPath,[OutputName,'_SNR_',num2str(SNR(SNRid)),'.dns']), '-struct', 'DNS'); %% Save .dns file

end 

fprintf('.dns file exported\n')

