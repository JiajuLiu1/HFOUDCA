function [V,D,HFO_Count,filtData,rawData,set,info] = Full_Analysis(path,locPath,eegPath,chan,detect) 
%Analysis of EEG data from EDF file format. Inlcudes 60 line noise removal,
%surface laplacian, 3 (and growing) different detectors, and unsupervised
%component diffusion analysis.
%--------------------------------------------------------------------------
%Parameters:
%Out:
%V - Eigenvectors with 20 largest eigenvalues (listed in descending order)
%D - Respective 20 largest eigenvalues
%HFO_Count - HFO_Count(i) is the number of HFO detections in channel i.
%Therefore, sum(HFO_Count)=size(V,1)
%filtData - Filtered detected HFOs with baseline removed.
%rawData - Raw detected HFO data.
%set - Lists start and stop times of all detections
%info - Detector info to help identify files.
%In:
%path - File path of EDF file to be analyzed
%locPath - File path of channels to use. Equivalent of EEG.chanlocs
%eegPath - Path to folder that contains eeglab
%chan - index of channels to use ex:[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 19 29 30]
%detect - detector to be used. Currently allows "STE", "HIL", and "MNI"
%--------------------------------------------------------------------------
%NOTE: Currently, surface Laplacian function is from MikeXCohen at
%http://mikexcohen.com and detectors are modified from RippleLab.
%Otherwise, entirely coded by Jiaju Liu at jiajuliu1@gmail.com as part of
%EPIBIOS4Rx project.
chans = load(locPath);
chanlocs= chans{1};
HFO_Count =zeros(size(chanlocs,1),1);
filtData = [];
rawData = [];
set=[];
addpath(eegPath);
[~ , ~, ~, ~] = eeglab;
%HFO detection and preproccessing are done together
for i = 1%:size(fileLocs,1)
    EEG = pop_biosig(strcat(path,fileLocs),'channels',chan);
    EEG = pop_select( EEG,'time',[1 EEG.xmax-2] );
    if EEG.srate >1000
        EEG = pop_resample( EEG, 1000);
    end
    %Get rid of 60 line noise
    EEG = pop_eegfiltnew(EEG, 58,62,1650,1,[],1);
    EEG = pop_eegfiltnew(EEG, 118,122,1650,1,[],1);
    EEG = pop_eegfiltnew(EEG, 178,182,1650,1,[],1);
    EEG.chanlocs = chanlocs3_13_0062;
    EEG=pop_chanedit(EEG, 'lookup',[eegPath '/eeglab14_1_2b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp']);
    EEG.data = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    %Record filtered data from those intervals
    for j = 1:19
        switch detect
            case 'STE'
                [sets, intervals,info] = findHFOxSTE(EEG,j);
            case 'MNI'
                [sets, intervals,info] = findHFOxMNI(EEG,j);
            case 'HIL'
                [sets, intervals,info] = findHFOxHIL(EEG,j);
            otherwise
                error('Detector does not exist.')
        end
        %Bad programming practices :/
        avg=floor((sets(:,1)+sets(:,2))/2);
        for k = 1:size(sets,1)
            rawInts = (avg(k)-24):(avg(k)+25);
            rawData=[rawData; EEG.data(j,rawInts)];
        end
        filtData = [filtData; intervals];
        set = [sets; set];
        HFO_Count(j,i) = size(intervals,1);
    end
end
%UDCA
[V,D] = UDCA(filtData);

