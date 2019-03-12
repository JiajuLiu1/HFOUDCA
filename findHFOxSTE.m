function [m_HFOEvents, intervals,st_HFOSettings] = findHFOxSTE(EEG,chanid)
%Low freq, high freq, and min number of oscilations to be relevant
st_HFOSettings.s_FreqIni= 80;
st_HFOSettings.s_FreqEnd= floor(EEG.srate/2-1);
st_HFOSettings.s_NumOscMin= 6;
%Generally should not change
st_HFOSettings.s_EpochTime = 600;
st_HFOSettings.s_RMSThres = 4;
st_HFOSettings.s_RMSWindow = 3;
st_HFOSettings.s_BPThresh = 3;
st_HFOSettings.s_MinTime= 10;
st_HFOSettings.s_MinWind = 6;
%% Variable declarations
    %do 1 chan at a time
    ist = chanid;
s_SampleFrec = EEG.srate;
pv_Signal       = EEG.data(ist,:)';
v_Freqs         = [st_HFOSettings.s_FreqIni st_HFOSettings.s_FreqEnd];% Filter freqs
s_Window        = st_HFOSettings.s_RMSWindow * 1e-3;           % RMS window time (ms)
s_RMSThresSD    = st_HFOSettings.s_RMSThres;                   % Threshold for RMS in standard deviation
s_MinWind       = st_HFOSettings.s_MinWind * 1e-3;             % Min window time for an HFO (ms)
s_MinTime       = st_HFOSettings.s_MinTime * 1e-3;             % Min Distance time Betwen two HFO candidates
s_NumOscMin     = st_HFOSettings.s_NumOscMin;                  % Minimum oscillations per interval
s_BPThresh      = st_HFOSettings.s_BPThresh;                   % Threshold for finding peaks
s_EpochLength   = st_HFOSettings.s_EpochTime;                  % Cycle Time

%% Preprocessing Filter

s_Filter        = f_GetIIRFilter(s_SampleFrec,v_Freqs);
v_SigFilt       = f_IIRBiFilter(pv_Signal,s_Filter);
clear s_Filter

%% RMS Calculus

s_Window        = round(s_Window * s_SampleFrec);
if mod(s_Window, 2) == 0
    s_Window = s_Window + 1;
end
v_Temp                      = v_SigFilt.^2;
v_Temp                      = filter(ones(1,s_Window),1,v_Temp)./s_Window;
v_RMS                       = zeros(numel(v_Temp), 1);
v_RMS(1:end - ceil(s_Window / 2) + 1) = v_Temp(ceil(s_Window / 2):end);
v_RMS                       = sqrt(v_RMS);

clear v_Temp

%% Thresholding

    s_MinWind       = round(s_MinWind * s_SampleFrec);
    s_MinTime       = round(s_MinTime * s_SampleFrec);
    s_EpochLength   = round(s_EpochLength * s_SampleFrec);
    v_EpochTemp     = (1:s_EpochLength:length(pv_Signal))';
    
    if v_EpochTemp(end) < length(pv_Signal)
        v_EpochTemp(end+1)  = length(pv_Signal);
    end
    
    m_EpochLims     = [v_EpochTemp(1:end-1) v_EpochTemp(2:end)-1];
    m_HFOEvents = [];
    
    clear v_EpochTemp s_EpochLength
    
    for ii = 1:size(m_EpochLims,1)
        
        v_Window        = zeros(numel(v_RMS),1);
        v_Window(m_EpochLims(ii,1):m_EpochLims(ii,2)) = 1;
        
        v_RMSEpoch      = v_RMS.*v_Window;
        v_RMSInterval   = v_RMS(m_EpochLims(ii,1):m_EpochLims(ii,2));
        v_EpochFilt     = v_SigFilt(m_EpochLims(ii,1):m_EpochLims(ii,2));

        v_RMSThres      = v_RMSEpoch > ...
                            (mean(v_RMSInterval)+ ...
                                s_RMSThresSD*std(v_RMSInterval));
                            
        v_WindThres     = [0;v_RMSThres;0];
        v_WindJumps     = diff(v_WindThres);
        v_WindJumUp     = find(v_WindJumps==1);
        v_WindJumDown   = find(v_WindJumps==-1);
        v_WinDist       = v_WindJumDown - v_WindJumUp;
        v_WindIni       = v_WindJumUp(v_WinDist > s_MinWind);  
        v_WindEnd       = v_WindJumDown(v_WinDist > s_MinWind)-1;
        
        clear v_WindThres v_WindJumps v_WindJumUp v_WindJumDown 
        clear v_DistSelect v_WindSelect   
    
        while 1
            v_NextIni   = v_WindIni(2:end);
            v_LastEnd   = v_WindEnd(1:end-1);
            v_WinIdx	= (v_NextIni - v_LastEnd) < s_MinTime;
            if sum(v_WinIdx)==0
                break
            end
            v_NewEnd    = v_WindEnd(2:end);
            
            v_LastEnd(v_WinIdx) = v_NewEnd(v_WinIdx);
            v_WindEnd(1:end-1)  = v_LastEnd;
            
            v_Idx       = diff([0;v_WindEnd])~=0;
            v_WindIni   = v_WindIni(v_Idx);
            v_WindEnd   = v_WindEnd(v_Idx);        
        end
        
        m_WindIntervals = [v_WindIni v_WindEnd];
        
        clear v_WindSelect v_WindIni v_WindEnd v_WinDist
   
        s_Count             = 1;
        m_WindSelect        = zeros(size(m_WindIntervals));

        s_Threshold         = mean(abs(v_EpochFilt)) + ...
                                        s_BPThresh.*std(abs(v_EpochFilt));
        s_TotalWindInterv	= size(m_WindIntervals,1);


        for jj=1:s_TotalWindInterv

            v_Temp          = abs(v_SigFilt(m_WindIntervals(jj,1):...
                                                    m_WindIntervals(jj,2)));
                                                
            if numel(v_Temp) < 3
                continue
            end
            
            s_NumPeaks      = findpeaks(v_Temp,'minpeakheight',s_Threshold);
            clear v_Temp

            if isempty(s_NumPeaks) || length(s_NumPeaks) < s_NumOscMin
                continue;
            end

            m_WindSelect(s_Count,:) = [m_WindIntervals(jj,1)...
                                                    m_WindIntervals(jj,2)];
            s_Count                 = s_Count + 1;

        end
        
        if any(m_WindSelect(:))
            m_HFOEvents     = vertcat(m_HFOEvents,...
                                    m_WindSelect(1:s_Count-1,:)); %#ok<AGROW>
        end

    end

    %Get filtered intervals
    sizeint=50;
    intervals = zeros(size(m_HFOEvents,1),sizeint);
    for i = 1:size(m_HFOEvents)
        avg = floor((m_HFOEvents(i,1)+m_HFOEvents(i,2))/2);
        try
            v_SigFilt((avg-floor(sizeint/2)-1):(avg+floor(sizeint/2)));
            intervals(i,:) = v_SigFilt((avg-floor(sizeint/2)-1):(avg+floor(sizeint/2)));
        catch
            try
            intervals(i,:) = v_SigFilt(m_HFOEvents(i,1):(m_HFOEvents(i,1)+sizeint-1));
            catch
                intervals(i,:) = v_SigFilt((m_HFOEvents(i,2)-(sizeint-1)):(m_HFOEvents(i,1)));
            end
        end
    end
end

