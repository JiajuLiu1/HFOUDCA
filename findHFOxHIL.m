%   f_findHFOxSTF.m [As a part of HFO Detection Project]
%   Written by:
%   Miguel G. Navarrete Mejia
%   Electrical Engineering MS candidate
%   UNIVERSIDAD DE LOS ANDES
%   Colombia, 2012
%   mnavarretem@gmail.com

function [m_HFOEvents, intervals,st_Data] = findHFOxHIL(EEG,chanID)
st_Data.s_Type='HIL';
st_Data.s_FreqIni= 80;
st_Data.s_FreqEnd= floor(EEG.srate/2);
st_Data.s_EpochTime= EEG.xmax-EEG.xmin;
st_Data.s_SDThres= 5;
%WHY IS THIS SO DUMB CANT IT BE MIN OSC
st_Data.s_MinWind= 5;     
pv_Signal       = EEG.data(chanID,:);

v_Freqs         = [st_Data.s_FreqIni st_Data.s_FreqEnd];% Filter freqs
s_SDThres   	= st_Data.s_SDThres;                    % Threshold in standard deviation
s_MinWind       = st_Data.s_MinWind * 1e-3;             % Min window time for an HFO (ms)
s_EpochLength   = st_Data.s_EpochTime;                  % Cycle Time

%% Preprocessing Filter

s_Filter        = f_GetIIRFilter(EEG.srate,v_Freqs);
v_SigFilt       = f_IIRBiFilter(pv_Signal,s_Filter);
clear s_Filter


%% Hilbert transform Calculus

v_SigFilt       = abs(hilbert(v_SigFilt));

     
%% Thresholding


    s_EpochLength   = round(s_EpochLength * EEG.srate);
    v_EpochTemp     = (1:s_EpochLength:length(pv_Signal))';
    s_MinWind       = round(s_MinWind * EEG.srate);
    
    if v_EpochTemp(end) < length(pv_Signal)
        v_EpochTemp(end+1)  = length(pv_Signal);
    end
    
    m_EpochLims     = [v_EpochTemp(1:end-1) v_EpochTemp(2:end)-1];
   
    clear v_EpochTemp s_EpochLength
    
    m_HFOEvents = [];  

    
    for ii = 1:size(m_EpochLims,1)
        

        v_EpochFilt     = v_SigFilt(m_EpochLims(ii,1):m_EpochLims(ii,2));
        
        v_WinThres      = v_EpochFilt > ...
                            (mean(v_EpochFilt)+ s_SDThres*std(v_EpochFilt));

                            
        v_WindThres     = [0;v_WinThres;0];
        v_WindJumps     = diff(v_WindThres);
        v_WindJumUp     = find(v_WindJumps==1);
        v_WindJumDown   = find(v_WindJumps==-1)-1;        
        v_WinDist       = v_WindJumDown - v_WindJumUp;

        v_DistSelect    = (v_WinDist > s_MinWind);
        v_WindJumUp     = v_WindJumUp(v_DistSelect);  
        v_WindJumDown   = v_WindJumDown(v_DistSelect)-1;

    
        m_WindSelect	= [v_WindJumUp v_WindJumDown] + m_EpochLims(ii,1)-1;
        
        if any(m_WindSelect(:))
            m_HFOEvents     = vertcat(m_HFOEvents,m_WindSelect); %#ok<AGROW>
        end
        

    end

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


%% Last
% % m_WindSelect                = m_WindSelect(1:s_Count-1,:);
% v_WindInter         = zeros(numel(v_EpochFilt),1);
% for kk=1:length(m_WindIntervals)-1
%     v_WindInter(m_WindIntervals(kk,1):m_WindIntervals(kk,2))   = 1;
% end

% v_WindInterPks    = zeros(numel(pv_Signal),1);
% for kk=1:length(m_WindSelect)
% %     v_WindInterPks(m_WindSelect(kk,1):m_WindSelect(kk,2))   = 1;
% % end
% figure
% % plot(v_SigFilt,'b')
% hold on;
% plot(v_EpochFilt,'g')
% plot(v_WindInter.*v_EpochFilt,'r')

% plot(abs(v_RMS),'r')
% plot(0.5*v_RMSThres,'y')
% plot(v_WindInter,'m')
% plot(50*v_WindInterPks,'g')
% plot((mean(v_RMS)+ s_RMSThresSD*std(v_RMS)).*ones(1,numel(v_RMS)),'y');

% hold off

end