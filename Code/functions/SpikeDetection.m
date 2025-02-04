%% %%%%%%%%%%%%%%%%%%   Electrode Spike Detection    %%%%%%%%%%%%%%%%%%% %%
% Written and maintained by Connor Beck
%                  contact: Connorbeck1997@gmail.com
% Updated June 2023
%%%%%%%%%%%%%%%%%%%%%%%%%      OVERVIEW      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Takes filtered Electrode Data and detects events based on falling edge
%   and the refractory period: 
%   base - 7 standard deviations & 3 ms refractory period
%   
%
%   Recommended Call Format:
%   [Parameters,Data]=SpikeDetection(Parameters,Data);
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%      INPUTS       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   REQUIRED INPUT ARGUMENTS
%   Data & Parameters
%   
%   Data and Parameters must be output from the load_MEA() and
%   filterElectrodes() functions before being used here.
%
%   Parameters (can) include attributes:
%
%   Parameters.standard_deviation=standard deviations; 
%       where standard deviations can be a value of how many standard
%       deviations the voltage must be away from 0 for it to be considered
%       a spike. Base is 7.
%   Parameters.refractory period =  refractory period;
%       where refractory period is a time in ms restricting how long a signal
%       must go after an event before the next event can be detected.
%       Biological limits of neurons often require 3 ms delay.
%
%  OPTIONAL NAME VALUE PAIR ARGUMENTS
%   The base SpikeDetection function can be used with these optional input
%   arguements to fine tune the spike detection protocol

%   The optional argument name/value pairs are:
%
%    NAME                   VALUE
%
%   ''Window Detection''    integer (Default no Window Detection)
%                           -Recommended (20) in seconds
%                           Input a value in seconds to adjust the standard
%                           deviation window. The function slides along the
%                           timeframe creating a smaller window of data
%                           to compute the standard deviation for spike
%                           detection, centered around the timepoint.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%      OUTPUTS       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Invoking SpikeDetection() returns:
%
%               Name             | Type          | Description 
%   Parameters
%               All Previously Contained Values
%               
%               if standard_deviation not contained in parameters on input
%               standard_deviation  | double     | standard deviations the 
%                                                  voltage must be away 
%                                                  from 0 for it to be 
%                                                  considered a spike.
%
%               if refractory_period not contained in parameters on input
%               refractory_period   | double     | time in ms restricting 
%                                                  the time in-between spikes
%
%   Data
%               All Previously Contained Values
%               
%               Electrodes
%                   Spikes  | cell array    | time locations of all detected
%                                             events sorted into a cell for 
%                                             each electrode.
%
%               SpikeOutput | double array  | a Nx2 array where column 1
%                                             represents the times
%                                             associated with an event and
%                                             column 2 represents the
%                                             electrode number associated
%                                             with the spike.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%      CODE       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Parameters,Data] = SpikeDetection(Parameters,Data,varargin)
    
    %Initialize parameters
    Sliding_Window=false;

    if nargin>2 && ~mod(nargin,2)
        for i=1:2:length(varargin)
            if strcmp(varargin{i},'Window Detection')
                Sliding_Window=true;
                Window=varargin{i+1}*Parameters.samplingFrequency;
            end
        end
    end
    if ~isfield(Parameters,'standard_deviation') || isempty(Parameters.standard_deviation)
        Parameters.standard_deviation=7;
    end
    if ~isfield(Parameters,'refractory_period') || isempty(Parameters.refractory_period)
        Parameters.refractory_period=3;
    end
    if ~isfield(Parameters,'electrode_removal') || isempty(Parameters.electrode_removal)
        Parameters.electrode_removal={};
    end

    Data.SpikeOutput=[];
    H = waitbar(0,'Detecting Electrode Spikes...'); 

     
    for i=1:Parameters.n_electrodes
        if ~strcmp(Parameters.ElectrodeLabel{i},'Ref') && ~any(ismember(Parameters.ElectrodeLabel{i},Parameters.electrode_removal))
            waitbar(i/Parameters.n_electrodes)
            %Select an individual electrode for analysis
            Electrode = Data.Electrodes(i).FilteredElectrode;
            if ~Sliding_Window
%%%%%%%%%%%%%%%%%%%%%% Standard detection code %%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                % Calculate standard deviation of electrode
                electrodeDeviation = std(Electrode);
                Parameters.threshold(i) = (-1) * Parameters.standard_deviation * electrodeDeviation; %Calculate threshold voltage
                
                %find times where the electrode went below the threshold
                SpikeTimes=find(Electrode <=Parameters.threshold(i));
            else
%%%%%%%%%%%%%%%%%%% Sliding deviation detection code %%%%%%%%%%%%%%%%%%%%%%
                electrodeDeviation = zeros(Parameters.t_max,1);
                startDev=std(Electrode(1:Window));
                endDev=std(Electrode(end-Window:end));
                for j=1:Parameters.t_max
                    if i>round(Window/2) && i<Parameters.t_max-round(Window/2)
                        electrodeDeviation(j)=std(Electrode(i-round(Window/2):i+round(Window/2)));
                    elseif i<=round(Window/2)
                        electrodeDeviation(j)=startDev;
                    else
                        electrodeDeviation(j)=endDev;
                    end
                end
                Parameters.threshold(i)=mean(electrodeDeviation);
                %find times where the electrode went below the threshold
                SpikeTimes=find(Electrode <=((-1) * Parameters.standard_deviation) .*electrodeDeviation);
            end
%%%%%%%%%%%%%%%%%%%% Refractory period spike removal %%%%%%%%%%%%%%%%%%%%%%                
            % Comb through spike times, remove events that occured within
            % refractory period
            for j = 2:length(SpikeTimes)
                if (SpikeTimes(j)-SpikeTimes(j-1))<Parameters.refractory_period*1E3/Parameters.samplingFrequency
                    SpikeTimes(j)=[];
                end
            end
            %Insert Spike Times into the electrodes object as Spikes.
            Data.Electrodes(i).Spikes=SpikeTimes;
            ISI=diff(SpikeTimes);
            Data.Electrodes(i).ISI=ISI;
            
            %Burst Detection
            maxISI=100*(10^-3)*Parameters.samplingFrequency;
            minSpikes=10;
            
            burst=zeros(length(ISI),1);
            burstID=1;
            burstCount=0;
            burstISI=[];

            burstActive=false;
            Data.Electrodes(i).burst=[];
            for j=2:length(ISI)
                if ISI(j)<maxISI
                    burstActive=true;
                    burst(j-1)=1;
                    
                    burstISI=cat(1,burstISI,ISI(j));
                    burstCount=burstCount+1;

                    %Data.Electrodes(i).burst(burstID).burstISI;
                    Data.Electrodes(i).burst(burstID).burstCount=burstCount;

                    
                   
                elseif burstActive     
                    burstID=burstID+1;
                    burstCount=0;
                    burstActive=false;
                end
            end
            
            j=1;
            if ~isempty(Data.Electrodes(i).burst)
                nburst=length(Data.Electrodes(i).burst);
                while j<nburst
                    if Data.Electrodes(i).burst(j).burstCount<minSpikes
                        Data.Electrodes(i).burst(j)=[];
                        nburst=nburst-1;
                    else
                        j=j+1;
                    end
                    
                end
            else
                Data.Electrodes(i).burst=[];
            end



            if ~isempty(Data.Electrodes(i).burst)
                burst=[Data.Electrodes(i).burst.burstCount];
                Data.Electrodes(i).burstSpikePercent=sum(burst)./length(Data.Electrodes(i).Spikes);
                Data.Electrodes(i).meanSpikesPerBurst=mean(burst);

            else
                Data.Electrodes(i).burstSpikePercent=NaN;
                Data.Electrodes(i).meanSpikesPerBurst=NaN;
            end

            Data.Electrodes(i).burstCount=length(Data.Electrodes(i).burst);
            Data.Electrodes(i).meanISI=mean(Data.Electrodes(i).ISI)/Parameters.samplingFrequency;

            %Create an (x,y) [=] (Cell ID, Spike Time) array for plotting
            Spikes=cat(2,SpikeTimes,i*ones(length(SpikeTimes),1));
            if ~isempty(Spikes)
                Data.SpikeOutput=cat(1,Data.SpikeOutput,Spikes);
            end
        else
            %No spikes on reference electrode
            Data.Electrodes(i).Spikes=[];
        end
    end
    delete(H)
end