%% %%%%%%%%%%%%%%%%%%%%   Filter Electrode Data    %%%%%%%%%%%%%%%%%%%%% %%
% Written and maintained by Connor Beck
%                  contact: Connorbeck1997@gmail.com
% Updated June 2023
%%%%%%%%%%%%%%%%%%%%%%%%%      OVERVIEW      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Takes Raw Electrode Data from the load_MEA() function and filters out
%   key noise elements. This can be used as either a low pass, high pass or
%   bandpass filter.
%   
%   Uses an input array with 2 doubles: [High Pass, Low Pass] where values
%   are in Frequency (Hz)
%
%
%   Recommended Call Format:
%   [Parameters,Data]=filterElectrodes(Parameters,Data);
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%      INPUTS       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   REQUIRED INPUT ARGUMENTS
%   Data & Parameters
%   
%   Data and Parameters must be output from the load_MEA function before 
%   running through filtering.
%
%   Parameters (can) include attribute:
%   Parameters.filter_frequencies=[high_pass,low_pass]; where high pass 
%   and low pass values are frequencies for filtering. 
%   
%   If you only wish to include 1 type of filtering, leave the other blank as:
%   Parameters.filter_frequencies[high_pass,];
%%%%%%%%%%%%%%%%%%%%%%%%%%      OUTPUTS       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Invoking filterElectrodes() returns:

%               Name             | Type          | Description 
%   Parameters
%               All Previously Contained Values
%               
%               if filter_frequencies not contained in parameters on input
%               filter_frequencies  | double array  | first value is high
%                                                     pass frequency and 
%                                                     second value is low 
%                                                     pass frequency

%   Data
%               All Previously Contained Values
%               
%           Electrodes
%               filteredElectrode  | double array  | Voltage values for
%                                                    electrodes after 
%                                                    filtering

%%%%%%%%%%%%%%%%%%%%%%%%%%%%      CODE       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Parameters,Data] = filterElectrodes(Parameters,Data)
    %Check if user input filter frequencies, if not, input standard values
    % High pass: 300 Hz, Low pass: 4000 Hz
    if ~isfield(Parameters,'filter_frequencies') || isempty(Parameters.filter_frequencies)
        high_pass=300;
        low_pass=4000;
        Parameters.filter_frequencies=[high_pass,low_pass];
    else
        %If user input frequencies:
        %Check if high_pass filter value input
        if ~isempty(Parameters.filter_frequencies(1))
            high_pass=Parameters.filter_frequencies(1);
        else %if not, set the high pass value to false
            high_pass=false;
        end
        %Check if low_pass filter value input
        if ~isempty(Parameters.filter_frequencies(2))
            low_pass=Parameters.filter_frequencies(2);
        else %if not set the low pass value to false
            low_pass=false;
        end
    end
    %Filter
    if ~islogical(high_pass) ||  ~islogical(low_pass) 
        %check if both high pass and low pass are false, if so, the signal
        %is not filtered
        H = waitbar(0,'Filtering electrode data...');
        for i = 1:Parameters.n_electrodes %individually loop through each electrode
            Electrode=Data.Electrodes(i).RawElectrode; %grab the signal
            if ~islogical(high_pass) %if there is a high pass value, filter
                Electrode = highpass(Electrode, high_pass, Parameters.samplingFrequency);
            else
                warning('Warning: No High Pass filter used: Signal may be noisy')
            end
            if ~islogical(low_pass) %if there is a low pass value, filter
                Electrode = lowpass(Electrode, low_pass, Parameters.samplingFrequency);
            end
            %Output the signal in the data electrodes array.
            Data.Electrodes(i).FilteredElectrode = Electrode;
            Data.Electrodes(i).Deviation = std(Electrode);
            waitbar(i/Parameters.n_electrodes) %waitbar for impatient people.
        end 
        delete(H)
    end
end