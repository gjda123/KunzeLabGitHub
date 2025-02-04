%% %%%%%%%%%%%%%%%%%%%%   Read Electrode Data    %%%%%%%%%%%%%%%%%%%%%%% %%
% Written and maintained by Connor Beck
%                  contact: Connorbeck1997@gmail.com
% Updated June 2023
%%%%%%%%%%%%%%%%%%%%%%%%%      OVERVIEW      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Imports '.h5' files (exported from multichannel systems via
%   DataManager). 
%   Prompts the user for the .h5 file containing the MEA data if not provided
%   with a filepath to appropriate data. Input 'ask' as the argument.
%   Returns a cell array of the electrode data, an associated time vector,
%   and the map from labels to their indices in the cell array.
%   
%   IMPORTANT: All input files must be saved directly from ImageJ with the
%   mean intensity after multimeasure - see ImageJ analysis for
%   clarification.
%
%   Recommended Call Format:
%   [Parameters,Data]=load_MEA();
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%      INPUTS       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   REQUIRED INPUT ARGUMENT
%   The input argument is either 
%   A)  Nothing - you will be prompted to select a file to analyze.
%   B)  Parameters, containing attribute: Filename 
%       E.g. Parameters.Filename='Filepath/Filename.h5'
%       A char array identifying a '.h5' text file containing the data to be
%       analyzed.
%
%   The '.h5' file must follow Multichannel Systems formatting which can be
%   saved via Multichannel DataManager software with an export as '.hdf5'.
%
%   Invoke load_MEA() with no arguments to search for a Dataset on your computer.
%%%%%%%%%%%%%%%%%%%%%%%%%%      OUTPUTS       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Invoking load_MEA() returns:

%               Name             | Type          | Description 
%   Parameters
%               Filename            | String        | Name and path of imported
%                                                     file
%               n_electrodes        | double        | Number of electrodes
%               t_max               | double        | Maximum timesteps
%                                                     (samples)
%               ElectrodeLabel      | Cell array    | Map of MATLAB Index to
%                                                     Multichannel Index
%               SamplingFrequency   | double        | Frequency of
%                                                     recording
%   Data
%               t                   | double array  | time sequence array
%                                                     in seconds
%               Electrodes          | Object-
%                                   | double array  | Voltage measurements
%                                                     of each electrode

%   Everything computed within the function is passed through the data
%   structure.
%
%%  NAME VALUE PAIR ARGUMENTS
%   The base load_MEA() function can be used with these input
%   arguements to import the data into MATLAB.

%   The optional argument name/value pairs are:
%
%    NAME                   VALUE
%
%   no input                no input   
%                           Opens user interface to select a '.h5' file to
%                           analyze.
%
%   ''Parameters'           String
%                           Select an '.h5' file to analyze by inputting a
%                           filename into Parameters.Filename 
%                           - the file must either be:
%                               (1) located within the folder of the code
%                               or
%                               (2) A full path must be included in the
%                               text.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%      CODE       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Parameters,Data] = load_MEA(varargin)
    %INPUT ARGUEMENTS
    switch nargin
        case 0
            [mea_file, mea_path] = uigetfile('*.h5','Select MEA data');
            Parameters.Filename = strcat(mea_path, mea_file);
        case 1
            Parameters = varargin{1};
    end
    %In case the mea_filename is empty,
    if isempty(Parameters.Filename)
        [mea_file, mea_path] = uigetfile('*.h5','Select MEA data');
        Parameters.Filename = strcat(mea_path, mea_file);
    end

    % Set paths for specific data
    %   ChannelPath contains the location of all of the channel data (e.g.
    %   the voltages over time across all 60 channels).
    ChannelPath = '/Data/Recording_0/AnalogStream/Stream_0/ChannelData';

    %InfoChannel contains the information about which indexed electrode
    %corresponds to which Multichannel ID (e.g. index 1 = electrode '47')
    InfoChannelPath = '/Data/Recording_0/AnalogStream/Stream_0/InfoChannel';
    
    % Load in electrode parameters for information
    InfoChannel = h5read(Parameters.Filename,InfoChannelPath);

    keySet=InfoChannel.Label;
    valueSet=1:length(keySet);
    Parameters.L=containers.Map(keySet,valueSet);

    % Load MEA data
    %imports the Channel Data as an Array where each row contains the
    %voltage at a specified time, and each column is an electrode
    H = waitbar(0,'Reading electrode data...');   %Opens a waitbar for the impatient 
    ChannelData = h5read(Parameters.Filename,ChannelPath); %reads the data

    % Identify how many electrodes are in the Data
    Parameters.n_electrodes=size(InfoChannel.ChannelID,1);
    Parameters.t_max=size(ChannelData,1);

    % Create MEA index Map
    %Locate electrode labels to map index electrode ID to multichannel
    %label ID (e.g. index 1 = electrode '47')
    Parameters.ElectrodeLabel=InfoChannel.Label;
    
    % Create time vector for MEA data and extract the data of the EOI
    fs=InfoChannel.Tick;
    if all(fs==fs(1)) %Ensure each electrode has the same sampling frequency
        %Tick gives time between samples in microseconds, so to measure 
        % sampling frequency, scale to seconds and take the inverse.
        Parameters.samplingFrequency=1/(double(fs(1))*(1E-6)); 
    else
        error('Electrodes sample at different frequencies')
    end
    t = 0:1/Parameters.samplingFrequency:(Parameters.t_max-1)/Parameters.samplingFrequency;
    Data.t = t';

    % Obtain conversion factor of data (depends on compression of hdf), then
    % apply it and convert the electrode data from picovolts to microvolts
    ADzero = double(InfoChannel.ADZero);  %ADC-Step that represents the 0-point of the measuring range of the ADC
    ConversionFactor = double(InfoChannel.ConversionFactor); %Conversion factor for the mapping ADC-Step to Measured Value 
    Exponent= double(InfoChannel.Exponent); %Exponent n for the 1En responsible for which the channel values magnitue is measured (e.g. kilo,micro,milli)
    for i = 1:Parameters.n_electrodes
        Electrode=(double(ChannelData(:,i))-ADzero(i))*ConversionFactor(i)*10^Exponent(i);
        Data.Electrodes(i).RawElectrode = Electrode; %Insert Electrode into Data Array
        waitbar(i/60)
    end
    delete(H)
end