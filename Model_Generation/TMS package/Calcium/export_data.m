function status = export_data()
% Export the SWC file, and the voltage data in the format compatible with
% the calcium concentration simulations.
% Results data are placed in Results\Calcium
% Returns the status: missing input files = 0, incomplete voltage = 1, successful = 2
%% Create output folder
output_folder = ['..' filesep '..' filesep 'Results' filesep ...
    'Calcium' filesep 'Converted_Voltage_Traces'];
if(~exist(output_folder,'dir'))
    mkdir(output_folder); % create output directory
end
%% Generate SWC file
if exist(['..' filesep '..' filesep 'Results' filesep 'NEURON' filesep 'locs' filesep 'locs_all_seg.txt'])
    generate_SWC;
else
    status = 0;
    return
end
%% Generate voltage data in the correct format
if ~exist(['..' filesep '..' filesep 'Results' filesep 'NEURON' filesep 'voltage_trace.dat'])
    status = 0;
    return
elseif convert_voltage()
    status = 1;
    return
end
%% Copy the time vector
t_in = ['..' filesep '..' filesep 'Results' filesep 'NEURON' filesep 'tvec.dat'];
t_out = ['..' filesep '..' filesep 'Results' filesep 'Calcium' filesep];
copyfile(t_in, t_out, 'f');
status = 2;
end