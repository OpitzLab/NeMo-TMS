function export_data()
% Export the SWC file, and the voltage data in the format compatible with
% the calcium concentration simulations.
% 'neuron_out.swc' file containing the neuron morphology is placed in the results folder
% Voltage data are placed in results\voltage_data_calcium
%% Create output folder
output_folder = ['voltage_data_calcium'];
if(~exist(output_folder,'dir'))
    mkdir(output_folder); % create output directory
end
%% Generate SWC file
generate_SWC;
%% Generate voltage data in the correct format
convert_voltage;
end