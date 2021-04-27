function incomplete = convert_voltage()
%% Read files
locs = load(['..' filesep '..' filesep 'Results' filesep 'NEURON' filesep 'locs' filesep 'locs_all_seg.txt']);
try
    voltage = load(['..' filesep '..' filesep 'Results' filesep 'NEURON' filesep 'voltage_trace.dat'])';
    tvec = load(['..' filesep '..' filesep 'Results' filesep 'NEURON' filesep 'tvec.dat']);
catch
    incomplete = 1;
    return
end
%% Check
if size(voltage,2) ~= length(tvec)
    incomplete = 1;
    return
end
%% Parse
output_folder = ['..' filesep '..' filesep 'Results' filesep ...
    'Calcium' filesep 'Converted_Voltage_Traces' filesep];
for time_step = 1:size(voltage,2)
    output = [locs(:,1:3) voltage(:,time_step)];
    writematrix(output,[output_folder 'vm_' sprintf('%07d', time_step-1) '.dat'], 'Delimiter', 'space');
end
incomplete = 0;
end