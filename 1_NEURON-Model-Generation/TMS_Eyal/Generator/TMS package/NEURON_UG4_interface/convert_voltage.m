%% Read files
locs = load(['..' filesep 'results' filesep 'locs' filesep 'locs_all_seg.txt']);
voltage = load(['..' filesep 'results' filesep 'voltage_trace.dat'])';
%% Parse
output_folder = ['voltage_data_calcium' filesep];
for time_step = 1:size(voltage,2)
    output = [locs(:,1:3) voltage(:,time_step)];
    writematrix(output,[output_folder 'vm_' sprintf('%07d', time_step-1) '.dat'], 'Delimiter', 'space');
end
