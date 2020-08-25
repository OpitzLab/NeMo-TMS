function visualize_neuron(folder_in,folder_out)
%% Load files
locs = load([folder_in filesep 'locs' filesep 'locs_all_seg.txt']);
v_trace = readmatrix([folder_in filesep 'voltage_trace.dat'])';
tvec = load([folder_in filesep 'tvec.dat']);
%% Extract pairs
nodes = locs(:,1:3);
line_list = zeros(size(locs,1)-1,2);
for ii = 2:size(locs,1)
    [~, temp] = ismember(locs(ii,4:6),locs(:,1:3),'rows');
    line_list(ii-1,:) = [ii temp];
end
%% Create folders
if(~exist([folder_out filesep 'gmsh'],'dir'))
    mkdir([folder_out filesep 'gmsh']);
end
if(~exist([folder_out filesep 'gmsh' filesep 'png_neuron'],'dir'))
    mkdir([folder_out filesep 'gmsh' filesep 'png_neuron']);
end
if(~exist([folder_out filesep 'video_neuron'],'dir'))
    mkdir([folder_out filesep 'video_neuron']);
end
%% Export gmsh views
m.points = nodes;
m.lines = line_list;

for ii = 1:size(v_trace,2)
    m.element_data{ii,1}.data = v_trace(2:end,ii);
    m.element_data{ii,1}.name = ['t = ' num2str(tvec(ii),'%.1f') ' ms'];
end

mesh_save_gmsh(m,[folder_out filesep 'gmsh' filesep 'spike_neuron.msh']);
%% Export png files from gmsh
fid = fopen('spike_movie_template_lines.geo','rt');
geo_file = fread(fid);
fclose(fid);
geo_file = char(geo_file.');
% replace SIMTYPE with proper number
geo_file_mod = strrep(geo_file, 'SIMTYPE', 'neuron');
% replace SPIKE_LEN with proper number
geo_file_mod = strrep(geo_file_mod, 'SPIKE_LEN', num2str(size(v_trace,2)));
% replace FILESEP with proper file separator
geo_file_mod = strrep(geo_file_mod, 'FILESEP', filesep);
% replace the value range with proper numbers
geo_file_mod = strrep(geo_file_mod, 'MINVAL', num2str(min(v_trace(:))));
geo_file_mod = strrep(geo_file_mod, 'MAXVAL', num2str(max(v_trace(:))));
fid = fopen([folder_out filesep 'gmsh' filesep 'spike_movie_neuron.geo'],'wt');
fwrite(fid,geo_file_mod);
fclose (fid);
system(['gmsh ' folder_out filesep 'gmsh' filesep 'spike_movie_neuron.geo']);
%% Save figures as a movie
video = VideoWriter([folder_out filesep 'video_neuron' filesep 'voltage'],'MPEG-4');
video.FrameRate = 20;
open(video);
for ii=1:size(v_trace,2)
  Im = imread([folder_out filesep 'gmsh' filesep 'png_neuron' filesep 'view_' int2str(ii) '.png']);
  writeVideo(video,Im); %write the image to file
end
close(video); %close the file
end
