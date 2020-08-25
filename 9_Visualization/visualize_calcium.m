function visualize_calcium(folder_in,folder_out)
%% Load files
data = readmatrix([folder_in filesep 'fullCalciumData.txt'])';
nodes = load([folder_in filesep 'outDom.txt']);
swc = load([folder_in filesep 'neuron_out.swc']);
swc(:,3:5) = swc(:,3:5)*1e6; % convert to um
tstep = load([folder_in filesep 'timeSteps.txt'])*1000; % convert to ms
%% Find the indecis for equispace time points
dt = median(diff(tstep));
tol = dt*1e-8;
temp = dt*(0:(length(tstep)-1));
tindex = find(ismembertol(tstep,temp,tol));
%% Extract pairs
line_list = zeros(size(nodes,1)-1,2);
for ii = 2:size(nodes,1)
    dist = sum(bsxfun(@minus, swc(:,3:5), nodes(ii,2:4)).^2,2);
    [~,node_idx] = min(dist);
    line_list(ii-1,:) = [node_idx swc(node_idx,7)];
end
%% Create folders
if(~exist([folder_out filesep 'gmsh'],'dir'))
    mkdir([folder_out filesep 'gmsh']);
end
if(~exist([folder_out filesep 'gmsh' filesep 'png_calcium'],'dir'))
    mkdir([folder_out filesep 'gmsh' filesep 'png_calcium']);
end
if(~exist([folder_out filesep 'video_calcium'],'dir'))
    mkdir([folder_out filesep 'video_calcium']);
end
%% Export gmsh views
m.points = nodes(:,2:4);
m.lines = line_list;

for ii = 1:length(tindex)
    m.element_data{ii,1}.data = data(2:end,tindex(ii));
    m.element_data{ii,1}.name = ['t = ' num2str(tstep(tindex(ii)),'%.1f') ' ms'];
end

mesh_save_gmsh(m,[folder_out filesep 'gmsh' filesep 'spike_calcium.msh']);
%% Export png files from gmsh
fid = fopen('spike_movie_template_lines.geo','rt');
geo_file = fread(fid);
fclose(fid);
geo_file = char(geo_file.');
% replace SIMTYPE with proper number
geo_file_mod = strrep(geo_file, 'SIMTYPE', 'calcium');
% replace SPIKE_LEN with proper number
geo_file_mod = strrep(geo_file_mod, 'SPIKE_LEN', num2str(size(data,2)));
% replace FILESEP with proper file separator
geo_file_mod = strrep(geo_file_mod, 'FILESEP', filesep);
% replace the value range with proper numbers
geo_file_mod = strrep(geo_file_mod, 'MINVAL', num2str(5e-8));%min(data(:))));
geo_file_mod = strrep(geo_file_mod, 'MAXVAL', num2str(1e-6));%max(data(:))));
fid = fopen([folder_out filesep 'gmsh' filesep 'spike_movie_calcium.geo'],'wt');
fwrite(fid,geo_file_mod);
fclose (fid);
system(['gmsh ' folder_out filesep 'gmsh' filesep 'spike_movie_calcium.geo']);
%% Save figures as a movie
video = VideoWriter([folder_out filesep 'video_calcium' filesep 'calcium'],'MPEG-4');
video.FrameRate = 20;
open(video);
for ii=1:length(tindex)
  Im = imread([folder_out filesep 'gmsh' filesep 'png_calcium' filesep 'view_' int2str(ii) '.png']);
  writeVideo(video,Im); %write the image to file
end
close(video); %close the file
end
