function visualize_adjacent(neuron_data,calcium_data,folder_out)
%% Check time steps
neuron_tvec = load([neuron_data filesep 'tvec.dat']);
calcium_tvec = load([calcium_data filesep 'timeSteps.txt'])*1000; % convert to ms
dt_neuron = median(diff(neuron_tvec));
dt_calcium = median(diff(calcium_tvec));
tol = 1e-8;
if(abs(dt_neuron-dt_calcium) > tol )
    error('The time steps should be equal for NEURON and calcium simulations.');
end
%% Create folders
visualize_neuron(neuron_data,folder_out)
visualize_calcium(calcium_data,folder_out);
%% Check output
calcium_figs = dir([folder_out filesep 'gmsh' filesep 'png_calcium' filesep '*.png']);
neuron_figs = dir([folder_out filesep 'gmsh' filesep 'png_neuron' filesep '*.png']);
if (isempty(calcium_figs) || isempty(neuron_figs))
    error('PNG files not found!');
end
if ~isequal(length(calcium_figs),length(neuron_figs))
    warning('The number of time points is not equal for NEURON and calcium simulations. Data will be visualized for the duration they both exist.');
end
%% Create folder
if(~exist([folder_out filesep 'video_adjacent'],'dir'))
    mkdir([folder_out filesep 'video_adjacent']);
end
%% Save adjacent figures as a movie
video = VideoWriter([folder_out filesep 'video_adjacent' filesep 'video'],'MPEG-4');
video.FrameRate = 20;
open(video);

nFigs = min([length(calcium_figs),length(neuron_figs)]);
for ii=1:nFigs
  img_neuron = imread([folder_out filesep 'gmsh' filesep 'png_neuron' filesep 'view_' int2str(ii) '.png']);
  img_calcium = imread([folder_out filesep 'gmsh' filesep 'png_calcium' filesep 'view_' int2str(ii) '.png']);
  combined_img = cat(1,img_neuron, img_calcium);
  writeVideo(video,combined_img); %write the image to file
end
close(video); %close the file
end