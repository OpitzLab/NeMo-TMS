function visualize_adjacent(neuron_data,calcium_data,results_folder)
%% Create folders
visualize_neuron(neuron_data,results_folder)
visualize_calcium(calcium_data,results_folder);
%% Check output
calcium_figs = dir([folder_out filesep 'gmsh' filesep 'png_calcium' filesep '*.png']);
neuron_figs = dir([folder_out filesep 'gmsh' filesep 'png_neuron' filesep '*.png']);
if (isempty(calcium_figs) || isempty(neuron_figs))
    error('No PNG file found!');
end
if ~isequal(length(calcium_figs),length(neuron_figs))
    error('The number of samples should be equal for NEURON and calcium simulations.');
end
%% Create folder
if(~exist([folder_out filesep 'video_adjacent'],'dir'))
    mkdir([folder_out filesep 'video_adjacent']);
end
%% Save adjacent figures as a movie
video = VideoWriter([folder_out filesep 'video_adjacent' filesep 'video'],'MPEG-4');
open(video);

for ii=1:length(calcium_figs)
  img_neuron = imread([folder_out filesep 'gmsh' filesep 'png_neuron' filesep 'view_' int2str(ii) '.png']);
  img_calcium = imread([folder_out filesep 'gmsh' filesep 'png_calcium' filesep 'view_' int2str(ii) '.png']);
  combined_img = cat(1,img_neuron, img_calcium);
  writeVideo(video,combined_img); %write the image to file
end
close(video); %close the file
end