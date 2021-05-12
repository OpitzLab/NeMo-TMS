function couple_gui()
%% Load mesh and remove unnecessary attributes
[meshfile,meshpath] = uigetfile('*.msh', 'Select the mesh file output from SimNIBS simulation');
mesh = mesh_load_gmsh4([meshpath meshfile]);

% Only keep the vector E-field
for idx = 1:length(mesh.element_data)
    if strcmp(mesh.element_data{idx}.name, 'E')
        break
    end
end
if ~strcmp(mesh.element_data{idx}.name, 'E')
    message = (['The input mesh file from SimNIBS does not include any electric fields. ' ...
        'Make sure to select "vector E" in "Simulation Options" in SimNIBS.']);
    errordlg(message,'Error: Missing Electric Fields');
    return
end
mesh.element_data = mesh.element_data(idx);%{mesh.element_data{idx}};
%% Plot mesh
gm_surf = mesh_extract_regions(mesh,'elemtype','tri','region_idx',1002);
p = gm_surf.nodes;
t = gm_surf.triangles;
fig = figure('units','normalized','outerposition',[0 0 1 1]);
% tr = triangulation(t,p(:,1),p(:,2),p(:,3));
gm_normE = normd(gm_surf.element_data{1}.tridata);
trisurf(t,p(:,1),p(:,2),p(:,3),gm_normE);
colorbar;
cbar = colorbar;
cbar.Label.String = '|E| (V/m)';
title('FEM Model');
axis equal
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
%% Enable datatip
dcm_obj = datacursormode(fig);
set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
%% Receive neuron location
prompt = {'\bfEnter the coordinates of the center of neuron on the grey matter surface in mm:\rm (X Y Z with space in between)',...
    '\bf Enter the depth of the center of neuron from the grey matter surface in mm:'};
dlgtitle = 'Neuron location';
dims = [1 92];
opts.Interpreter = 'tex';
opts.WindowStyle = 'normal';
definput = {'','1'};
while 1
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    nrnloc = str2double(split(answer{1}));
    nrndpth = str2double(answer{2});
    if length(nrnloc) == 3 && sum(isnan(nrnloc)) == 0 && ~isnan(nrndpth)
        break
    end
    definput = answer;
    if length(nrnloc) ~= 3 || sum(isnan(nrnloc)) ~= 0
        definput{1} = 'Wrong format!';
    end
    if isnan(nrndpth)
        definput{2} = 'Wrong format!';
    end
end
nrnloc = nrnloc';

try
    close(fig);
catch
end
%% Load neuron and plot
nrnfile = 'locs_all_seg.txt';
nrnpath = [fullfile('..','..','Results','NEURON','locs') filesep];
try
    locs = load([nrnpath nrnfile]);
catch
    msgbox('Export NEURON model segments first.','Error');
end
% Compartment coordinates
Xc = locs(:,1);
Yc = locs(:,2);
Zc = locs(:,3);
% Parent coordinates
Xp = locs(:,4);
Yp = locs(:,5);
Zp = locs(:,6);
% Plot
fig = figure('units','normalized','outerposition',[0 0 1 1]);
h = plot3([Xc Xp]',[Yc Yp]',[Zc Zp]','k','LineWidth',1.5);
view(0,90)
title('Neuron Morphology');
xlabel('X (\mum)');
ylabel('Y (\mum)');
zlabel('Z (\mum)');
axis equal
%% Receive neuron axis and desired orientation
prompt = {'\bfEnter the cartesian orientation of the somato-dendritic axis:\rm from soma to apical dendrites (X Y Z with space in between)',...
    '\bfEnter the cartesian orientation of the desired neuron direction:\rm if left empty, the neuron is placed perpendicular to the surface (X Y Z with space in between)'};
dlgtitle = 'Neuron orientation';
dims = [1 100];
opts.Interpreter = 'tex';
opts.WindowStyle = 'normal';
definput = {'0 1 0',''};
while 1
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    nrnaxs = str2double(split(answer{1}));
    if isempty(answer{2})
        nrnori = [];
    else
        nrnori = str2double(split(answer{2}));
    end
    if length(nrnaxs) == 3 && sum(isnan(nrnaxs)) == 0 && (length(nrnori) == 3 || isempty(nrnori)) && sum(isnan(nrnori)) == 0
        break
    end
    definput = answer;
    if length(nrnaxs) ~= 3 || sum(isnan(nrnaxs)) ~= 0
        definput{1} = 'Wrong format!';
    end
    if ~(length(nrnori) == 3 || isempty(nrnori)) || sum(isnan(nrnori)) > 0
        definput{2} = 'Wrong format!';
    end
end
nrnaxs = nrnaxs';
nrnori = nrnori';

try
    close(fig);
catch
end
%% Translation and rotation
% Translates the neuron to the desired location (nrnloc) and depth (nrndpth)
% and rotates neuron from its dendro-somatic axis (nrnaxs) to align with the
% desired orientation if specified (nrnori), or otherwise the surface normal (v2)
% at the desired location.
v1 = nrnaxs./normd(nrnaxs);
normals = mesh_get_triangle_normals_new(gm_surf);
[~,tri_idx] = mesh_get_closest_triangle_from_point(gm_surf,nrnloc);
v2 = normals(tri_idx,:);
v2 = v2./normd(v2);

nrnloc2 = nrnloc - v2*nrndpth;

if ~isempty(nrnori)
    v2 = nrnori;
    v2 = v2./normd(v2);
end

theta = acos(v1*v2'); % Angle between v1 and v2
axis_rot = cross(v1,v2)/norm(cross(v1,v2)); % axis around which v1 should turn to align with v2

% A skew symmetric representation of the normalized axis
axis_skewed = [0 -axis_rot(3) axis_rot(2); axis_rot(3) 0 -axis_rot(1); -axis_rot(2) axis_rot(1) 0];

% Rodrigues' formula for the rotation matrix
R = eye(3) + sin(theta)*axis_skewed + (1-cos(theta))*axis_skewed^2;
rotmat = [R, nrnloc2'; 0 0 0 1];
c_trans = apply_4D_transformation_matrix([Xc Yc Zc]/1000, rotmat); % Transformed child nodes coordinates (in mm)
p_trans = apply_4D_transformation_matrix([Xp Yp Zp]/1000, rotmat); % Transformed parent nodes coordinates (in mm)
%% Interpolate E-field at the neuron components
E_int = get_fields_at_coordinates(mesh,c_trans);
% Extract the values for vector E-field
% for idx = 1:length(E_int)
%     if strcmp(E_int{idx}.name, 'E')
%         break
%     end
% end
% if ~strcmp(E_int{idx}.name, 'E')
%     error(['The input mesh file from SimNIBS does not include electric field. ' ...
%         'Make sure to select "vector E" in "Simulation Options" in SimNIBS.']);
% end
% E_int = E_int{idx}.data;
E_int = E_int{1}.data;
%% Scale E-field
prompt = {'\bfScaling coefficient for the electric field:\rm Leave as 1 if you do not want to scale'};
dlgtitle = 'Scale E-field';
dims = [1 100];
opts.Interpreter = 'tex';
opts.WindowStyle = 'normal';
definput = {'1'};
while 1
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    scale_E = str2double(answer);
    if ~isnan(scale_E)
        break
    end
    definput = {'Wrong format!'};
end

E_int = scale_E * E_int;
%% Calculate quasipotentials
% Wang et al., 2018, "Coupling Magnetically Induced Electric Fields to Neurons: Longitudinal
% and Transverse Activation"

psi = nan(size(E_int,1),1); % Quasipotentials
psi(1) = 0;

[~, parent_id] = ismember(p_trans,c_trans,'rows');

% Since the compartments may not be in order, it may be necessary to loop
% over the compartments several times, until quasipotential is calculated
% for all of them
while any(isnan(psi))
    for ii = 2:size(c_trans,1)
        if ~isnan(psi(parent_id(ii)))
            Spc = c_trans(ii,:) - p_trans(ii,:); % Displacement vectors from parents to children
            psi(ii) = psi(parent_id(ii)) - dot((E_int(ii,:)+E_int(parent_id(ii),:))/2,Spc);
        end
    end
end

%% Trim mesh
% For Faster and easier visualization, trim the mesh to the region of
% interest
radius = 20; % Keep the mesh within this radius (mm) of the neuron location
dist = get_distance(gm_surf.nodes, nrnloc);
nodes_idx = dist < radius;
gm_surf_trim = mesh_extract_points(gm_surf,nodes_idx);
%% Export results
% Export the quasipotentials in a txt file for NEURON and the mesh file for
% visualization in gmsh
respath = fullfile('..','..','Results','E-field_Coupling');
if ~exist(respath,'dir')
    mkdir(respath);
end

lines = zeros(size(c_trans,1),2);
for ii = 1:size(c_trans,1)
    [~, idx] = ismember(p_trans(ii,:),c_trans,'rows');
    lines(ii,:) = [ii idx];
end
nrn_mesh.points = c_trans;
nrn_mesh.lines = lines(2:end,:);

nrn_mesh.element_data{1}.data = psi(2:end);
nrn_mesh.element_data{1}.name = 'Quasipotentials';

nrn_mesh.node_data{1}.data = E_int;
nrn_mesh.node_data{1}.name = 'E';
mesh_save_gmsh(nrn_mesh,[respath filesep 'quasipotentials.msh']);

save([respath filesep 'quasipotentials.txt'], 'psi','-ascii');
% write_values_in_mesh(psi, c_trans, [respath file]);

mesh_save_gmsh4(gm_surf_trim,[respath filesep 'mesh_trim.msh']);
%% Copy SimNIBS results
if ~exist([respath filesep 'SimNIBS_Simulation'],'dir')
    mkdir([respath filesep 'SimNIBS_Simulation']);
else
    delete([respath filesep 'SimNIBS_Simulation' filesep '*'])
end
copyfile(meshpath, [respath filesep 'SimNIBS_Simulation'], 'f');
%% Export parameters
fid = fopen([respath filesep 'parameters.txt'],'wt');
fprintf(fid, '%% FEM mesh file name\n');
fprintf(fid, ['meshfile = ''' meshfile ''';\n']);
fprintf(fid, '%% FEM mesh pathway\n');
fprintf(fid, ['meshpath = ''' strrep(meshpath,'\','\\') ''';\n']);
fprintf(fid, '%% neuron location\n');
fprintf(fid, ['nrnloc = [' num2str(nrnloc) '];\n']);
fprintf(fid, '%% neuron depth\n');
fprintf(fid, ['nrndpth = ' num2str(nrndpth) ';\n']);
fprintf(fid, '%% neuron segment coordinates file name\n');
fprintf(fid, ['nrnfile = ''' strrep(nrnfile,'\','\\') ''';\n']);
fprintf(fid, '%% neuron segment coordinates pathway\n');
fprintf(fid, ['nrnpath = ''' strrep(nrnpath,'\','\\') ''';\n']);
fprintf(fid, '%% neuron axis\n');
fprintf(fid, ['nrnaxs = [' num2str(nrnaxs) '];\n']);
fprintf(fid, '%% neuron desired oriention\n');
if isempty(nrnori)
    fprintf(fid, 'nrnori = [];\n');
else
    fprintf(fid, ['nrnori = [' num2str(nrnori) '];\n']);
end
fprintf(fid, '%% E-field scaling factor\n');
fprintf(fid, ['scale_E = ' num2str(scale_E) ';\n']);
fprintf(fid, '%% results directory\n');
fprintf(fid, ['respath = ''' strrep(respath,'\','\\') ''';\n']);
fclose(fid);

msgbox('Coupling successful.','Success');
end
