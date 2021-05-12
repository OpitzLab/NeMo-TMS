function couple_script(paramfile)
%% Import parameters
fid = fopen(paramfile);
tline = fgetl(fid);
while ischar(tline)
    eval(tline);
    tline = fgetl(fid);
end
fclose(fid);
%% Load mesh
mesh = mesh_load_gmsh4([meshpath meshfile]);
gm_surf = mesh_extract_regions(mesh,'elemtype','tri','region_idx',1002);
%% Load neuron
try
    locs = load([nrnpath nrnfile]);
catch
    error('Export NEURON model segments first.');
end
% Compartment coordinates
Xc = locs(:,1);
Yc = locs(:,2);
Zc = locs(:,3);
% Parent coordinates
Xp = locs(:,4);
Yp = locs(:,5);
Zp = locs(:,6);
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
for idx = 1:length(E_int)
    if strcmp(E_int{idx}.name, 'E')
        break
    end
end
if ~strcmp(E_int{idx}.name, 'E')
    error(['The input mesh file from SimNIBS does not include electric field. ' ...
        'Make sure to select "vector E" in "Simulation Options" in SimNIBS.']);
end
E_int = E_int{idx}.data;
%% Scale E-field
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

% write_values_in_mesh(psi, c_trans,'D:\Neuron\Coupling\test.txt');
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

if(~exist(respath,'dir'))
    mkdir(respath); % create results folder
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

disp('Coupling successful.');
end
