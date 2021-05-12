function [m]=mesh_extract_points(m,chunk_points_logical,varargin)
  % Returns a mesh which contains only the selected points (i.e., nodes) 
  % and affected simplices and data.
  % USAGE: [m] = MESH_EXTRACT_POINTS(m,chunk_points_logical)
  % Mirko Windhoff, 2009
  % updated by Alexander Opitz, 2016
  % modified by Ivan Alekseichuk, 2018
  
  if nargin > 2
      extract_incompletely_selected_elements = varargin{1}; 
  else
      extract_incompletely_selected_elements = false; 
  end
  
  if size(chunk_points_logical,1) ~= size(m.nodes,1) 
      error('Dimension mismatch'); 
  end
  
  % generate map for points indices to the reduce number of points
  ab_map = zeros(size(m.nodes,1),1);
  ab_map(chunk_points_logical) = 1:sum(chunk_points_logical);
  
  % get selection of simplices, that are formed by the selected points
  if extract_incompletely_selected_elements
    chunk_triangles_logical = any(ab_map(m.triangles)>0,2);
    chunk_tetrahedra_logical = any(ab_map(m.tetrahedra)>0,2);
    chunk_points_logical = linear2logical(unique(m.tetrahedra(chunk_tetrahedra_logical,:)), [size(m.nodes,1),1]);
    ab_map = zeros(size(m.nodes,1),1);
    ab_map(chunk_points_logical) = 1:sum(chunk_points_logical);
  else
    chunk_triangles_logical = all(ab_map(m.triangles)>0,2);
    chunk_tetrahedra_logical = all(ab_map(m.tetrahedra)>0,2);
  end
  
  if ~isempty(m.element_data)
    if length(m.tetrahedra) == length(m.element_data{1,1}.tetdata)
          istet = 1;
    else istet = 0;
    end
  
    if length(m.triangles) == length(m.element_data{1,1}.tridata)
          istri = 1;
    else istri = 0;
    end
  end
      
  m.nodes = m.nodes(chunk_points_logical,:);
  m.triangles = ab_map(m.triangles(chunk_triangles_logical,:));
  m.triangle_regions = m.triangle_regions(chunk_triangles_logical);
  if numel(m.triangles) == 3
      m.triangles = m.triangles';
  end
  m.tetrahedra = ab_map(m.tetrahedra(chunk_tetrahedra_logical,:));
  m.tetrahedron_regions = m.tetrahedron_regions(chunk_tetrahedra_logical);
  if numel(m.tetrahedra) == 4
      m.tetrahedra = m.tetrahedra';
  end
  if isfield(m, 'hexahedra')
    if extract_incompletely_selected_elements
      chunk_hexahedra_logical = any(ab_map(m.hexahedra)>0,2);
    else
      chunk_hexahedra_logical = all(ab_map(m.hexahedra)>0,2);
    end
    m.hexahedra = ab_map(m.hexahedra(chunk_hexahedra_logical,:));
    m.hexahedron_regions = m.hexahedron_regions(chunk_hexahedra_logical);
    if numel(m.hexahedra) == 8
        m.hexahedra = m.hexahedra'; 
    end
  end
  
  for j=1:size(m.node_data,1)
    if ~isfield(m.node_data{j},'ids')
      m.node_data{j}.data(~chunk_points_logical,:)=[];
    else
      [m.node_data{j}.ids, ia]=intersect(m.node_data{j}.ids,find(chunk_points_logical));
      m.node_data{j}.data=m.node_data{j}.data(ia,:);
    end
  end
  
  for j=1:size(m.element_data,1)
    if isfield(m, 'hexahedra'), warning('element values on hexahedra not supported by this function, make sure that the element_data is correctly extracted!'); end;
    if ~isfield(m.element_data{j},'ids')
      if istet
        m.element_data{j}.tetdata(~chunk_tetrahedra_logical,:)=[];
      end
      if istri
        m.element_data{j}.tridata(~chunk_triangles_logical,:)=[];
      end
    else
      [m.element_data{j}.ids, ia]=intersect(m.element_data{j}.ids,find([false(size(m.triangles,1),1);chunk_tetrahedra_logical]));
      m.element_data{j}.data=m.element_data{j}.data(ia,:);
    end
  end
  
end