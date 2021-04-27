function mesh_save_gmsh(m,filename)
  % Write the mesh struct m to a Gmsh mesh file of.
  % USAGE:
  % MESH_SAVE_GMSH(mesh_struct, filename)
  % mesh_struct can contain fields:
  %   - points (np x 3)
  %   - lines (nl x 2)
  %   - line_regions (nl x 1)
  %   - triangles (ntr x 3)
  %   - triangle_regions (ntr x 1)
  %   - tetrahedra (nth x 4)
  %   - tetrahedron_regions (nth x 1)
  %   - hexahedra (nh x 8)
  %   - hexahedron_regions (nh x 1)
  %   - node_data (np x m, m = 1 or 3; optional, data for each point)
  %         node_data={struct('name', name_for_data_set, 'data',
  %           values_for_each_node[, 'ids', node_ids_for_data]); ...};
  %         values_for_each_node = numNodes x n, values can be scalar, vector...
  %   - element_data (nth x m or ntr+nth x m, m = 1 or 3; optional, data for each tetrahedron or for each triangle and tetrahedron)
  %         element_data={struct('name', name_for_data_set, 'data',
  %           values_for_each_element[, 'ids', element_ids_for_data]); ...};
  %         values_for_each_element = numElements x {1}
  %         if ids is given, the number of values must match the number of
  %         ids.
  %   - meta (optional)
  % See http://www.geuz.org/gmsh/doc/texinfo/gmsh-full.html#SEC56 for
  % documentation of the file format.
  %
  % Mirko Windhoff, 2009
  % $Id: mesh_save_gmsh.m 513 2010-09-29 11:35:55Z mwindhoff $
  if ~isfield(m, 'points'), error('Can''t write mesh without points'); end;
  
  if ~isfield(m, 'lines'), m.lines=[]; end;
  if ~isfield(m, 'line_regions'), m.line_regions=zeros(size(m.lines,1),1); end;
  
  if ~isfield(m, 'triangles'), m.triangles=[]; end;
  if ~isfield(m, 'triangle_regions'), m.triangle_regions=zeros(size(m.triangles,1),1); end;
  
  if ~isfield(m, 'tetrahedra'), m.tetrahedra=[]; end;
  if ~isfield(m, 'tetrahedron_regions'), m.tetrahedron_regions=zeros(size(m.tetrahedra,1),1); end;

  if ~isfield(m, 'hexahedra'), m.hexahedra=[]; end;
  if ~isfield(m, 'hexahedron_regions'), m.hexahedron_regions=zeros(size(m.hexahedra,1),1); end;
  
  if ~isfield(m, 'node_data'), m.node_data={}; end;
  if ~isfield(m, 'element_data'), m.element_data={}; end;
  if ~isfield(m, 'meta'), m.meta=struct(); end;
  
  % append .msh if needed
  [fpath,fname,fext]=fileparts(filename);
  if ~strcmp(fext,'.msh'), filename = [filename,'.msh']; end;

  fid = fopen(filename,'wt');
  if fid<0, error(['Could not open file ', filename, ' for writing']); end;

  % header information
  if ~fprintf(fid, '$MeshFormat\n2.0 0 8\n$EndMeshFormat\n$Nodes\n')
    error(['Could not write to file ', filename]);
  end;
  fprintf(fid, '%d\n', size(m.points, 1));

  np =size(m.points,1);
  nl =size(m.lines,1);
  ntr=size(m.triangles,1);
  nth=size(m.tetrahedra,1);
  nh =size(m.hexahedra,1);
  % points
  fprintf(fid, '%d %f %f %f\n', [1:np; m.points']);
  fprintf(fid, '$EndNodes\n$Elements\n');
  fprintf(fid, '%d\n', nl + ntr + nth + nh);
  n=0;
  % lines
  if nl>0
    fprintf(fid, '%d 1 2 %d %d %d %d\n', [n+(1:nl); m.line_regions'; m.line_regions'; m.lines']);
  end;
  n=n+nl;
  % triangles
  if ntr>0
    fprintf(fid, '%d 2 2 %d %d %d %d %d\n', [n+(1:ntr); m.triangle_regions'; m.triangle_regions'; m.triangles']);
  end;
  n=n+ntr;
  % tetrahedra
  if nth>0
    m.tetrahedra=double(m.tetrahedra);
    m.tetrahedron_regions=double(m.tetrahedron_regions);
    fprintf(fid, '%d 4 2 %d %d %d %d %d %d\n', [double(n+(1:nth)); m.tetrahedron_regions'; m.tetrahedron_regions'; m.tetrahedra']);
  end;
  n=n+nth;
  % hexahedra
  if nh>0
    m.hexahedra=double(m.hexahedra); % fix that should be elswhere...
    m.hexahedron_regions=double(m.hexahedron_regions); % fix that should be elswhere...
    fprintf(fid, '%d 5 2 %d %d %d %d %d %d %d %d %d %d\n', [double(n+(1:nh)); m.hexahedron_regions'; m.hexahedron_regions'; m.hexahedra']);
  end;
  fprintf(fid, '$EndElements\n');

  if isfield(m, 'node_data')
    for i=1:size(m.node_data,1)
      values=m.node_data{i}.data;
      if isfield(m.node_data{i},'ids')
        if size(m.node_data{i}.ids,1)~=size(values,1), size(m.node_data{i}.ids,1), size(values,1), error('sizes(m.node_data{%d}.data,1)~=sizes(m.node_data{%d}.ids,1)',i,i); end;
        ids=m.node_data{i}.ids.';
      else
        ids=1:size(values,1);
      end;
      fprintf(fid, '$NodeData\n1\n"%s"\n1\n0.0\n3\n0\n%d\n%d\n', m.node_data{i}.name, size(values, 2), size(values, 1));
      pattern='%d';
      for l=1:size(values, 2), pattern=[pattern,' %e']; end;
      pattern=[pattern,'\n'];
      fprintf(fid, pattern, [ids; values']);
      fprintf(fid, '$EndNodeData\n');
    end;
  end;

  if isfield(m, 'element_data')
    for i=1:size(m.element_data,1)
      if isfield(m.element_data{i},'ids')
        if size(m.element_data{i}.ids,1)~=size(m.element_data{i}.data,1), error('size(m.element_data{%d}.ids,1)~=size(m.element_data{%d}.data,1)',i,i); end;
        ids=m.element_data{i}.ids.';
      elseif size(m.element_data{i}.data,1)==nth % only values on tetrahedra given
        ids=(1:nth)+nl+ntr;
      elseif size(m.element_data{i}.data,1)==nh % only values on hexahedra given
        ids=(1:nh)+nl+ntr+nth;
      elseif size(m.element_data{i}.data,1)==nth+nh % only values on tetrahedra and hexahedra given
        ids=(1:nth+nh)+nl+ntr;
      elseif size(m.element_data{i}.data,1)==ntr+nth+nh  % only values on triangles, tetrahedra and hexahedra given
        ids=1:(ntr+nth+nh)+nl;
      else
        if size(m.element_data{i}.data,1)~=nl+ntr+nth+nh, error('element_data has wrong size'); end;
        ids=1:(nl+ntr+nth+nh);
      end;
      values=m.element_data{i}.data;
%       fprintf(fid, '$ElementData\n1\n"%s"\n1\n0.0\n3\n0\n1\n%d\n', m.element_data{i}.name, size(values,1));
      fprintf(fid, '$ElementData\n1\n"%s"\n1\n0.0\n3\n0\n%d\n%d\n', m.element_data{i}.name, size(values,2), size(values,1));
      pattern='%d';
      for l=1:size(values, 2), pattern=[pattern,' %e']; end;
      pattern=[pattern,'\n'];
      fprintf(fid, pattern, [ids; values']);
      fprintf(fid, '$EndElementData\n');
    end;
  end;
  fclose(fid);
end
