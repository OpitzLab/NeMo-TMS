function [m, chunk_points_logical]=mesh_extract_triangle_regions(m, triangle_regions, varargin)
  % Removes all points, triangles and data that are not contained in that
  % triangle_region. Also all tetrahedra are removed. You can use absolute 
  % region numbers by disabling the relative_numbers flag. 
  % relative_numbers are the indices of the ascending sorted region number list.
  % USAGE: m=MESH_EXTRACT_TRIANGLE_REGIONS(m, triangle_regions[, relative_numbers=true])
  % Mirko Windhoff, 2009
  % $Id: mesh_extract_triangle_regions.m 261 2009-09-21 10:13:02Z mirko $
  if nargin>2, relative_numbers=logical(varargin{1}); else relative_numbers=true; end;
  if relative_numbers
    regions=unique(m.triangle_regions);
    triangle_regions=regions(triangle_regions);
  end;
  chunk_triangles_logical=false(size(m.triangle_regions,1),1);
  for i=1:numel(triangle_regions)
    chunk_triangles_logical = chunk_triangles_logical | m.triangle_regions==triangle_regions(i);
  end;
  chunk_points_logical=false(size(m.points,1),1);
  chunk_points_logical(unique(m.triangles(chunk_triangles_logical,:)))=true;
  m=mesh_extract_points(m,chunk_points_logical);
  m.tetrahedra=[];
  m.tetrahedron_regions=[];
end