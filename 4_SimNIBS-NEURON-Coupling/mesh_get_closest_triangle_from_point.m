function [minval,idx] = mesh_get_closest_triangle_from_point(m, point, varargin)
  % Gives the linear index to the closest triangle of the surface specified by triangle_regions to the point.
  % USAGE: idx=MESH_GET_CLOSEST_TRIANGLE_FROM_POINT(m, point[, triangle_regions=1])
  % Mirko Windhoff, 2010
  % $Id: mesh_get_closest_triangle_from_point.m 625 2010-12-23 18:24:22Z mwindhoff $
  % modified by Ivan Alekseichuk, 2018
  
  if nargin > 2 
      triangle_region = varargin{1};
      m = mesh_extract_triangle_regions(m, triangle_region);
  end  
  
  center       = mesh_get_triangle_centers(m);
  [minval,idx] = min(normd(center-repmat(point, size(center,1), 1)));
end
