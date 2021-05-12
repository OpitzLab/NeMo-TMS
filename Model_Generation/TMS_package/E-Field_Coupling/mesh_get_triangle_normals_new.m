function [triangle_normals]=mesh_get_triangle_normals_new(m)
  % Calcuates the normals for each triangles. If the orientation of the
  % triangle edges is ccw to a viewer, then the normal points to the viewer.
  % USAGE: triangle_normals=MESH_GET_TRIANGLE_NORMALS(m)
  % Mirko Windhoff, 2009
  % $Id: mesh_get_triangle_normals.m 255 2009-09-14 17:19:25Z mirko $
  % Modified by Sina Shirinpour 1/16/2019
  a=m.nodes(m.triangles(:,1),:)-m.nodes(m.triangles(:,2),:);
  b=m.nodes(m.triangles(:,1),:)-m.nodes(m.triangles(:,3),:);
  triangle_normals=cross(a,b,2);
  triangle_normals = bsxfun(@rdivide, triangle_normals, normd(triangle_normals));
end