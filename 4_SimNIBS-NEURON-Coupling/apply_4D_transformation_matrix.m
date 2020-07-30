function p = apply_4D_transformation_matrix(points, matrix)
  % Applys a 4x4 transformation matrix (rotation + translation) to the
  % points (http://en.wikipedia.org/wiki/Homogeneous_coordinates).
  % points have Nx3 dim
  % Alexander Opitz, 2010
  % $Id$
  if any(size(matrix)~=4)
      error('matrix has to be 4x4 sized.'); 
  end
  p=[points,ones(size(points,1),1)]';
  p=(matrix*p)';
  p=p(:,1:3);