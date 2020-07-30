function dist = get_distance(points, refpoints)
% this function calculates the distance of every point to a given set of
% reference_points
% Alexander Opitz, 2018
% modified by Ivan Alekseichuk, 2018
%
% points = list of points, refpoints = reference points

dist = zeros(size(points,1),1);
for i = 1: size(points,1)
  dist(i) = min(sum((refpoints - repmat(points(i,:),size(refpoints,1),1)).^2,2));
%   dist(i) = min(normd(m.points - repmat(points(i,:),size(m.points,1),1)));
end
dist = sqrt(dist);