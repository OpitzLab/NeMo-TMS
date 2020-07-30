function [norm_values] = normd(in_vector, varargin)
  % calculates the square p-norm along norm_dimension
  % With default values for it is like sqrt(sum(in_vector.^2,2)).
  % Mirko Windhoff, 2009
  % $Id: normd.m 432 2010-08-05 11:06:14Z mwindhoff $
  % USAGE: [norm_values]=normd(in_vector[, p=2, norm_dimension=2])
  
  p = 2;
  norm_dimension = 2;
  
  if nargin > 1
      p = varargin{1};
  end
  if nargin > 2
      norm_dimension = varargin{2};
  end
  
  norm_values = (sum(in_vector .^ p, norm_dimension)) .^ (1/p);
  
end