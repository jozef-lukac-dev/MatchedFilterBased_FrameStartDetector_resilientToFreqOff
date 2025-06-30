function allTuples = getAllTuples(varargin)
% function GETALLTUPLES
% allTuples = getAllTuples(rng1,rng2,...)
% Get all tuples from the ranges given as inputs.
% Inputs:
%  <rng_i> -- vector of (n_i) values.
%
% n = n_1 * ... * n_N 
%
% Outputs:
%  <allTuples> -- n x N matrix of tuples.
%

N = nargin;
n_i_vect = cell2mat(cellfun(@(x) numel(x), varargin,'UniformOutput',false));
n = prod(n_i_vect);

allTuples = zeros(n,N);
n_cumProd = flip(cumprod([1,flip(n_i_vect)]));
for i_src = 1:N
    N_rep = n/n_i_vect(i_src);
    n_rep_row = n_cumProd(i_src+1);
    auxVect = repmat( varargin{i_src}(:)', n_rep_row, N_rep/n_rep_row);
    allTuples(:,i_src) = auxVect(:);
end

end