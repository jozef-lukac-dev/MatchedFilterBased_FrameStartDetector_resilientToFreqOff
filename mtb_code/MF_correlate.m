function corrMetric = MF_correlate(h_vec, rx_samps, N_seq,  abs_or_absSq)
% function MF_CORRELATE
% corrMetric = MF_correlate(h_vec, rx_samps, N_seq, abs_or_absSq)
% Do a crosscorrelation of the received samps <rx_samps> with FIR filter
% coefficients <h_preamb_samps>. Divide the preamble to <N_seq> sequences 
% (if neccessary, the last sequence is filled with zeros).
% Inputs:
%  <h_vec> -- N_h-element vector, preamble sequence.
%  <rx_samps> -- N_x x N_frames matrix, column vectors of received samples.
%  <N_seq> -- 1 x 1 number, number of subsequences.
%  <abs_or_absSq> -- char-vector, 'abs' or 'absSq' (default). Phase
%  marginalization operation.
% Outputs:
%  <corrMetric> -- N_x x N_frames  vector of correlation metric values.
%

if nargin < 4
    abs_or_absSq = 'absSq';
end
assert( strcmp(abs_or_absSq,'absSq') ...
      || strcmp(abs_or_absSq,'abs'), ...
   ['abs_or_absSq should be either',...
   ' ''abs'' or ''absSq''.\n']);
if strcmp(abs_or_absSq,'absSq')
    sqFun = @(x) abs(x).^2;
else
    sqFun = @(x) abs(x);
end

N_h = numel(h_vec);
L_seq = ceil(N_h / N_seq);
N_h_ext = N_seq * L_seq; % extended length
N_add = N_h_ext - N_h; %number of added zeros at the end of h[n]
h_parts = reshape([h_vec(:); zeros(N_add,1)], L_seq, N_seq);

[N_x, N_frames] = size(rx_samps);

sz_rx_N_h = [N_h_ext-1, N_frames]; %no spaces around +- !!! 
corrMetric = zeros([N_x+N_h_ext-1, N_frames]);
for iPart = 1:N_seq
    metric_part = sqFun( ...
        filter( conj(flip(h_parts(:,iPart))), 1, ...
          [rx_samps; zeros(sz_rx_N_h)]) ...
      );
    N_init_zeros = (N_seq-iPart) * L_seq;
    corrMetric = corrMetric + ...
        [zeros([N_init_zeros, N_frames]);
         metric_part(1:end-N_init_zeros, :)];
end
corrMetric = corrMetric(N_add+1:N_add+N_x, :);
end
