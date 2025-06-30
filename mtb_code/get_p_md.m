function [p_md1, p_md2] = get_p_md(snr_dB_rng, h, N_seq, p_fa, F_eps, use_perAprox)
% function GET_P_MD
% [p_md1, p_md2] = get_p_md(snr_dB_rng, h, N_seq, p_fa, F_eps [,use_perAprox])
% Computes probability of misdetection using y'_1[n] metric.
% Inputs:
%  <snr_dB_rng> -- N_snr x 1 vector, SNR in dB
%  <h>  -- N_h element vector, preamble
%  <N_seq> -- 1 x 1 number, number of sub-sequences.
%  <p_fa> -- 1 x 1 number, false-alarm probability.
%  <F_eps> -- 1 x 1 number, relative frequency offset (relative to sampling
%  rate)
%  <use_perAprox> -- 1x1 boolean, whether to use perdiodic approximation
%  (default is false)
% Outputs:
%  <p_md1> -- N_snr x 1 vector, probability of misdetection p_{md,1}.
%  <p_md2> -- N_snr x 1 vector, probability of misdetection p_{md,2}.
%
if nargin < 6
    use_perAprox = false;
end
phi_norm = @(x) 1/2*(1 + erf(x/sqrt(2)));
lHalf = @(x) exp(x/2).*( ...
  (1-x).*besseli(0,-x/2) - x.* besseli(1,-x/2) );

h = h ./ sqrt( sum(abs(h(:)).^2) );
L_seq = ceil(numel(h)/N_seq);
N_h = numel(h);
N_add = L_seq * N_seq - N_h;
h_parts = reshape([h(:); zeros(N_add,1)], L_seq, N_seq); %L_seq x N_seq

F_h_parts_sq = abs(...
    sum(abs(h_parts).^2 .* exp(1i*2*pi*F_eps *(0:L_seq-1).'), 1)...
    ).^2; % 1 x N_seq vector
if use_perAprox
    F_h_parts_sq = repmat(F_h_parts_sq(1), 1, N_seq);
end
sigma_w2 = 10.^(-snr_dB_rng/10); %noise variance

arg_aux = ( ...
  sqrt( (gamma(N_seq+1)/gamma(N_seq+ 1/2))^2 * ...
    gammaincinv(p_fa, N_seq,'upper') )...
   -sum(lHalf(-F_h_parts_sq*N_seq./sigma_w2),2)...
  )  ./   ...
  ( sqrt( 4/pi*N_seq   + 4/pi ./sigma_w2 .* ...
     sum(F_h_parts_sq * N_seq, 2) ...
    -sum(lHalf(-F_h_parts_sq*N_seq./sigma_w2).^2,2)...
    ) ...
  );
p_md1 = phi_norm(arg_aux); %probability of misdetection

p_md2 = 1 - marcumq(...
  sqrt(2*N_seq./sigma_w2*sum(F_h_parts_sq,2)),...
  sqrt(2 * gammaincinv(p_fa, N_seq,'upper')),...
  N_seq);
end