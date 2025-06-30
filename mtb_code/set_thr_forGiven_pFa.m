function [r_thr, p_fa_final, iIter] = set_thr_forGiven_pFa(...
    h, N_seq, p_fa_req, abs_absSq, mAv_mMed, N_samps, L_av)
%[r_thr, p_fa_final, iIter] = set_thr_forGiven_pFa(...
 %   h, N_seq, p_fa_req, abs_absSq, mAv_mMed, N_samps, L_av)
% Set the noise-independent threshold to ensure the given/required
% <p_fa_req> -- false-alarm probability.
% Inputs:
%  <h> -- N_h element vector. h[n] preamble
%  <N_seq> -- number subsequences
%  <p_fa_req> -- 1>real number>0. Required false-alarm probability.
%  <abs_absSq> -- char-array 'abs' or 'absSq'. Taking just abs(.) or
%   abs(.).^2.
%  <mAv_mMed> -- char-array 'mAv' or 'mMed' -- moving average or moving
%   median.
%  <N_samps> -- natural number > 0. Number of testing samples.
%  <L_av> -- natural number >0. Averaging length.
% Outputs:
%  <r_thr_final> -- real number > 0,

N_h = numel(h);
snr_dB = 7; %dB   a random value, result doesn't depend on it
L_seq = ceil(N_h / N_seq);
idx_off = N_h + L_seq; %offset index
y_ = MF_correlate(h, ...
    sqrt(1 * 10^(-snr_dB/10)/2) * ...
        (randn(N_samps,1) + 1i*randn(N_samps,1)),...
    N_seq, abs_absSq); %remove several first samples, transient

if strcmp(mAv_mMed,'mAv')
    mAv_mMed_filt = @(x) filter(1/(L_av)*ones(L_av,1),1, x);
else
    mAv_mMed_filt = @(x) ...
      movmedian([zeros(L_av-1,1); x],L_av,'Endpoints','discard');
end
z = mAv_mMed_filt( y_ ); 

r_thr = 50; %a sufficiently large number
r_thr_delta = r_thr/2;
iIter = 1;
p_fa_final = 1;
tol = 1e-3;
N_iter_max = 25;
while true
    p_fa_final = sum( y_(idx_off:end)...
       > r_thr * z(idx_off-1:end-1), 1)...
        / (N_samps - idx_off +1);
    if (abs(log10(p_fa_final) - log10(p_fa_req)) < tol) ...
            || (iIter > N_iter_max)
        break
    end
    if p_fa_final > p_fa_req
        r_thr = r_thr + r_thr_delta;
    else
        r_thr = r_thr - r_thr_delta;
    end
    r_thr_delta = r_thr_delta /2;
    iIter = iIter + 1;
end
end
