function p_d = p_d_y1_nEst(snr_dB_rng, p_fa, h, N_seq, fOff, N_frames, N_rep, L_av)

% nEst -- noise (variance) estimation
N_h = numel(h);
N_h0 = ceil(N_h / N_seq);

if L_av == 0
    L_av = N_h0 -1;
end

assert(abs(sum(abs(h).^2)- 1) < 100*eps,'h should be normalized to unit energy.');
assert(all(diff(snr_dB_rng) > 0),'snr_dB range should be increasing vector.');

N_snr = numel(snr_dB_rng);
sigma_w2 = 1 * 10.^(-snr_dB_rng/10);
r_opt_coeff = sqrt(  (gamma(N_seq +1)/gamma(N_seq + 1/2))^2 * pi/4 *...
    1/N_seq  * gammaincinv(p_fa, N_seq, 'upper')  );
N_all_tests = N_rep * N_frames;
p_d = N_all_tests * ones(N_snr,2);

idx_mid = 2*N_h;
idx_rng = (idx_mid-(N_h-1)):(idx_mid+(N_h-1));

idx_mid2 = N_h + L_av;
idx_rng2 = (idx_mid2-(N_h-1)):(idx_mid2+(N_h-1));

k_idc = (0:N_h-1).';
freq_off_factor = exp(1i*2*pi*fOff * k_idc);
for iSnr = 1:N_snr
    p_d_mid_val = 0;
    p_d_all_vals = 0;
    sigmaw2_curr = sigma_w2(iSnr);
    parfor iRep = 1:N_rep
        y_ = MF_correlate(h, ...
                    [zeros(N_h,N_frames); ...
  freq_off_factor .* repmat(h, 1, N_frames); ...
                    zeros(N_h,N_frames)] + ...
                    sqrt(sigmaw2_curr/2) * ...
                        (randn(3*N_h,N_frames) + 1i*randn(3*N_h,N_frames)),...
                N_h0, 'abs');
        z_nVar_est = filter(1/(L_av)*ones(L_av,1),1, ...
            abs( ...
                sqrt(sigmaw2_curr/2) * ...
                 (randn(2*N_h + L_av,N_frames) + 1i*randn(2*N_h + L_av,N_frames)) ...
               ).^2 ...
           );
        p_d_mid_val = p_d_mid_val + sum(y_(idx_mid,:) > r_opt_coeff * sqrt(z_nVar_est(idx_mid2, :)),2);
        p_d_all_vals = p_d_all_vals + sum(any(y_(idx_rng,:) > r_opt_coeff * sqrt(z_nVar_est(idx_rng2, :)),1),2);
    end
    if p_d_mid_val == N_all_tests && p_d_all_vals == N_all_tests
        break;
    end
    p_d(iSnr,:) = [p_d_mid_val, p_d_all_vals];
end
p_d = p_d / N_all_tests;
end
