outDir = fullfile('..','data_out_probs');

h_preambs_params = {... %---------------
    {'wifi',16,10},...
    {'ZC',100,1,1},{'ZC',100,1,2},{'PN',100,1},...
    {'ZC',200,1,1},{'ZC',200,1,2},{'PN',200,1},...
    {'ZC',400,1,1},{'ZC',400,1,2},{'PN',400,1}...
    ...{'wifi',N_h0,5},{'ZC',N_h0,5,1},{'ZC',N_h0,5,2},{'ZC',N_h0,5,5},{'PN',N_h0,5};
    };
% h_preambs_params = {{'PN',10,1}}; 
h_preambs = cellfun(@(cellInput) get_preamble(cellInput{:}),...
    h_preambs_params,'UniformOutput',false);

p_fa_rng = 10.^(-3:-1:-4).';

snr_dB_rng = [8, 10];
% snr_dB_rng = (8:1.5:30.5).'; %---------------

freq_rng = 0.001;
% freq_rng = (0:0.001:0.01).';
N_freq = numel(freq_rng);

% N_seq_rng = [1, 2, 3, 10].'; %---------------
N_seq_rng = [1, 2].'; 
N_seq_rng_numel = numel(N_seq_rng);


% L_rng = [10, 20, 30, 40, 0]; %---------------
L_rng = [10];
N_L = numel(L_rng);

N_cases = numel(h_preambs);
N_snr = numel(snr_dB_rng);
N_fa = numel(p_fa_rng);

% [snr_dB_grid1, p_fa1_grid] = ndgrid(snr_dB_rng,p_fa_1);
% [snr_dB_grid2, p_fa2_grid] = ndgrid(snr_dB_rng,p_fa_2);

p_d_1_ideal = zeros(N_snr, 2*N_fa, N_seq_rng_numel, N_cases, N_freq);
p_d_1_mAv = zeros(N_snr, 2*N_fa, N_seq_rng_numel, N_cases, N_freq, N_L);
p_d_1_mMed = zeros(N_snr, 2*N_fa, N_seq_rng_numel, N_cases, N_freq, N_L);
p_d_1_nEst = zeros(N_snr, 2*N_fa, N_seq_rng_numel, N_cases, N_freq, N_L);

p_d_2_ideal = zeros(N_snr, 2*N_fa, N_seq_rng_numel, N_cases, N_freq);
p_d_2_mAv = zeros(N_snr, 2*N_fa, N_seq_rng_numel, N_cases, N_freq, N_L);
p_d_2_mMed = zeros(N_snr, 2*N_fa, N_seq_rng_numel, N_cases, N_freq, N_L);
p_d_2_nEst = zeros(N_snr, 2*N_fa, N_seq_rng_numel, N_cases, N_freq, N_L);

fprintf(sprintf(['start: ', char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss')), '\n']))
N_frames = 1e1;
N_rep = 2;
for iCase = 1:N_cases
    h = h_preambs{iCase};
    tStart = tic;
    for iN_seq = 1:N_seq_rng_numel
        for iFreq = 1:N_freq
            p_d_sim = cell2mat(arrayfun(@( p_fa) ...
                p_d_y1_ideal(snr_dB_rng, p_fa, h, N_seq_rng(iN_seq), freq_rng(iFreq), N_frames, N_rep), ...
                p_fa_rng.','UniformOutput',false)); %(snr_dB, N_frames, h0, p_fa, N_reph0, N_rep)
            p_d_1_ideal(:,:,iN_seq, iCase, iFreq) = p_d_sim;

            p_d_sim = cell2mat(arrayfun(@( p_fa) ...
                p_d_y2_ideal(snr_dB_rng, p_fa, h, N_seq_rng(iN_seq), freq_rng(iFreq), N_frames, N_rep), ...
                p_fa_rng.','UniformOutput',false)); %(snr_dB, N_frames, h0, p_fa, N_reph0, N_rep)
            p_d_2_ideal(:,:,iN_seq, iCase, iFreq) = p_d_sim;
            
            for iL = 1:N_L
                p_d_sim = cell2mat(arrayfun(@( p_fa) ...
                    p_d_y1_mAv(snr_dB_rng, p_fa, h, N_seq_rng(iN_seq), freq_rng(iFreq), N_frames, N_rep, L_rng(iL)), ...
                    p_fa_rng.','UniformOutput',false)); %(snr_dB, N_frames, h0, p_fa, N_reph0, N_rep)
                p_d_1_mAv(:,:,iN_seq, iCase, iFreq, iL) = p_d_sim;
                
                p_d_sim = cell2mat(arrayfun(@( p_fa) ...
                    p_d_y1_mMed(snr_dB_rng, p_fa, h, N_seq_rng(iN_seq), freq_rng(iFreq), N_frames, N_rep, L_rng(iL)), ...
                    p_fa_rng.','UniformOutput',false)); %(snr_dB, N_frames, h0, p_fa, N_reph0, N_rep)
                p_d_1_mMed(:,:,iN_seq, iCase, iFreq, iL) = p_d_sim;

                p_d_sim = cell2mat(arrayfun(@( p_fa) ...
                    p_d_y1_nEst(snr_dB_rng, p_fa, h, N_seq_rng(iN_seq), freq_rng(iFreq), N_frames, N_rep, L_rng(iL)), ...
                    p_fa_rng.','UniformOutput',false)); %(snr_dB, N_frames, h0, p_fa, N_reph0, N_rep)
                p_d_1_nEst(:,:,iN_seq, iCase, iFreq, iL) = p_d_sim;



                p_d_sim = cell2mat(arrayfun(@( p_fa) ...
                    p_d_y2_mAv(snr_dB_rng, p_fa, h, N_seq_rng(iN_seq), freq_rng(iFreq), N_frames, N_rep, L_rng(iL)), ...
                    p_fa_rng.','UniformOutput',false)); %(snr_dB, N_frames, h0, p_fa, N_reph0, N_rep)
                p_d_2_mAv(:,:,iN_seq, iCase, iFreq, iL) = p_d_sim;
                
                p_d_sim = cell2mat(arrayfun(@( p_fa) ...
                    p_d_y2_mMed(snr_dB_rng, p_fa, h, N_seq_rng(iN_seq), freq_rng(iFreq), N_frames, N_rep, L_rng(iL)), ...
                    p_fa_rng.','UniformOutput',false)); %(snr_dB, N_frames, h0, p_fa, N_reph0, N_rep)
                p_d_2_mMed(:,:,iN_seq, iCase, iFreq, iL) = p_d_sim;

                p_d_sim = cell2mat(arrayfun(@( p_fa) ...
                    p_d_y2_nEst(snr_dB_rng, p_fa, h, N_seq_rng(iN_seq), freq_rng(iFreq), N_frames, N_rep, L_rng(iL)), ...
                    p_fa_rng.','UniformOutput',false)); %(snr_dB, N_frames, h0, p_fa, N_reph0, N_rep)
                p_d_2_nEst(:,:,iN_seq, iCase, iFreq, iL) = p_d_sim;
            end
        end
    end
    tElapsed = toc(tStart);
    fprintf('iCase: %d, Elapsed time: %.3f\n', iCase, tElapsed);
end

if ~exist(outDir,'dir')
    mkdir(outDir);
end
fOutFile = fullfile(outDir,...
    sprintf([char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss')),...
    '_p_d_vs_snr_freqOff_hGeneral_N_framesN_rep%d.mat'],N_frames*N_rep));

save(fOutFile, ...
    'p_d_1_ideal', 'p_d_1_mAv', 'p_d_1_mMed', 'p_d_1_nEst',...
    'p_d_2_ideal', 'p_d_2_mAv', 'p_d_2_mMed', 'p_d_2_nEst',...
    'h_preambs_params', 'h_preambs', 'N_cases',...
    'snr_dB_rng', 'N_snr', 'p_fa_rng', 'N_fa', 'L_rng', 'N_L',...
    'N_seq_rng', 'N_seq_rng_numel', 'freq_rng', 'N_freq',...
    'N_frames', 'N_rep');
 fprintf(sprintf(['end: ', char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss')), '\n']))