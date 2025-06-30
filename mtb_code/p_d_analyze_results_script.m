outDir = 'data_out_probs';
% fName_in = '2025-06-08_21-17-54_p_d_vs_snr_freqOff_hGeneral_N_framesN_rep12000.mat';
fName_in = '2025-06-15_04-04-51_p_d_vs_snr_freqOff_hGeneral_N_framesN_rep120000.mat'; %only ideal case
load(fullfile('..', outDir, fName_in));

% colors = ...
% [240,249,232;
% 186,228,188;
% 123,204,196;
% 67,162,202;
% 8,104,172]./255;

% colors = ...
% [27,158,119;
% 217,95,2;
% 117,112,179;
% 231,41,138;
% 102,166,30;
% 100,100,100]./255;

% colors = ...
% [166,206,227;
% 31,120,180;
% 178,223,138;
% 51,160,44;
% 251,154,153;
% 227,26,28;
% 253,191,111;
% 255,127,0;
% 202,178,214;
% 106,61,154]./255;

% colors = ...
% [166,206,227;
% 31,120,180;
% 178,223,138;
% 51,160,44;
% 251,154,153;
% 227,26,28;
% 253,191,111;
% 255,127,0;
% 202,178,214;
% 106,61,154;
% 255,255,153]./255;

colors = ...
[166,206,227;
31,120,180;
178,223,138;
51,160,44;
251,154,153;
227,26,28;
253,191,111;
255,127,0;
202,178,214;
106,61,154;
0,0,0;
255,0,0]./255;


lnStyle = {'-','-.',':','--'};


ltx  = {'Interpreter','LaTeX'};
%% check the ideal case, compare it with the analytical expression
lnWid = {'LineWidth',1.7};
% semilogy(snr_dB_rng, 1-p_d_1_ideal(:,1, 1, 1))
p_d_midVal_vs_true = 1; %1 or 2
p_fa_id = 2; %
legCell = arrayfun(@(x) sprintf('$f_{\\varepsilon}$=%.3f', x), freq_rng,...
    'UniformOutput', false);
xlims = [min(snr_dB_rng), max(snr_dB_rng)];
ylims = [1/(1.2e5), 1];

xlbl = 'SNR [dB]';
ylbl = '$p_{\rm md} [-]$';

p_fa = p_fa_rng(p_fa_id);
p_fa_id_ = (p_fa_id - 1)*2 + p_d_midVal_vs_true;
N_snr_an = 40;
snr_dB_rng_an = linspace(min(snr_dB_rng), max(snr_dB_rng), N_snr_an).';
for iCase = 1:N_cases
    h = h_preambs{iCase};
    for iN_seq = 1:N_seq_rng_numel
        subplot(2,2,iN_seq);
        % *** analytical y'_1 ***
%         for iFreq = 1:N_freq
% %             p_md = p_md_y1_ideal_analytical_hGeneral(snr_dB_rng_an, h, ...
% %                 N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
% %             p_md = p_md1(snr_dB_rng_an, h, ...
% %                 N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
%             [p_md,~] = get_p_md(snr_dB_rng_an, h, ...
%                 N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
%             semilogy(snr_dB_rng_an, p_md, 'Color', colors(iFreq,:),...
%                     lnWid{:}, 'LineStyle','-');
%             if iFreq == 1
%                 hold on
%             end
%         end
        % *** analytical y'_2 ***
        for iFreq = 1:N_freq
%             p_md = p_md_y2_ideal_analytical_hGeneral(snr_dB_rng_an, h, ...
%                 N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
%             p_md = p_md2(snr_dB_rng_an, h, ...
%                 N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
            [~, p_md] = get_p_md(snr_dB_rng_an, h, ...
                N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
            semilogy(snr_dB_rng_an, p_md, 'Color', colors(iFreq,:),...
                    lnWid{:}, 'LineStyle','--');
            if iFreq == 1
                hold on
            end
        end

        max_nonZero_val_id = 1;
%         for iFreq = 1:N_freq
%             % *** sim y'_1 ***
%             p_md = 1 - p_d_1_ideal(:, p_fa_id_, iN_seq, iCase, iFreq);
%             last_nz_id = find(p_md,1,'last');
%             if last_nz_id > max_nonZero_val_id 
%                 max_nonZero_val_id = last_nz_id;
%             end
%             semilogy(snr_dB_rng, p_md, 'Color', colors(iFreq,:),...
%                     lnWid{:}, 'LineStyle','-.', 'Marker','x');
%         end
        % *** sim y'_2 ***
        for iFreq = 1:N_freq
            p_md = 1 - p_d_2_ideal(:, p_fa_id_, iN_seq, iCase, iFreq);
            last_nz_id = find(p_md,1,'last');
            if last_nz_id > max_nonZero_val_id 
                max_nonZero_val_id = last_nz_id;
            end
            semilogy(snr_dB_rng, p_md, 'Color', colors(iFreq,:),...
                    lnWid{:}, 'LineStyle',':', 'Marker','o');
        end
        hold off

        grid on
        grid minor

%         xlim([xlims(1), snr_dB_rng(max_nonZero_val_id)]);
        xlim(xlims);
        ylim(ylims);
        xlabel(xlbl, ltx{:});
        ylabel(ylbl, ltx{:});
%             if idxNrep == 1 && iFreq == 11
        legend(legCell{:},ltx{:}, 'Location','southwest');
%             end
        
        tit_str = sprintf(['$p_{\\rm fa}=10^{%d}$, $N_{\\rm seq}=%d$,',...
            'iCase: %02d'], round(log10(p_fa)), N_seq_rng(iN_seq),...
            iCase);
        title(tit_str, ltx{:});
    end
    pause
end

%% compare ideal (analytical) case and mAv, mMed, nEst cases
fName_in = '2025-06-08_21-17-54_p_d_vs_snr_freqOff_hGeneral_N_framesN_rep12000.mat';
load(fullfile('..', outDir, fName_in));

p_d_midVal_vs_true = 1; %1 or 2
p_fa_id = 2; %
xlims = [min(snr_dB_rng), max(snr_dB_rng)];
ylims = [1/(1.2e4), 1];

xlbl = 'SNR [dB]';
ylbl = '$p_{\rm md} [-]$';

p_fa = p_fa_rng(p_fa_id);
p_fa_id_ = (p_fa_id - 1)*2 + p_d_midVal_vs_true;
N_snr_an = 40;
snr_dB_rng_an = linspace(min(snr_dB_rng), max(snr_dB_rng), N_snr_an).';
freq_rng_id = [1, 6, 11];

legCell = arrayfun(@(x) sprintf('$f_{\\varepsilon}$=%.3f', x),...
    freq_rng(freq_rng_id),...
    'UniformOutput', false);
one_or_two = 1;
if one_or_two == 1
    p_md_an_fcn = @p_md_y1_ideal_analytical_hGeneral;
    p_d_ideal = p_d_1_ideal;
    p_d_mAv = p_d_1_mAv;
    p_d_mMed = p_d_1_mMed;
    p_d_nEst = p_d_1_nEst;
else
    p_d_ideal = p_d_2_ideal;
    p_d_mAv = p_d_2_mAv;
    p_d_mMed = p_d_2_mMed;
    p_d_nEst = p_d_2_nEst;
end
phi_norm = @(x) 1/2*(1 + erf(x/sqrt(2)));
cdf_rNorm = @(x, r, sigma2) phi_norm((x - r)./sqrt(sigma2));
iL = 4;
% p_d_2_mAv(:,:,iN_seq, iCase, iFreq, iL)
for iCase = 1:N_cases
    h = h_preambs{iCase};
    for iN_seq = 1:N_seq_rng_numel
        subplot(2, 2, iN_seq);
        last_nz_id = 1;
        % *** analytical y'_1 ***
        for iFreq = freq_rng_id %1:N_freq
            [p_md,~] = get_p_md(snr_dB_rng_an, h, ...
                N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
            semilogy(snr_dB_rng_an, p_md, 'Color', colors(iFreq,:),...
                    lnWid{:}, 'LineStyle','-');
            if iFreq == 1
                hold on
            end
        end

         % *** analytical y'_2_mAv ***
        for iFreq = freq_rng_id %1:N_freq
            [~,p_md] = get_p_md(snr_dB_rng_an, h, ...
                N_seq_rng(iN_seq), p_fa, freq_rng(iFreq), L_rng(iL));
            semilogy(snr_dB_rng_an, p_md, 'Color', colors(iFreq,:),...
                    lnWid{:}, 'LineStyle',':');
            if iFreq == 1
                hold on
            end
        end
        
        % mAv
        for iFreq = freq_rng_id %1:N_freq
            p_md = 1 - p_d_mAv(:, p_fa_id_, iN_seq, iCase, iFreq, iL);
            last_nz_id = find(p_md,1,'last');
            if last_nz_id > max_nonZero_val_id 
                max_nonZero_val_id = last_nz_id;
            end
            semilogy(snr_dB_rng, p_md, 'Color', colors(iFreq,:),...
                    lnWid{:}, 'LineStyle','--', 'Marker','x');
        end

        % mMed
        for iFreq = freq_rng_id %1:N_freq
            p_md = 1 - p_d_mMed(:, p_fa_id_, iN_seq, iCase, iFreq, iL);
            last_nz_id = find(p_md,1,'last');
            if last_nz_id > max_nonZero_val_id 
                max_nonZero_val_id = last_nz_id;
            end
            semilogy(snr_dB_rng, p_md, 'Color', colors(iFreq,:),...
                    lnWid{:}, 'LineStyle','-.', 'Marker','x');
        end

        % nEst
        for iFreq = freq_rng_id %1:N_freq
            p_md = 1 - p_d_nEst(:, p_fa_id_, iN_seq, iCase, iFreq, iL);
            last_nz_id = find(p_md,1,'last');
            if last_nz_id > max_nonZero_val_id 
                max_nonZero_val_id = last_nz_id;
            end
            semilogy(snr_dB_rng, p_md, 'Color', colors(iFreq,:),...
                    lnWid{:}, 'LineStyle',':', 'Marker','x');
        end

        hold off

        grid on
        grid minor

        xlim([xlims(1), snr_dB_rng(max_nonZero_val_id)]);
        ylim(ylims);
        xlabel(xlbl, ltx{:});
        ylabel(ylbl, ltx{:});
%             if idxNrep == 1 && iFreq == 11
        legend(legCell{:},ltx{:}, 'Location','southwest');
%             end
        
        tit_str = sprintf(['$p_{\\rm fa}=10^{%d}$, $N_{\\rm seq}=%d$,',...
            'iCase: %02d'], round(log10(p_fa)), N_seq_rng(iN_seq),...
            iCase);
        title(tit_str, ltx{:});
    end
    pause
end

%% check the ideal case, compare it with the analytical expression, 
% p_d as a function of F_eps
fName_in = '2025-06-15_04-04-51_p_d_vs_snr_freqOff_hGeneral_N_framesN_rep120000.mat'; %only ideal case
load(fullfile('..', outDir, fName_in));

lnWid = {'LineWidth',1.7};
p_d_midVal_vs_true = 1; %1 or 2
p_fa_id = 2; %
% legCell = arrayfun(@(x) sprintf('$f_{\\varepsilon}$=%.3f', x), freq_rng,...
%     'UniformOutput', false);
snr_dB_idc = [6, 11, 16];
N_snr_idc = numel(snr_dB_idc);
legCell = arrayfun(@(x) sprintf('SNR=$%.3f$ dB', x), snr_dB_rng(snr_dB_idc),...
    'UniformOutput', false);
xlims = [min(freq_rng), max(freq_rng)];
ylims = [1/(1.2e5), 1];

xlbl = '$F_{\varepsilon}$ [-]';
ylbl = '$p_{\rm md}$ [-]';

p_fa = p_fa_rng(p_fa_id);
p_fa_id_ = (p_fa_id - 1)*2 + p_d_midVal_vs_true;
N_freq_an = 100;
freq_rng_an = linspace(min(freq_rng), max(freq_rng), N_freq_an).';
for iCase = 1:N_cases
    h = h_preambs{iCase};
    for iN_seq = 1:N_seq_rng_numel
        subplot(2,2,iN_seq);
        % *** analytical y'_1 ***
        p_md = zeros(N_snr_idc, N_freq);
        for iFreq = 1:N_freq_an
        [p_md(:,iFreq), ~] = get_p_md(snr_dB_rng(snr_dB_idc), h, ...
            N_seq_rng(iN_seq), p_fa, freq_rng_an(iFreq));
        end
        for iSnr = 1:N_snr_idc
            semilogy(freq_rng_an, p_md(iSnr,:), 'Color', colors(iSnr,:),...
                    lnWid{:}, 'LineStyle',':');
            if iSnr == 1
                hold on
            end
        end

        % *** analytical y'_2 ***
%         p_md = zeros(N_snr_idc, N_freq);
%         for iFreq = 1:N_freq_an
%         [~, p_md(:,iFreq)] = get_p_md(snr_dB_rng(snr_dB_idc), h, ...
%             N_seq_rng(iN_seq), p_fa, freq_rng_an(iFreq));
%         end
%         for iSnr = 1:N_snr_idc
%             semilogy(freq_rng_an, p_md(iSnr,:), 'Color', colors(iSnr,:),...
%                     lnWid{:}, 'LineStyle',':');
%             if iSnr == 1
%                 hold on
%             end
%         end

        max_nonZero_val_id = 1;

        p_md = zeros(N_snr_idc, N_freq);
        for iFreq = 1:N_freq
            p_md_aux = 1 - p_d_1_ideal(snr_dB_idc, p_fa_id_, iN_seq, iCase, iFreq);
            p_md(:,iFreq) = p_md_aux;
            last_nz_id = find(p_md_aux,1,'last');
            if last_nz_id > max_nonZero_val_id 
                max_nonZero_val_id = last_nz_id;
            end
        end
        for iSnr = 1:N_snr_idc
            semilogy(freq_rng, p_md(iSnr,:), 'Color', colors(iSnr,:),...
                    lnWid{:}, 'LineStyle', '--', 'Marker','o');
        end

        % *** sim y'_2 ***
%         p_md = zeros(N_snr_idc, N_freq);
%         for iFreq = 1:N_freq
%             p_md_aux = 1 - p_d_2_ideal(snr_dB_idc, p_fa_id_, iN_seq, iCase, iFreq);
%             p_md(:,iFreq) = p_md_aux;
%             last_nz_id = find(p_md_aux,1,'last');
%             if last_nz_id > max_nonZero_val_id 
%                 max_nonZero_val_id = last_nz_id;
%             end
%         end
%         for iSnr = 1:N_snr_idc
%             semilogy(freq_rng, p_md(iSnr,:), 'Color', colors(iSnr,:),...
%                     lnWid{:}, 'LineStyle', '--', 'Marker','o');
%         end
        hold off

        grid on
        grid minor

%         xlim([xlims(1), snr_dB_rng(max_nonZero_val_id)]);
        xlim(xlims);
        ylim(ylims);
        xlabel(xlbl, ltx{:});
        ylabel(ylbl, ltx{:});
%             if idxNrep == 1 && iFreq == 11
        legend(legCell{:},ltx{:}, 'Location','southwest');
%             end
        
        tit_str = sprintf(['$p_{\\rm fa}=10^{%d}$, $N_{\\rm seq}=%d$,',...
            'iCase: %02d'], round(log10(p_fa)), N_seq_rng(iN_seq),...
            iCase);
        title(tit_str, ltx{:});
    end
    pause
end

%% V2. Check the ideal case -- the analytical expression
% p_d as a function of F_eps
lnWid = {'LineWidth',1.7};
p_fa_id = 2; %
snr_dB_rng_aux = [15]; %, 22, 30
N_snr_aux = numel(snr_dB_rng_aux);
[N_seq_grid_aux, snr_grid_aux,] = ndgrid(N_seq_rng, snr_dB_rng_aux);
legCell = arrayfun(@(N_seq_, snr_) ...
    sprintf('$N_{\\rm seq}=%d$, SNR=$%.1f$ dB', N_seq_, snr_), ...
     N_seq_grid_aux, snr_grid_aux, ...
    'UniformOutput', false);
xlbl = '$F_{\varepsilon}$ [-]';
ylbl = '$p_{\rm d}$ [-]';
N_freq_an = 800;
freq_rng_an = linspace(0, 0.105, N_freq_an).';

clf

% xlims = [min(freq_rng), max(freq_rng)];
xlims = [min(freq_rng_an), max(freq_rng_an)];
ylims = [1/(1.2e5), 1];
% ylims = [3e-7, 1];

p_fa = p_fa_rng(p_fa_id);
% p_fa_id_ = (p_fa_id - 1)*2 + p_d_midVal_vs_true;
for iCase = 1:N_cases
    h = h_preambs{iCase};
    for iSnr = 1:N_snr_aux
        % *** analytical y'_2 ***
        p_md = zeros(N_seq_rng_numel, N_freq);
        for iFreq = 1:N_freq_an
            for iN_seq = 1:N_seq_rng_numel
                [~, p_md(iN_seq,iFreq)] = get_p_md(snr_dB_rng_aux(iSnr), h, ...
                N_seq_rng(iN_seq), p_fa, freq_rng_an(iFreq));
            end
        end
        for iN_seq = 1:N_seq_rng_numel
            iColor = (iSnr-1)*N_seq_rng_numel + iN_seq;
            semilogy(freq_rng_an, 1-p_md(iN_seq,:), 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle','--');
            if iSnr == 1 && iN_seq == 1
                hold on
            end
        end
    end

%     for iSnr = 1:N_snr_aux
%         % *** analytical y'_1 ***
%         p_md = zeros(N_seq_rng_numel, N_freq);
%         for iFreq = 1:N_freq_an
%             for iN_seq = 1:N_seq_rng_numel
%                 [p_md(iN_seq,iFreq), ~] = get_p_md(snr_dB_rng_aux(iSnr), h, ...
%                 N_seq_rng(iN_seq), p_fa, freq_rng_an(iFreq));
%             end
%         end
%         for iN_seq = 1:N_seq_rng_numel
%             iColor = (iSnr-1)*N_seq_rng_numel + iN_seq;
%             semilogy(freq_rng_an, 1-p_md(iN_seq,:), 'Color', colors(iColor,:),...
%                     lnWid{:}, 'LineStyle',':');
%             if iSnr == 1 && iN_seq == 1
%                 hold on
%             end
%         end      
%     end

    
    % write auxiliary lines -- the lobe border
    N_h = numel(h);
    for iN_seq = 1:N_seq_rng_numel
        x_coord = repmat( N_seq_rng(iN_seq)/N_h, 2,1);
        y_coord = [3e-7, 1];
        semilogy(x_coord, y_coord, 'Color', colors(iN_seq,:),...
                'LineWidth',1.0, 'LineStyle','-');
    end

    hold off

    grid on
    grid minor

%         xlim([xlims(1), snr_dB_rng(max_nonZero_val_id)]);
    xlim(xlims);
    ylim(ylims);
    xlabel(xlbl, ltx{:});
    ylabel(ylbl, ltx{:});
%             if idxNrep == 1 && iFreq == 11
    legend(legCell{:},ltx{:}, 'Location','southwest');
%             end
    
    tit_str = sprintf(['$p_{\\rm fa}=10^{%d}$,',...
        'iCase: %02d'], round(log10(p_fa)), ...
        iCase);
    title(tit_str, ltx{:});
    

    pause
end