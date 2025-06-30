outDir = 'data_out_probs';
% fName_in = '2025-06-08_21-17-54_p_d_vs_snr_freqOff_hGeneral_N_framesN_rep12000.mat';
fName_in = '2025-06-15_04-04-51_p_d_vs_snr_freqOff_hGeneral_N_framesN_rep120000.mat'; %only ideal case
load(fullfile('..', outDir, fName_in));


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
% 0,0,0]./255;

% 6 colors, https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=6
colors = [...
166,206,227;
31,120,180;
178,223,138;
51,160,44;
251,154,153;
227,26,28]./255;


lnStyle = {'-','-.',':','--'};

xlbl = 'SNR [dB]';
ylbl = '$p_{\rm md} [-]$';
ltx  = {'Interpreter','LaTeX'};
%% compare ideal case, sim vs analytical formula
lnWid = {'LineWidth',1.7};
fntSz_leg = 11;
fntSz_axLabel = 14;
fntSz_title = 14;
% semilogy(snr_dB_rng, 1-p_d_1_ideal(:,1, 1, 1))
p_d_midVal_vs_true = 2; %1 or 2
p_fa_id = 2; %
xlims = [min(snr_dB_rng), max(snr_dB_rng)];
ylims = [1/(1.2e5), 1];

p_fa = p_fa_rng(p_fa_id);
p_fa_id_ = (p_fa_id - 1)*2 + p_d_midVal_vs_true;
N_snr_an = 300;
snr_dB_rng_an = linspace(min(snr_dB_rng), max(snr_dB_rng), N_snr_an).';
%
% prmbID_Nseq__cases = ...
%     [1, 1;
%      1, 4;
%      2, 2;
%      7, 3;
%      8, 2];
prmbID_Nseq__cases = ...
    [1, 4;
     7, 3;
     8, 2];
N_pair_cases = size(prmbID_Nseq__cases,1);
cases_mtx = getAllTuples(1:N_pair_cases, [3, 10]);
N_plotCases = size(cases_mtx,1);
p_md_id = 2;


legCell = arrayfun(@(id1, id_Freq) ...
    sprintf('id:$%02d$, $N_{\\rm seq}=%02d$, $F_{\\varepsilon}=%.3f$', ...
    prmbID_Nseq__cases(id1,1),...
    N_seq_rng(prmbID_Nseq__cases(id1,2)),...
    freq_rng(id_Freq)), ...
    cases_mtx(:,1), cases_mtx(:,2), ...
    'UniformOutput', false);
% plot simulation results
for iPlotCase = 1:N_plotCases
    prmb_id = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),1);
    iN_seq = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),2);
    iFreq = cases_mtx(iPlotCase,2);

    h = h_preambs{ prmb_id };
    N_seq = N_seq_rng( iN_seq );
    F_eps = freq_rng( iFreq);

%     figure(1);
    if p_md_id == 1
    p_md = 1 - p_d_1_ideal(:, p_fa_id_, iN_seq, prmb_id, iFreq);
    else
    p_md = 1 - p_d_2_ideal(:, p_fa_id_, iN_seq, prmb_id, iFreq);
    end
    semilogy(snr_dB_rng, p_md, 'Color', colors(iPlotCase,:),...
                    lnWid{:}, 'LineStyle','-', 'Marker','x');
    if iPlotCase == 1
        hold on
    end
end
% plot analytical results
for iPlotCase = 1:N_plotCases
    prmb_id = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),1);
    iN_seq = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),2);
    iFreq = cases_mtx(iPlotCase,2);

    h = h_preambs{ prmb_id };

%     figure(1);
    [p_md1, p_md2] = get_p_md(snr_dB_rng_an, h, ...
        N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
    if p_md_id == 1
        p_md = p_md1;
    else
        p_md = p_md2;
    end
    semilogy(snr_dB_rng_an, p_md, 'Color', colors(iPlotCase,:),...
            lnWid{:}, 'LineStyle',':');
%     if iPlotCase == 1
%         hold on
%     end
end

hold off

grid on
grid minor

set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',fntSz_axLabel-2);
%         xlim([xlims(1), snr_dB_rng(max_nonZero_val_id)]);
xlim(xlims);
ylim(ylims);
xlabel(xlbl, ltx{:}, 'FontSize',fntSz_axLabel);
ylabel(sprintf('$p_{{\\rm md},%d}$ [-]', p_md_id), ltx{:}, 'FontSize', fntSz_axLabel);
if p_d_midVal_vs_true == 2
    ylabel(sprintf('$p_{{\\rm md}}$ [-]'), ltx{:});
end
%             if idxNrep == 1 && iFreq == 11
legend(legCell{:},ltx{:}, 'Location','southeast','FontSize',fntSz_leg);
%             end

tit_str = sprintf('$p_{\\rm fa}=10^{%d}$, $p_{{\\rm md},%d}$ -- analytical formula vs. simulation',...
    round(log10(p_fa)), p_md_id );
if p_d_midVal_vs_true == 2
tit_str = sprintf('$p_{\\rm fa}=10^{%d}$, $p_{{\\rm md},%d}=1-\\Pr\\{{\\cal E}_{%d,{\\rm approx}}\\}$ (an.) vs. $1-\\Pr\\{{\\cal E}_{%d}\\}$ (sim.)',...
    round(log10(p_fa)), p_md_id, p_md_id, p_md_id );
end
title(tit_str, ltx{:}, 'FontSize', fntSz_title);

fName_out = fullfile('..', 'text', 'radio_eng_article', 'figs',...
    sprintf('p_md_%d_cases.pdf', p_md_id) );
if p_d_midVal_vs_true == 2
fName_out = fullfile('..', 'text', 'radio_eng_article', 'figs',...
    sprintf('p_md_%d_cases_truePmd.pdf', p_md_id) );
end
% export_fig(fName_out,'-pdf','-transparent');
%% compare pmd_1 and p_md2, analytical formulas
% plot analytical results p_md2
for iPlotCase = 1:N_plotCases
    prmb_id = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),1);
    iN_seq = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),2);
    iFreq = cases_mtx(iPlotCase,2);

    h = h_preambs{ prmb_id };

    [~, p_md2] = get_p_md(snr_dB_rng_an, h, ...
        N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
    p_md = p_md2;

    semilogy(snr_dB_rng_an, p_md, 'Color', colors(iPlotCase,:),...
            lnWid{:}, 'LineStyle','--');
    if iPlotCase == 1
        hold on
    end
end
% plot analytical results p_md1
for iPlotCase = 1:N_plotCases
    prmb_id = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),1);
    iN_seq = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),2);
    iFreq = cases_mtx(iPlotCase,2);

    h = h_preambs{ prmb_id };
    N_seq = N_seq_rng( iN_seq );
    F_eps = freq_rng( iFreq);

    [ p_md1, ~] = get_p_md(snr_dB_rng_an, h, ...
        N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
    p_md = p_md1;
    semilogy(snr_dB_rng_an, p_md, 'Color', colors(iPlotCase,:),...
                    lnWid{:}, 'LineStyle',':');
%     if iPlotCase == 1
%         hold on
%     end
end


hold off

grid on
grid minor

set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',fntSz_axLabel-2);
%         xlim([xlims(1), snr_dB_rng(max_nonZero_val_id)]);
xlim(xlims);
ylim(ylims);
xlabel(xlbl, ltx{:}, 'FontSize', fntSz_axLabel);
ylabel(sprintf('$p_{{\\rm md}}$ [-]'), ltx{:}, 'FontSize', fntSz_axLabel);
%             if idxNrep == 1 && iFreq == 11
legend(legCell{:},ltx{:}, 'Location','southeast','FontSize',fntSz_leg);
%             end

tit_str = sprintf('$p_{\\rm fa}=10^{%d}$, $p_{{\\rm md},1}$ (analytical) vs $p_{{\\rm md},2}$  (analytical)',...
    round(log10(p_fa)) );
title(tit_str, ltx{:}, 'FontSize', fntSz_title);

fName_out = fullfile('..', 'text', 'radio_eng_article', 'figs',...
    sprintf('p_md_1_vs_2_analytical_cases.pdf') );
% export_fig(fName_out,'-pdf','-transparent');

%% compare ideal-case, and nEst, mAv, mMed variant
fName_in = '2025-06-08_21-17-54_p_d_vs_snr_freqOff_hGeneral_N_framesN_rep12000.mat';
load(fullfile('..', outDir, fName_in));
lnWid = {'LineWidth',1.7};
% fntSz_leg = 8;
% semilogy(snr_dB_rng, 1-p_d_1_ideal(:,1, 1, 1))
p_d_midVal_vs_true = 1; %1 or 2
p_fa_id = 2; %
xlims = [min(snr_dB_rng), max(snr_dB_rng)];
ylims = [1/(1.2e4), 1];

p_fa = p_fa_rng(p_fa_id);
p_fa_id_ = (p_fa_id - 1)*2 + p_d_midVal_vs_true;
N_snr_an = 300;
snr_dB_rng_an = linspace(min(snr_dB_rng), max(snr_dB_rng), N_snr_an).';
%
% prmbID_Nseq__cases = ...
%     [1, 1;
%      1, 4;
%      2, 2;
%      7, 3;
%      8, 2];
prmbID_Nseq__cases = ...
    [1, 4;
     7, 3;
     8, 2];
N_pair_cases = size(prmbID_Nseq__cases,1);
cases_mtx = getAllTuples(1:N_pair_cases, [3, 10]);
N_plotCases = size(cases_mtx,1);

choose_ids = [2, 3];
choose_ids = reshape([2*choose_ids-1; 2*choose_ids-1+1], [],1);
N_chosen = numel(choose_ids);
p_md_id = 2;
iL = 3;

legCell = arrayfun(@(id1, id_Freq) ...
    sprintf('id:$%02d$, $N_{\\rm seq}=%02d$, $F_{\\varepsilon}=%.3f$', ...
    prmbID_Nseq__cases(id1,1),...
    N_seq_rng(prmbID_Nseq__cases(id1,2)),...
    freq_rng(id_Freq)), ...
    cases_mtx(:,1), cases_mtx(:,2), ...
    'UniformOutput', false);
% plot simulation results, nEst
for id = 1:numel(choose_ids)
    iColor = id + 2;
    iPlotCase = choose_ids(id);
    prmb_id = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),1);
    iN_seq = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),2);
    iFreq = cases_mtx(iPlotCase,2);

    h = h_preambs{ prmb_id };
    N_seq = N_seq_rng( iN_seq );
    F_eps = freq_rng( iFreq);

    if p_md_id == 1
    p_md = 1 - p_d_1_nEst(:, p_fa_id_, iN_seq, prmb_id, iFreq, iL);
    else
    p_md = 1 - p_d_2_nEst(:, p_fa_id_, iN_seq, prmb_id, iFreq, iL);
    end
    semilogy(snr_dB_rng, p_md, 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle','-', 'Marker','x');
    if id == 1
        hold on
    end
end
% plot analytical results
for id = 1:numel(choose_ids)
    iColor = id + 2;
    iPlotCase = choose_ids(id);
    prmb_id = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),1);
    iN_seq = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),2);
    iFreq = cases_mtx(iPlotCase,2);

    h = h_preambs{ prmb_id };

%     figure(1);
    [p_md1, p_md2] = get_p_md(snr_dB_rng_an, h, ...
        N_seq_rng(iN_seq), p_fa, freq_rng(iFreq));
    if p_md_id == 1
        p_md = p_md1;
    else
        p_md = p_md2;
    end
    semilogy(snr_dB_rng_an, p_md, 'Color', colors(iColor,:),...
            lnWid{:}, 'LineStyle',':');
%     if iPlotCase == 1
%         hold on
%     end
end
% plot simulation results, mAv
for id = 1:numel(choose_ids)
    iColor = id + 2;
    iPlotCase = choose_ids(id);
    prmb_id = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),1);
    iN_seq = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),2);
    iFreq = cases_mtx(iPlotCase,2);

    h = h_preambs{ prmb_id };
    N_seq = N_seq_rng( iN_seq );
    F_eps = freq_rng( iFreq);

    if p_md_id == 1
    p_md = 1 - p_d_1_mAv(:, p_fa_id_, iN_seq, prmb_id, iFreq, iL);
    else
    p_md = 1 - p_d_2_mAv(:, p_fa_id_, iN_seq, prmb_id, iFreq, iL);
    end
    semilogy(snr_dB_rng, p_md, 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle','--', 'Marker','x');
end
% plot simulation results, mMed
for id = 1:numel(choose_ids)
    iColor = id + 2;
    iPlotCase = choose_ids(id);
    prmb_id = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),1);
    iN_seq = prmbID_Nseq__cases(cases_mtx(iPlotCase,1),2);
    iFreq = cases_mtx(iPlotCase,2);

    h = h_preambs{ prmb_id };
    N_seq = N_seq_rng( iN_seq );
    F_eps = freq_rng( iFreq);

    if p_md_id == 1
    p_md = 1 - p_d_1_mMed(:, p_fa_id_, iN_seq, prmb_id, iFreq, iL);
    else
    p_md = 1 - p_d_2_mMed(:, p_fa_id_, iN_seq, prmb_id, iFreq, iL);
    end
    semilogy(snr_dB_rng, p_md, 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle','-.', 'Marker','x');
end


hold off

grid on
grid minor

set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',fntSz_axLabel-2);
%         xlim([xlims(1), snr_dB_rng(max_nonZero_val_id)]);
xlim(xlims);
ylim(ylims);
xlabel(xlbl, ltx{:}, 'FontSize', fntSz_axLabel);
% ylabel(sprintf('$p_{{\\rm md},%d}$ [-]', p_md_id), ltx{:});
% if p_d_midVal_vs_true == 2
    ylabel(sprintf('$p_{{\\rm md}}$ [-]'), ltx{:}, 'FontSize', fntSz_axLabel);
% end
%             if idxNrep == 1 && iFreq == 11
legend(legCell{choose_ids},ltx{:}, 'Location','southeast','FontSize',fntSz_leg);
%             end

tit_str = sprintf(['$p_{\\rm fa}=10^{%d}$, $p_{{\\rm md},%d}$ ',...
    '(an.) vs. nEst vs. mAv vs. mMed with $L=%d$'],...
    round(log10(p_fa)), p_md_id, L_rng(iL) );
if p_d_midVal_vs_true == 2
tit_str = sprintf(['$p_{\\rm fa}=10^{%d}$, $p_{{\\rm md},%d}=',...
    '1-\\Pr\\{{\\cal E}_{%d,{\\rm approx}}\\}$ (an.) ',...
    'vs. $1-\\Pr\\{{\\cal E}_{%d}\\}$ (sim.) [nEst, mAv, mMed] with $L=%d$'],...
    round(log10(p_fa)), p_md_id, p_md_id, p_md_id, L_rng(iL) );
end
title(tit_str, ltx{:}, 'FontSize', fntSz_title);

fName_out = fullfile('..', 'text', 'radio_eng_article', 'figs',...
    sprintf('p_md_%d_cases_nEst_mAv_mMed.pdf', p_md_id) );
if p_d_midVal_vs_true == 2
fName_out = fullfile('..', 'text', 'radio_eng_article', 'figs',...
    sprintf('p_md_%d_cases_truePmd_nEst_mAv_mMed.pdf', p_md_id) );
end
% export_fig(fName_out,'-pdf','-transparent');

%% plots of p_d for bigger frequency range --- p_d vs F_eps
lnWid = {'LineWidth', 2.2};
fName_in = '2025-06-30_12-33-11_p_d_vs_snr_freqOff_hGeneral_N_framesN_rep12000.mat';
load(fullfile('..', outDir, fName_in));
fntSz_leg = 11;
fntSz_axLabel = 14;
fntSz_title = 14;

p_d_midVal_vs_true = 1; %1 or 2
p_fa_id = 2; %
% legCell = arrayfun(@(x) sprintf('$f_{\\varepsilon}$=%.3f', x), freq_rng,...
%     'UniformOutput', false);
p_md_id = 2;

p_md_vs_pd = 1; %or 2
p_md_vs_pd_str = 'md';
if p_md_vs_pd == 2
    p_md_vs_pd_str = 'd';
end

snr_dB_idx = [1];
N_snr_idc = numel(snr_dB_idx);
legCell = arrayfun(@(x) sprintf('$N_{\\rm seq} = %d$', x), N_seq_rng,...
    'UniformOutput', false);

xlbl = '$F_{\varepsilon}$ [-]';
ylbl = '$p_{\rm d}$ [-]';

p_fa = p_fa_rng(p_fa_id);
p_fa_id_ = (p_fa_id - 1)*2 + p_d_midVal_vs_true;
N_freq_an = 300;
% freq_rng_an = linspace(min(freq_rng), max(freq_rng)+0.005, N_freq_an).';
freq_rng_an = linspace(min(freq_rng), 0.06, N_freq_an).';

% xlims = [min(freq_rng), max(freq_rng)];
xlims = [min(freq_rng_an), max(freq_rng_an)];
ylims = [(1.5e-5), 1];
for iCase = 7%1:N_cases
    h = h_preambs{iCase};

    for iN_seq = 1:N_seq_rng_numel
%         iColor = iN_seq+2;
        iColor = iN_seq;
        if p_md_id == 1
            % *** sim y'_1 ***
            p_md = zeros(1, N_freq);
            for iFreq = 1:N_freq
                p_md_aux = 1 - p_d_1_ideal(snr_dB_idx, p_fa_id_, iN_seq, iCase, iFreq);
                p_md(:,iFreq) = p_md_aux;
            end
            if p_md_vs_pd == 1
                semilogy(freq_rng, p_md(1,:), 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle', '--', 'Marker','x');
            else
                semilogy(freq_rng, 1-p_md(1,:), 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle', '--', 'Marker','x');
            end
        else
        % *** sim y'_2 ***
            p_md = zeros(1, N_freq);
            for iFreq = 1:N_freq
                p_md_aux = 1 - p_d_2_ideal(snr_dB_idx, p_fa_id_, iN_seq, iCase, iFreq);
                p_md(:,iFreq) = p_md_aux;
            end
            if p_md_vs_pd == 1
                semilogy(freq_rng, p_md(1,:), 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle', '--', 'Marker','x');
            else
                semilogy(freq_rng, 1-p_md(1,:), 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle', '--', 'Marker','x');
            end
        end
        if iN_seq == 1
            hold on
        end
    end

    % write auxiliary lines -- the lobe border
    N_h = numel(h);
    for iN_seq = 1:N_seq_rng_numel
%         iColor = iN_seq +2;
        iColor = iN_seq;
        x_coord = repmat( N_seq_rng(iN_seq)/N_h, 2,1);
        y_coord = [1e-6, 1];
        semilogy(x_coord, y_coord, 'Color', colors(iColor,:),...
                'LineWidth',1.5, 'LineStyle','-');
    end

    for iN_seq = 1:N_seq_rng_numel
%         iColor = iN_seq +2;
        iColor = iN_seq;
%         subplot(2,2,iN_seq);
        if p_md_id == 1
            % *** analytical y'_1 ***
            p_md = zeros(1, N_freq);
            for iFreq = 1:N_freq_an
            [p_md(:,iFreq), ~] = get_p_md(snr_dB_rng(snr_dB_idx), h, ...
                N_seq_rng(iN_seq), p_fa, freq_rng_an(iFreq));
            end
            if p_md_vs_pd == 1
                semilogy(freq_rng_an, p_md(1,:), 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle',':');
            else
                semilogy(freq_rng_an, 1-p_md(1,:), 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle',':');
            end
        else
            % *** analytical y'_2 ***
            p_md = zeros(1, N_freq);
            for iFreq = 1:N_freq_an
            [~, p_md(:,iFreq)] = get_p_md(snr_dB_rng(snr_dB_idx), h, ...
                N_seq_rng(iN_seq), p_fa, freq_rng_an(iFreq));
            end
            if p_md_vs_pd == 1
                semilogy(freq_rng_an, p_md(1,:), 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle',':');
            else
                semilogy(freq_rng_an, 1-p_md(1,:), 'Color', colors(iColor,:),...
                    lnWid{:}, 'LineStyle',':');
            end
        end
%         if iN_seq == 1
%             hold on
%         end
    end
    hold off

    grid on
    grid minor

%         xlim([xlims(1), snr_dB_rng(max_nonZero_val_id)]);
    xlim(xlims);
    ylim(ylims);
    xlabel(xlbl, ltx{:}, 'FontSize', fntSz_axLabel);
    if p_md_vs_pd == 1
        ylabel('$p_{\rm md}$ [-]', ltx{:}, 'FontSize', fntSz_axLabel);
    else
        ylabel('$p_{\rm d}$ [-]', ltx{:}, 'FontSize', fntSz_axLabel);
    end
%             if idxNrep == 1 && iFreq == 11
    legend(legCell{:},ltx{:}, 'Location','northeast', 'FontSize', fntSz_leg);
%             end
    
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',fntSz_axLabel-2);
    tit_str = sprintf(['$p_{\\rm fa}=10^{%d}$, ${\\rm SNR}=%d$ dB,',...
        ' id:%02d, $p_{{\\rm %s}, %d}$'], ...
        round(log10(p_fa)), snr_dB_rng(snr_dB_idx),...
        iCase, p_md_vs_pd_str, p_md_id);
    title(tit_str, ltx{:}, 'FontSize', fntSz_title);

%     pause
end

fName_out = fullfile('..', 'text', 'radio_eng_article', 'figs',...
    sprintf('p_%s_%d_case%d_ideal.pdf', p_md_vs_pd_str, p_md_id, iCase) );
if p_d_midVal_vs_true == 2
fName_out = fullfile('..', 'text', 'radio_eng_article', 'figs',...
    sprintf('p_%s_%d_case%d_ideal.pdf', p_md_vs_pd_str, p_md_id, iCase) );
end
% export_fig(fName_out,'-pdf','-transparent');