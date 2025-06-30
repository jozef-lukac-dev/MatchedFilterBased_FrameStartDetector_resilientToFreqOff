
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

h = h_preambs{4};
N_seq = 5;
p_fa_req = 1e-4;
abs_absSq = 'abs';
mAv_mMed = 'mAv';
N_samps = 1e7;
L_av = 30;

[r_thr, p_fa_final, iIter] = set_thr_forGiven_pFa(...
   h, N_seq, p_fa_req, abs_absSq, mAv_mMed, N_samps, L_av);

fprintf('thr: [%.3f], p_fa: [%0.5f], iIter: [%02d]\n',...
    r_thr, p_fa_final, iIter);
