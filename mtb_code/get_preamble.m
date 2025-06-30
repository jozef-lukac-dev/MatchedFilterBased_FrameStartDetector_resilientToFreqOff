function preamble_h = get_preamble(str_id, N_h0, N_rep, u)
% function GET_PREAMBLE
% preamble_h = get_preamble(str_id, N_h0, N_rep, u)
% Creates samples of a preamble sequence.
% Inputs:
%  <str_id> -- a string 'wifi' (for wifi preamble), 'zc' for Zadoff-Chu
%  sequence or 'pn' for pseudorandom complex gaussian sequence.
%  <N_h0> -- 1x1 number. Length of the one period of the preamble
%  <N_rep> -- 1x1 number. Number of repetitions of base sequence
%   (default = 1).
%  <u> -- 1x1 number. A parameter for ZC-seunce (default = 1).
% Outputs:
%  <preamble_h> -- (N_rep * N_h0) x 1 vector of preamble samples. It is
%  noramlized to unit energy.
%

    if nargin < 3
        N_rep = 1;
    end
    if nargin < 4
        u = 1;
    end
    switch upper(str_id)
        case 'WIFI'
            assert(N_h0 == 16,['N_h0 must be 16 for wifi preamble ',...
                'base sequence (N_h0 == %d).'],N_h0)
%             CP_PREAMB_LEN = 32;
            FFT_SIZE = 64;
            mod64 = @(x) mod(x,64);
            
            %short sequence
            s = sqrt(13/6);
            auxSymb = s+1i*s;
            symb = [auxSymb, -auxSymb];
            assign_id = [...
            -24, -20, -16, -12, ...
	         -8,  -4,   4,   8, ...
	         12,  16,  20,  24];
            assign_val = 1+ [...
            0, 1, 0, 1, ...
	        1, 0, 1, 1, ...
	        0, 0, 0, 0];
        
            one_period_short = zeros(FFT_SIZE,1);
            one_period_short(mod64(assign_id) +1) = symb(assign_val);
        
            one_period_short_ifft = ifft(one_period_short);
            h0 = one_period_short_ifft(1:16);
%             samps_shortSeq = [one_period_short_ifft(end-CP_PREAMB_LEN+1:end);...
%                 one_period_short_ifft; one_period_short_ifft; ...
%                  one_period_short_ifft(1)];%*0.5];
%             samps_shortSeq(1) = 0.5 * samps_shortSeq(1);

        case 'ZC'
            n = (0:(N_h0-1)).';
%             u = 1;
            cf = mod(N_h0,2);
            h0 = exp(-1i*pi*u*n.*(n + cf)/N_h0);
        case 'PN'
            h0 = 1/sqrt(2) * (randn(N_h0,1) + 1i*randn(N_h0,1));
        otherwise
            error('Error: Not known option: ''%s''\n', str_id);
    end

    preamble_h = repmat(h0, N_rep, 1);
    E_h0 = sum( abs(preamble_h).^2,1 ); %normalize energy
    preamble_h = preamble_h / sqrt(E_h0);
end