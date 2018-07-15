function check_WC_params(par)

if isfield(par,'features') && strcmp('features','wav')

    if isfield(par,'scales')&&  isfield(par,'par.w_post') &&  isfield(par,'par.w_pre')
         L = wmaxlev(par.w_pre+par.w_post,'haar');
         if L<scale
             error('[par.scale] exceeds maximum wavelet decomposition level for a waveform of length [par.w_pre+par.w_post]')
         end
    end
    if isfield(par,'par.w_post') &&  isfield(par,'par.w_pre')
        x = log2(par.w_pre+par.w_post);
        if floor(x)~=x
            error('The length [par.w_pre+par.w_post] should be a power of 2')
        end
    end
end
end