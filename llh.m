function a = llh(rxns2, ts2, sstart, send, tend, tstart, rxn_rates)
    if isempty(ts2)
        prev_lh = -sum(rxn_rates .* get_F(sstart))*(tend-tstart);
    else
       taus2 = diff([tstart, ts2, tend]);
       prev_lh = 0;
       im_state = sstart;
        for k=1:length(rxns2)
           haz = rxn_rates.* get_F(im_state);
           if haz(rxns2(k)) == 0
               error('log of zero');
           end
           prev_lh = prev_lh - sum(haz)*taus2(k) + log(haz(rxns2(k)));
           im_state = im_state +rxn_coeffs(rxns2(k),:);
        end
        if any(im_state ~= send)
            error ('state doesnt match')
        end
        prev_lh = prev_lh - sum(rxn_rates .* get_F(im_state))* taus2(end);
    end
end