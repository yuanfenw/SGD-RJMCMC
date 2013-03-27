function [grdtsum, state, rxns, trxns] = gradient(state, t, input_coeffs, output_coeffs, rxn_rates, inv_sets, i_max, rxns, trxns)
% gradient calculation func
s = length(t);

rxn_coeffs = output_coeffs - input_coeffs;   
[num_rxns, num_species] = size(rxn_coeffs);

burnin = 100;
num_sets = size(inv_sets, 1);
maxnum_set = 20*ones(1, num_sets);
inter_state = zeros(100, num_species);
state_his = zeros((i_max+burnin)/10, num_species);
% num_his = zeros(1,(i_max+burnin)/10); 
% i_max =  500; % MCMC sample size for each interval    
grdtsum = zeros(1, num_rxns); %should be num_paras in general
k2 = 1:num_rxns;

% function F = getF(s)
% %    combinations = @(n,k) prod((ones(1,k)+(n-1)) - [0:(k-1)]);
%     F = ones(1, num_rxns);
%     for r=1:num_rxns
%         for p = 1: length(s)
%             F(r) = F(r) * prod((ones(1,input_coeffs(r,p))+(s(p)-1)) - [0:(input_coeffs(r,p)-1)]);
% %                 F(r) = F(r) * nchoosek(s, input_coeffs(r,p)) * factorial(input_coeffs(r,p));
%         end
%     end
% end

% function a = llh(sstart, send, tstart, tend)
%     if isempty(ts2)
%         a = -sum(rxn_rates .* getF(sstart))*(tend-tstart);
%     else
%        taus2 = diff([tstart, ts2, tend]);
%        a = 0;
%        im_state = sstart;
%         for k=1:length(rxns2)
%            haz = rxn_rates.* getF(im_state);
%            if haz(rxns2(k)) == 0
%                error('log of zero');
%            end
%            a = a - sum(haz)*taus2(k) + log(haz(rxns2(k)));
%            im_state = im_state +rxn_coeffs(rxns2(k),:);
%         end
%         if any(im_state ~= send)
%             error ('state doesnt match')
%         end
%         a = a - sum(rxn_rates .* getF(im_state))* taus2(end);
%     end
% end


% prev_lh = zeros(1,s-2);
% prev_numset = 1;

grdt = zeros(1, num_rxns);
jg = 0;

inv_setp = [1     0     0     0     0     0     0     0;
     0     1     0     0     0     0     0     0;
     0     0     1     1     0     0     0     0;
     0     0     0     0     1     1     0     0;
     0     0     0     0     0     0     1     1;
     0     0     0     0     1     0     2     0;
     0     0     0     0     0     1     0     2];
i = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%  start sampling  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

while i <= (i_max+burnin)
    
%update the first state
    [nons idx2] = max(trxns>t(2));
    idx2 = idx2 -1;
    if idx2>=1
        rxns1 = rxns(1:idx2);
        ts1 = trxns(1:idx2);
    else
        rxns1 = [];
        ts1=[];
    end
    
    if isempty(ts1)
        prev_lh = -sum(rxn_rates .* getF(state(1,:), input_coeffs))*(t(2)-t(1));
    else
       taus1 = diff([t(1), ts1, t(2)]);
       prev_lh = 0;
       im_state = state(1,:);
        for k=1:length(rxns1)
           haz = rxn_rates.* getF(im_state, input_coeffs);
           if haz(rxns1(k)) == 0
               error('log of zero');
           end
           prev_lh = prev_lh - sum(haz)*taus1(k) + log(haz(rxns1(k)));
           im_state = im_state +rxn_coeffs(rxns1(k),:);
        end
%             if any(im_state ~= state(j+2,:))
%                 error ('state doesnt match')
%             end
        prev_lh = prev_lh - sum(rxn_rates .* getF(im_state, input_coeffs))* taus1(end);
    end
    
    set_idx = floor(rand*size(inv_setp, 1))+1;
    inv_set = inv_setp(set_idx,:);
    nzeroid = find(inv_set);

    updaterand = rand;
    if updaterand < 0.25 %& cond_del
        mtype = 1;
    elseif updaterand < 0.5
        mtype = 2;
    else
        mtype = 3;
    end

    if ~isempty(rxns1)
        for id = 1:num_rxns
        num_rxn(id) = sum(rxns1==id);
        end
    else
        num_rxn = zeros(1,8);
    end
    [newpath, new_numrxn] = grpath3(rxns1, ts1, t(1), t(2), mtype, inv_set, num_rxn); 

    if  min(new_numrxn)>=0
        tot_num = sum(new_numrxn);

        err = 0;
        inter_state(1,:) = state(1,:);
        if mtype == 1
            for k=1:num_rxns
              inter_state(1,:) = inter_state(1,:) + inv_set(k)*rxn_coeffs(k,:);
            end
            
        elseif mtype == 2
            for k=1:num_rxns
              inter_state(1,:) = inter_state(1,:) - inv_set(k)*rxn_coeffs(k,:);
            end
        end
%                 new_midstate=im_state;
        if min(inter_state(1,:))<0
            err =1;
        else
        for k=1:tot_num
            if k> length(newpath.rxns) || newpath.rxns(k)>8
                error('exceed dimension');
            end
            inter_state(k+1,:) = inter_state(k,:) + rxn_coeffs(newpath.rxns(k),:);
            if min(inter_state(k,:)- input_coeffs(newpath.rxns(k),:))<0
                err = 2;
                break
            end
        end    
        end
        
    end
    if min(new_numrxn)<0 || err ~= 0
%             a(j) = a(j) +1;
        acpt_ratio = 0;

    else
        rxns1 = newpath.rxns;
        ts1 = newpath.trxns;
       % calculte prob. density of the path (with two adjacent data points)
       if isempty(ts1)
            lh = -sum(rxn_rates .* getF(inter_state(1,:), input_coeffs))*(t(2)-t(1));
       else
            taus1 = diff([t(1), ts1, t(2)]);
            lh = 0;
%            inter_state(1,:) = state(j,:);
            for k=1:length(ts1)
               haz = rxn_rates.* getF(inter_state(k,:), input_coeffs);
               if haz(rxns1(k)) == 0
                   error('log of zero');
               end
               lh = lh - sum(haz)*taus1(k) + log(haz(rxns1(k)));
%                inter_state(k+1,:)= inter_state(k,:) +rxn_coeffs(rxns1(k),:);
            end
            lh = lh - sum(rxn_rates .* getF(inter_state(k+1,:), input_coeffs))* taus1(tot_num+1);
            if inter_state(k+1,:) ~= state(2,:)
                error('err: state unmatch')
            end
       end

       % calculate the acceptance ratio according to move types
       inv_setnz = inv_set(nzeroid);
       nr = new_numrxn(nzeroid);
       if mtype == 2
           acpt_ratio = exp(lh-prev_lh) * (t(2)-t(1))^sum(inv_setnz)/combn(nr, inv_setnz);
       elseif mtype == 1
           acpt_ratio = exp(lh-prev_lh) * combn(nr+inv_setnz, inv_setnz)/(t(2)-t(1))^sum(inv_setnz);
       else
%            error('mtype error')
           acpt_ratio = exp(lh-prev_lh);
       end

    end
    % accept or reject new sample with prob acpt_ratio
    if acpt_ratio > 1 || acpt_ratio > rand
%        prev_lh = lh;
       state(1,:) = inter_state(1,:);
       rxns = [rxns1 rxns(idx2+1:end)];
       trxns = [ts1 trxns(idx2+1:end)];
    end
        
    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update reaction within paired intervals
    for j=1:s-2
        
        [nons idx1] = max(trxns>t(j));
        if j == s-2
            idx2 = length(trxns);
        else
        [nons2 idx2] = max(trxns > t(j+2));
        idx2 = idx2 -1;
        end
        
        if idx2>=idx1
            rxns2 = rxns(idx1:idx2);
            ts2 = trxns(idx1:idx2);
        else
            rxns2 = [];
            ts2=[];
        end
%         prev_lh = llh(state(j,:), state(j+2,:), t(j), t(j+2));
        if isempty(ts2)
            prev_lh = -sum(rxn_rates .* getF(state(j,:), input_coeffs))*(t(j+2)-t(j));
        else
           taus2 = diff([t(j), ts2, t(j+2)]);
           prev_lh = 0;
           im_state = state(j,:);
            for k=1:length(rxns2)
               haz = rxn_rates.* getF(im_state, input_coeffs);
               prev_lh = prev_lh - sum(haz)*taus2(k) + log(haz(rxns2(k)));
               im_state = im_state +rxn_coeffs(rxns2(k),:);
            end
%             if any(im_state ~= state(j+2,:))
%                 error ('state doesnt match')
%             end
            prev_lh = prev_lh - sum(rxn_rates .* getF(im_state, input_coeffs))* taus2(end);
        end
  
        
%       for j2 = 1:5
%           %recalculate index of 
%         [nons idx1] = max(trxns>t(j));
%         if idx1~=idx_start
%             error('index 1')
%         end
%         
%         if j == s-2
%             idx2 = length(trxns);
%         else
%         [nons2 idx2] = max(trxns > t(j+2));
%         idx2 = idx2 -1;
%         end
%         if idx2<idx1
%             error('index 2 less than index 1')
%         end
        
        set_idx = floor(rand*num_sets)+1;
        inv_set = inv_sets(set_idx,:);
        nzeroid = find(inv_set);
               
    
%     if j<=s-2 % update two interval at a time
            
        % update path by either deleting or inserting null set
        
        updaterand = rand;
        if updaterand < 0.25 %& cond_del
            mtype = 1;
        elseif updaterand < 0.5
            mtype = 2;
        else
            mtype = 3;
        end

        
% generate new path between t(j) and t(j+2)
        
        if ~isempty(rxns2)
            for id = 1:num_rxns
            num_rxn(id) = sum(rxns2==id);
            end
        else
            num_rxn = zeros(1,8);
        end
            
        [newpath, new_numrxn] = grpath3(rxns2, ts2,  t(j), t(j+2), mtype, inv_set, num_rxn);                
        
        if  min(new_numrxn)>=0
            tot_num = sum(new_numrxn);
            if tot_num == 0  % if there is no rxns between t_j and t_j+2
                new_midstate = state(j,:);
                err = 0;
            else
                [nons idx3] = max(newpath.trxns>t(j+1)); % find the rxn index before t(j+1)
                if nons == 0
                    idx3 = tot_num;
                elseif idx3 == 1
                    idx3 = 0;
                else
                    idx3 = idx3-1;
                end

                err = 0;
                if idx3 ==0
                    if any(state(j,3:4) ~= state(j+1,3:4))
                        err = 1;
                    else
                        new_midstate = state(j,:);
                    end
                end
                if err ~= 1 
                    inter_state(1,:) = state(j,:);
    %                 new_midstate=im_state;
                    for k=1:tot_num
                        if k> length(newpath.rxns) || newpath.rxns(k)>8
                            error('exceed dimension');
                        end
                        inter_state(k+1,:) = inter_state(k,:) + rxn_coeffs(newpath.rxns(k),:);
                        if min(inter_state(k,:)- input_coeffs(newpath.rxns(k),:))<0
                            err = 2;
                            break
                        end           

                        %check the middle point
                        if k == idx3
                            if any(inter_state(k+1,3:4) ~= state(j+1,3:4))
                                err = 3;
                                break
                            else
                                new_midstate=inter_state(k+1,:);
                            end
                        end
                    end
                end
            end            
        end
        
        
        if min(new_numrxn)<0 || err ~= 0
%             a(j) = a(j) +1;
            acpt_ratio = 0;

        else
            rxns2 = newpath.rxns;
            ts2 = newpath.trxns;
            if length(ts2)~=tot_num
                error('err3')
            end
            
            
           % calculte prob. density of the path (with two adjacent data points)
           if isempty(ts2)
                lh = -sum(rxn_rates .* getF(state(j,:), input_coeffs))*(t(j+2)-t(j));
           else
                taus2 = diff([t(j), ts2, t(j+2)]);
                lh = 0;
    %            inter_state(1,:) = state(j,:);
               for k=1:length(ts2)
                   haz = rxn_rates.* getF(inter_state(k,:), input_coeffs);
                   if haz(rxns2(k)) == 0
                       error('log of zero');
                   end
                   lh = lh - sum(haz)*taus2(k) + log(haz(rxns2(k)));
    %                inter_state(k+1,:)= inter_state(k,:) +rxn_coeffs(rxns2(k),:);
               end
               lh = lh - sum(rxn_rates .* getF(inter_state(tot_num+1,:), input_coeffs))* taus2(tot_num+1);
           end
           % calculate the acceptance ratio according to move types
           inv_setnz = inv_set(nzeroid);
           nr = new_numrxn(nzeroid);
           if mtype == 2
               acpt_ratio = exp(lh-prev_lh) * (t(j+2)-t(j))^sum(inv_setnz)/combn(nr, inv_setnz);
           elseif mtype == 1
               acpt_ratio = exp(lh-prev_lh) * combn(nr+inv_setnz, inv_setnz)/(t(j+2)-t(j))^sum(inv_setnz);
           else
    %            error('mtype error')
               acpt_ratio = exp(lh-prev_lh);
           end

        end
    % accept or reject new sample with prob acpt_ratio
        if acpt_ratio > 1 || acpt_ratio > rand
%            prev_lh = lh;
           state(j+1,:) = new_midstate;
           rxns = [rxns(1:idx1-1) rxns2 rxns(idx2+1:end)];
           trxns = [trxns(1:idx1-1) ts2 trxns(idx2+1:end)];
        end
%       end 
    end
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    %update the last state of unobserved
    [nons idx1] = max(trxns>t(end-1));
    if nons == 0
        rxns1 = [];
        ts1=[];
    else
        idx2 = length(trxns);
        rxns1 = rxns(idx1:idx2);
        ts1 = trxns(idx1:idx2);
    end
    
    if isempty(ts1)
        prev_lh = -sum(rxn_rates .* getF(state(s-1,:), input_coeffs))*(t(s)-t(s-1));
    else
       taus1 = diff([t(s-1), ts1, t(s)]);
       prev_lh = 0;
       im_state = state(s-1,:);
        for k=1:length(rxns1)
           haz = rxn_rates.* getF(im_state, input_coeffs);
           if haz(rxns1(k)) == 0
               error('log of zero');
           end
           prev_lh = prev_lh - sum(haz)*taus1(k) + log(haz(rxns1(k)));
           im_state = im_state +rxn_coeffs(rxns1(k),:);
        end
%             if any(im_state ~= state(j+2,:))
%                 error ('state doesnt match')
%             end
        prev_lh = prev_lh - sum(rxn_rates .* getF(im_state, input_coeffs))* taus1(end);
    end
    
    set_idx = floor(rand*size(inv_setp, 1))+1;
    inv_set = inv_setp(set_idx,:);
    nzeroid = find(inv_set);

    updaterand = rand;
    if updaterand < 0.25 %& cond_del
        mtype = 1;
    elseif updaterand < 0.5
        mtype = 2;
    else
        mtype = 3;
    end

    if ~isempty(ts1)
        for id = 1:num_rxns
        num_rxn(id) = sum(rxns1==id);
        end
    else
        num_rxn = zeros(1,8);
    end
    [newpath, new_numrxn] = grpath3(rxns1, ts1, t(s-1), t(s), mtype, inv_set, num_rxn); 

    if  min(new_numrxn)>=0
        tot_num = sum(new_numrxn);

        err = 0;
        inter_state(1,:) = state(s-1,:);

        for k=1:tot_num
%             if k> length(newpath.rxns) || newpath.rxns(k)>8
%                 error('exceed dimension');
%             end
            inter_state(k+1,:) = inter_state(k,:) + rxn_coeffs(newpath.rxns(k),:);
            if min(inter_state(k,:)- input_coeffs(newpath.rxns(k),:))<0
                err = 2;
                break
            end
        end    
     
    end
    if min(new_numrxn)<0 || err ~= 0
%             a(j) = a(j) +1;
        acpt_ratio = 0;

    else
        rxns1 = newpath.rxns;
        ts1 = newpath.trxns;
       % calculte prob. density of the path (with two adjacent data points)
       if isempty(ts1)
            lh = -sum(rxn_rates .* getF(inter_state(1,:), input_coeffs))*(t(s)-t(s-1));
       else
            taus1 = diff([t(s-1), ts1, t(s)]);
            lh = 0;
%            inter_state(1,:) = state(j,:);
            for k=1:length(ts1)
               haz = rxn_rates.* getF(inter_state(k,:), input_coeffs);
               if haz(rxns1(k)) == 0
                   error('log of zero');
               end
               lh = lh - sum(haz)*taus1(k) + log(haz(rxns1(k)));
%                inter_state(k+1,:)= inter_state(k,:) +rxn_coeffs(rxns1(k),:);
            end
            lh = lh - sum(rxn_rates .* getF(inter_state(k+1,:), input_coeffs))* taus1(tot_num+1);
       end

       % calculate the acceptance ratio according to move types
       inv_setnz = inv_set(nzeroid);
       nr = new_numrxn(nzeroid);
       if mtype == 2
           acpt_ratio = exp(lh-prev_lh) * (t(s)-t(s-1))^sum(inv_setnz)/combn(nr, inv_setnz);
       elseif mtype == 1
           acpt_ratio = exp(lh-prev_lh) * combn(nr+inv_setnz, inv_setnz)/(t(s)-t(s-1))^sum(inv_setnz);
       else
%            error('mtype error')
           acpt_ratio = exp(lh-prev_lh);
       end

    end
    % accept or reject new sample with prob acpt_ratio
    if acpt_ratio > 1 || acpt_ratio > rand
%        prev_lh = lh;
        if ~isempty(ts1)
             state(s,:) = inter_state(k+1,:);
        else
            state(s,:) = state(s-1,:);
        end
       rxns = [rxns(1:idx1-1) rxns1];
       trxns = [trxns(1:idx1-1) ts1];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    
    if i>burnin %&& mod(i,10) == 0
        
        taus = diff([t(1), trxns, t(end)]);
        im_state = state(1,:); F = getF(state(1,:), input_coeffs);
        grdt = - taus(1) .* F; 
        for k3 = 1:length(trxns)
            im_state(:) = im_state + rxn_coeffs(rxns(k3),:);
            F(:) = getF(im_state, input_coeffs);
            grdt(:) = grdt + (k2 == rxns(k3))./rxn_rates - taus(k3+1) .* F;
        end
        grdtsum = grdtsum + grdt;
        
    end
    
    i=i+1;
    
end

grdtsum = grdtsum/(i_max);
end