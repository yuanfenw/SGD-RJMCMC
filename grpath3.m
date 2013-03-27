function [newpath, num_rxn] = grpath3(rxns, t, t1, t2, mtype, inv_set, num_rxn)
% generate new path by different type of move
%a = isstruct(spath)

tot_num = length(t);

rxnid = find(inv_set);
inv_setnz = inv_set(rxnid);
if mtype == 1
    for k = 1:length(rxnid)
        pos_rxnid = find(rxns==rxnid(k));  % index of corresponding reactions
        if length(pos_rxnid) ~= num_rxn(rxnid(k))
           error('err4');
        end
        if num_rxn(rxnid(k)) < inv_setnz(k)             
            break;
        else
    % randomly delete /inv_set(rxnid(k)) number of reaction k 
        idx2 = randperm(length(pos_rxnid));
        pos_rxnid = pos_rxnid(idx2);
        t(pos_rxnid(1:inv_setnz(k))) = [];
        rxns(pos_rxnid(1:inv_setnz(k))) = [];
        end
    end
    num_rxn = num_rxn - inv_set;
    
elseif mtype == 2

    for k =1:length(rxnid)
        pt = t1+(t2-t1)*rand(1,inv_set(rxnid(k)));
        pr = rxnid(k)*ones(1, inv_set(rxnid(k)));
        
        t = [t pt]; rxns = [rxns pr];
        [t id] = sort(t);
        rxns = rxns(id);
%         if any(t>pt)
%             [a, idx] = max(pt<t);
%         else
%             idx=tot_num+1;
%         end
%     
%         new_t = [t(1:idx-1), pt, t(idx:tot_num)];
%         new_rxns = [rxns(1:idx-1), pr, rxns(idx:tot_num)];
    end
    num_rxn = num_rxn + inv_set;

else
    
    if tot_num > 1
        idx = randperm(tot_num);
        rxns = rxns(idx);
    end

end

newpath.trxns = t;
newpath.rxns = rxns;

