function [path, num_rxn, err] = grpath(s1, s2, coeff_in, coeff_o, inv_sets)

[num_rxns, num_species] = size(coeff_in);
rxn_coeffs = coeff_o - coeff_in;

s = 1;
c = ones(1,num_rxns);
b = s2 - s1; b = b';
a = rxn_coeffs';

lb=zeros(num_rxns,1); ub=[];
ctype = repmat('S', num_species, 1);
vartype = repmat('I', num_rxns, 1);
param.msglev=1;
[xmin,fmin,status,extra]=glpk(c,a,b,lb,ub,ctype,vartype,s,param);

if status ==2 | status == 5
  num_rxn = xmin';
%   num_set = ceil(rand(1,size(inv_sets,1)));
%   num_rxn = num_rxn + num_set * inv_sets;
  tot_num = sum(num_rxn);
  
%   i_rxns = zeros(1,tot_num);
  
  i_rxns = []; 
  for i=1: num_rxns
      i_rxns = [i_rxns i*ones(1, num_rxn(i))];
  end
  perm_idx = randperm(tot_num);
  rxns = i_rxns(perm_idx);
  
%   im_state = zeros(tot_num+1, num_species);
%   inter_state(1,:) = s1;
  im_state = s1;
  
  err = 0;
  for k=1:tot_num    
       if any(im_state - coeff_in(rxns(k),:) <0) %|| inter_state(k+1,2)<0 || inter_state(k+1,3)<0
          err = 1;
          break
       end
       im_state  = im_state + rxn_coeffs(rxns(k),:);
  end
  t = rand(1,tot_num);
%   taus = dt .* diff(sort(t_rxn));
  path.trxns = sort(t);
  path.rxns = rxns;
else
    error('can''t initialze path')
end