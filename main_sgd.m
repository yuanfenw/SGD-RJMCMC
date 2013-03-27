% Parameter Estimation by Stochastic Gradient Descent
% Yuanfeng Wang (c) 10/4/09


%rd = matlabroot;
cd(['/home/wjb/src/SGD'])
addpath(['/home/wjb/src/SGD/metatool'])
addpath(['/home/wjb/src/SGD/glpkmex-2.4'])
%specify the stochiometry matrix

% usegraph=0;
% reply = input('Do you want to use graphic interface? Y/N [Y]: ', 's');
% if isempty(reply) || strcmp(reply,'Y') || strcmp(reply,'y')
%    usegraph=1;
% else
%    usegraph=0;
% end

% if usegraph
%    [modelfile modelpath] = uigetfile('.xml','Please specify the SBML model');
%    modelpath = [modelpath modelfile];
% else
%    modelpath = input('SBML model path: ', 's');
% end

usegraph=0;
reply = input('Do you want to import stoichiometry matrix from files? Y/N [Y]: ', 's');
if isempty(reply) || strcmp(reply,'Y') || strcmp(reply,'y')
   usegraph=1;
else
   usegraph=0;
end

if usegraph
   [Smat_file Smat_path] = uigetfile('.dat','Please specify the file for Stoichiometry matrix of reactants');
   Smat_rct = [Smat_path Smat_file];
   
   [Smat_file Smat_path] = uigetfile('.dat','Please specify the file for Stoichiometry matrix of products');
   Smat_prdt = [Smat_path Smat_file];
   
   input_coeffs = importdata(Smat_rct);
   output_coeffs = importdata(Smat_prdt);
else
    disp('Please Specify the stoichiometry matrix in "main_sgd.m"')
    % Specify the stoichiometry matrix here
    input_coeffs =  [1 0 0 0 1; 0 1 0 0 0; 1 0 0 0 0; 0 0 1 0 0; 0 0 0 2 0; 0 0 0 0 1; 0 0 1 0 0; 0 0 0 1 0];
    output_coeffs = [0 1 0 0 0; 1 0 0 0 1; 1 0 1 0 0; 0 0 0 0 0; 0 0 0 0 1; 0 0 0 2 0; 0 0 1 1 0; 0 0 0 0 0];
end

% Specify file of measurement Data
disp('Please Specify the file of observation data:')
if 1
    [datafile datapath] = uigetfile('.dat', 'Please specify the data file');
    datapath = [datapath datafile];
else
    datapath = input('Data file:', 's');
end
A = importdata(datapath);
t = A.data(:,1);
data= A.data(:,2:end);

% set initial step size
reply2 = input('Please input the step size of gradient descent:', 's');
if isempty(reply2)
    speed = 0.01;
    disp('Default step size: 0.01')
else
    speed = str2num(reply2);
end

% Model = TranslateSBML(modelpath);
k_max = 5;
ssize = 100;  % sample size of each RJMCMC sampling step

[num_rxns, num_species] = size(input_coeffs);
num_para = num_rxns;



%%% calculate null-set %%%
rxn_coeffs = output_coeffs - input_coeffs;   
net.st = rxn_coeffs';
net.irrev_react= ones(1,num_rxns);
net = metatool(net);


rxns_ini = [];
trxns_ini = [];
num_rxn = zeros(1,8);
% num_rxn =  zeros(s-1,8);

for j=1:length(t)-1
%     dt = t(j+1)-t(j);
    err = 1;
    while err
    [path, num_rxn, err] = grpath(data(j,:), data(j+1,:), input_coeffs, output_coeffs, net.sub);
    end
    rxns_ini = [rxns_ini  path.rxns];
    trxns_ini = [trxns_ini t(j)+(t(j+1)-t(j))*path.trxns];
end

rxn_rates = ones(1,8);
rrate_exp = log(rxn_rates);
state_ini = data; 
min_grad = 10000*ones(1, num_para);
[grdt, state, rxns, trxns] = gradient(state_ini, t, input_coeffs, output_coeffs, rxn_rates, net.sub, 10, rxns_ini, trxns_ini);

diary('diary.txt')
k=1;
rxn_rates_k(1,:) = rxn_rates; 
gradient_k(1,:) = grdt;
disp('Initial values of rate paramters:')
disp(num2str(rxn_rates))

while k<k_max
    k = k+1;
    ssize = 500*ceil(k/10);
%     adpatively change the step size

    cont_decrease_speed = 1;
    while cont_decrease_speed
        d_rate_exp = speed .* grdt .* rxn_rates;
        rrate_exp = rrate_exp + d_rate_exp;
        rxn_rates = exp(rrate_exp); 
        [grdt, state_new, rxns_new, trxns_new] = gradient(state, t, input_coeffs, output_coeffs, rxn_rates, net.sub, ssize, rxns, trxns);
        
        cont_decrease_speed = 0;
        if sum(abs(grdt) > abs(gradient_k(k-1,:)))>=4 & speed>=0.005 & max(abs(d_rate_exp)) > 5.e-03
            speed = speed /2;
            cont_decrease_speed = 0;
            
        elseif  abs(grdt) < abs(gradient_k(k-1,:))
            if rand<0.5
                speed = speed*2;
%             else
%                 speed = speed/2;
            end
            
        else
            if rand<1/2 & speed>=0.005
                speed = speed /2;
            else   speed = speed*2; 
            end
        end
    
%         grdt = gradient_new(data, t, input_coeffs, output_coeffs, rxn_rates, net.sub);
    end 
    disp(['current step: ' num2str(k-1)])
    disp(['step size is ' num2str(speed)])
    disp(' ')
    state = state_new;
    rxns = rxns_new;
    trxns = trxns_new;
    rxn_rates_k(k,:) = rxn_rates; 
    gradient_k(k,:) = grdt;
    
    if  sum(abs(grdt./min_grad))<num_rxns
        min_grad = grdt;
        min_rate = rxn_rates;                    % save rxn rate with minimum gradient so far in "min_rates"
    end
  
    if max(abs(d_rate_exp)) < 5.e-03
        sprintf('Final step size: %5.3f \n', speed)
        disp('Paramter values:')
        disp(num2str(rxn_rates))
        break
    end
    
    if k == k_max
        disp(' ')
        disp('maximum steps reached.')
    end
end
diary off

% toc

% min_rate
% min_grad
% min_rates(i,:) =  rxn_rates;
% min_grads(i,:) = grdt;

% end
figure(1)
plot (rxn_rates_k);
xlabel('Steps');
ylabel('rate parameters');
figure(2)
plot (gradient_k);
xlabel('Steps');
ylabel('gradients');
