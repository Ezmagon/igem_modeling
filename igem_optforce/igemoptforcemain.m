%% iGEM OptForce
% I edited the OptForce  tutorial script to use the yeastgem model [1],[2]
% the script is used to calculate interventions leading to overproducing
% phenylanaline

%% init CobraToolbox
addpath(genpath('cobratoolbox'));
addpath(genpath('yeast-GEM-8.1.2/yeast-GEM-8.1.2'));
initCobraToolbox();
%% loading model
changeCobraSolver('ibm_cplex', 'ALL');
model=loadYeastModel(); %function defined in yeastGEM. 
fprintf('loading model done');

%% Define constraints for both wild-type and mutant strain
%the wildtype's biomass flux (r_4041) is constrained to near it's max while the
%mutants phenylanaline flux (r_1903) is constrained to near it's max.
constrWT = struct('rxnList', {{'r_4041'}}, 'rxnValues', 0.07, 'rxnBoundType', 'b');
constrMT = struct('rxnList', {{'r_1903'}}, 'rxnValues', 0.5,  'rxnBoundType', 'b');

fprintf('constraints set \n');
%%  Flux Variability Analysis
fprintf('running FVA...');
[minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ~, ~] = FVAOptForce(model, ...
                                                                     constrWT, constrMT);
disp([minFluxesW, maxFluxesW, minFluxesM, maxFluxesM]);
fprintf('FVA done');

%% Find Must Sets

runID = 'TestOptForce5M';

% The structure of the folders will be as follows:
% 
% |??? CurrentFolder|
% 
% ||   ??? run ID|
% 
% ||   |   ??? Inputs|
% 
% ||   |   ??? Outputs|
% 
%% We define constraints.
 
constrOpt = struct('rxnList', {{'r_1903'}}, 'values', 0.5);

%% MUSTL
fprintf('running mustL...');
[mustLSet, pos_mustL] = findMustL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                  'runID', runID, 'outputFolder', 'OutputsFindMustL', ...
                                  'outputFileName', 'MustL' , 'printExcel', 1, 'printText', 1, ...
                                  'printReport', 1, 'keepInputs', 1, 'printLevel', 3);
                              
fprintf('MustL done\n');

disp(mustLSet)
%% MUST U
fprintf('running mustU...');
[mustUSet, pos_mustU] = findMustU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                  'runID', runID, 'outputFolder', 'OutputsFindMustU', ...
                                  'outputFileName', 'MustU' , 'printExcel', 1, 'printText', 1, ...
                                  'printReport', 1, 'keepInputs', 1, 'printLevel', 3);
                              
                                    
fprintf('MustU done\n');

disp(mustUSet)
%% second order sets

% First, we define the reactions that will be excluded from the analysis. 
% As tt is suggested to include in this list the reactions found in the previous 
% step as well as exchange reactions we do this

[SelExc, SelUpt] = findExcRxns(model);
idx = find(SelExc);
exchangeRxns = model.rxns(idx);
excludedRxns = unique([mustUSet; mustLSet; exchangeRxns]);

%% MUST UU
fprintf('running mustUU...');
[mustUU, pos_mustUU, mustUU_linear, pos_mustUU_linear] = ...
    findMustUU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUU', 'outputFileName', 'MustUU', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);
                                               
fprintf('MustUU done\n');

disp(mustUU);

%% MUST LL
fprintf('running mustLL...');
[mustLL, pos_mustLL, mustLL_linear, pos_mustLL_linear] = ...
    findMustLL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustLL', 'outputFileName', 'MustLL', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);
fprintf('MustLL done\n');

disp(mustLL);
%%  OptForce
% We define constraints and we define |K| the number of interventions allowed, 
% |nSets| the maximum number of sets to find, and |targetRxn| the reaction producing 
% the metabolite of interest (in this case, succinate). 
% 
% Additionally, we define the |mustU| set as the union of the reactions that 
% must be upregulated in both first and second order must sets; and |mustL| set 
% as the union of the reactions that must be downregulated in both first and second 
% order must sets .k = 2;


fprintf("running optforce...");
constrOpt = struct('rxnList', {{'r_4041'}}, 'values', 0);
mustU = unique(union(mustUSet, mustUU));
mustL = unique(union(mustLSet, mustLL));
targetRxn = 'r_1903';
biomassRxn = 'r_4041';
k = 1;
nSets = 20;

[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
    optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
             minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
             'k', k, 'nSets', nSets, 'constrOpt', constrOpt, ...
             'runID', runID, 'outputFolder', 'OutputsOptForce', ...
             'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
             'printReport', 1, 'keepInputs', 1, 'verbose', 1);
         
         fprintf("done\n");

disp(optForceSets)
%% second run excluding reactions found in previous run as well as echange reactions
[SelExc, SelUpt] = findExcRxns(model);
idx = find(SelExc);
exchangeRxns = model.rxns(idx);
excludedRxns = unique([exchangeRxns; 'r_1215'; 'r_0148']);
mustU=setdiff(mustU,excludedRxns);
mustL=setdiff(mustL,excludedRxns);
k = 2;
nSets = 5;
constrOpt = struct('rxnList', {{'r_4041'}}, 'values', 0);

runID = 'Test_local6';


[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
    optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
             minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
             'k', k, 'nSets', nSets, 'constrOpt', constrOpt, ...
             'runID', runID, 'outputFolder', 'OutputsOptForce', ...
             'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
             'printReport', 1, 'keepInputs', 1, 'printLevel', 1);
         
         fprintf('done\n');
disp(optForceSets)
%% Acknowledgments
% I would to thanks to the research group of Costas D. Maranas who were
% willing to answer all questions I had while running this algorithm
%% References
% [1] Ranganathan S, Suthers PF, Maranas CD (2010) OptForce: An Optimization 
% Procedure for Identifying All Genetic Manipulations Leading to Targeted Overproductions. 
% PLOS Computational Biology 6(4): e1000744. https://doi.org/10.1371/journal.pcbi.1000744.
%
% [2]Aung HW, Henry SA, Walker LP. Revising the Representation of Fatty Acid, Glycerolipid, and 
% Glycerophospholipid Metabolism in the Consensus Model of Yeast Metabolism. Ind Biotechnol 
% (New Rochelle N Y). 2013;9(4):215-228. doi:10.1089/ind.2013.0013
