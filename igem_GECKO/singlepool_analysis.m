%% singlepool_analysis
% in this script the yeastgem model with extra constraints based on the
% GECKO method [1] is used to analyse the interplay between biomass growth and styrene production
%% initialize cobra toolbox
initCobraToolbox();
%% use cplex
changeCobraSolver('ibm_cplex');
%% model creation
% first the single-pool gecko model was exported to .mat file using a python script:
    % from geckopy import GeckoModel
    % import cobra, pandas, os
    % model = GeckoModel('single-pool')
    % cobra.io.save_matlab_model(model, "ec_model_singlepool.mat")
% the .mat model is then loaded in matlab:
modelFileName = 'ec_model_singlepool.mat';
cd 'C:\Users\Matthijs\Desktop\igem_optforce'
modelDirectory = pwd;
modelFileName= [modelDirectory filesep modelFileName];
ecModel = readCbModel(modelFileName);
fprintf("loading ecModel done\n"); 
%% add new reactions and metabolites
%YNL294C: PAL2, P45724, mw: 77.860kd, Kc = 3.2
%YDR539W: FDC1, Q03034, mw: 56.164kd, Kc = 4.6
%s_0456_c: co2
%s_1032_c: phenylalanine
%s_0419_c: amonium
%s_0794_c: H

%add TCA, proteine and styrene as metobalites
stModel = ecModel;
stModel = addMetabolite(stModel,'TCA_c','Trans-cinnamate','C9H8O2');
stModel = addMetabolite(stModel,'Styr_c','Styrene','C8H8');
stModel = addMetabolite(stModel,'prot_P45724_c', 'Phenylalanine ammonia lyase');   %PAL2: P45724, mw: 77.860kd 
stModel = addMetabolite(stModel,'prot_Q03034_c', 'Ferulic acid decarboxylase 1');   %FDC1: 'Q03034, mw: 56.164kd

%add reactions, in gecko models enzymes/proteins are added added as metabolites
%with 1/kcat as their stochiometric coefficent

stModel = addReaction(stModel, 'TCA_prod',...
'metaboliteList', {'s_0794_c' , 's_1032_c' , 'prot_P45724_c', 's_0419_c','TCA_c' },...
'stoichCoeffList', [-1; -1; -1/3.2; 1; 1;], 'reversible', false);
%kcat =	3.2 from; The Arabidopsis phenylalanine ammonia lyase gene family: kinetic characterization of the four PAL isoform
stModel = addReaction(stModel, 'TCAtoPhe',...
'metaboliteList', {'s_0794_c' , 's_1032_c' , 'prot_P45724_c', 's_0419_c','TCA_c' },...
'stoichCoeffList', [1; 1; 1/3.2; -1; -1;], 'reversible', false);
%FDC1
stModel = addReaction(stModel, 'Styr_prod',...
'metaboliteList', {'TCA_c' , 'prot_Q03034_c', 'Styr_c' , 's_0456_c' },...
'stoichCoeffList', [-1; -1/4.6; 1; 1;], 'reversible', false);
%kcat = 4.6 from: https://www.uniprot.org/uniprot/Q03034

stModel = addReaction(stModel, 'Styr_sink', 'metaboliteList', {'Styr_c'} ,...
'stoichCoeffList', -1, 'reversible', false);

%proteins are drawn from the protein pool
%PAL2
stModel = addReaction(stModel, 'draw_prot_P45724',...
'metaboliteList', {'prot_pool_c', 'prot_P45724_c'},...
'stoichCoeffList', [-77.860; 1; ], 'reversible', false);
%FDC1
stModel = addReaction(stModel, 'draw_prot_Q03034',...
'metaboliteList', {'prot_pool_c', 'prot_Q03034_c'},...
'stoichCoeffList', [-56.164; 1; ], 'reversible', false);


% add associated genes
 stModel = changeGeneAssociation(stModel, 'draw_prot_Q03034', 'YDR539W'); %FDC1
 stModel = changeGeneAssociation(stModel, 'draw_prot_P45724', 'YNL294C'); %PAL2

fprintf('added reaction and genes \n');
%% run FVA
%optimized for Styrene with growth constrained to 50% of max growth
stModel = changeObjective(stModel, 'r_2111');
growthRate = optimizeCbModel(stModel); 
fprintf('The maximum growth rate stModel is %1.8f \n', growthRate.f);
stModel = changeRxnBounds(stModel, 'r_2111',growthRate.f*0.1, 'l'); %constrain grwoth
stModel = changeObjective(stModel, 'Styr_sink');
styrProd = optimizeCbModel(stModel); 
fprintf('The maximum styrene production is %1.8f \n', styrProd.f);
[minflux, maxflux] = fluxVariability(stModel);

%% write models to excel file for easy reference
cd 'C:\Users\Matthijs\Desktop\igem_optforce'; 
 writeCbModel(model, 'format','xls', 'fileName', 'optmodel');

%% fba
stModel = changeObjective(stModel, 'r_2111');
growthRate = optimizeCbModel(stModel); 
fprintf('The maximum growth rate stModel is %1.8f \n', growthRate.f);
stModel = changeObjective(stModel, 'Styr_sink');
styrProd = optimizeCbModel(stModel); 
fprintf('The maximum styrene production is %1.8f \n', styrProd.f);
%% fba2
stModel = changeObjective(stModel, 'TCA_prod');
growthRate = optimizeCbModel(stModel); 
fprintf('The maximum growth rate stModel is %1.8f \n', growthRate.f);
stModel = changeObjective(stModel, 'Styr_sink');
styrProd = optimizeCbModel(stModel); 
fprintf('The maximum styrene production is %1.8f \n', styrProd.f);
%% make graph of growth vs Styrene
%optimize for growth whil gradually increasing styrene constraints
stModel = changeObjective(stModel, 'r_2111');
j=0;
GrowthGr= []; 
for i = 0:0.01:1
    j=j+1;
    stModel = changeRxnBounds(stModel, 'Styr_sink', styrProd.f*i , 'l');   %increase lowerbound of Styrene prod
    opt=optimizeCbModel(stModel);
    GrowthGr(j,1) = opt.obj;
    GrowthGr(j,2) = styrProd.f*i;
    GrowthGr(j,3) = opt.full(1573);                                     %Store resulting growthrate as percentage of max growthrate
    fprintf('%1.3f ', GrowthGr(j,1));
    
end

%reset constraints
stModel = changeRxnBounds(stModel, 'Styr_sink', 0 , 'l');

%% plot
plot(GrowthGr(:,1), '-');
title('Growthrate vs Styrene Production');
xlabel('Percentage of max Styrene');
ylabel('Percentage of max Growth');

%% references
% [1]Aung HW, Henry SA, Walker LP. Revising the Representation of Fatty Acid, Glycerolipid, and 
% Glycerophospholipid Metabolism in the Consensus Model of Yeast Metabolism. Ind Biotechnol 
% (New Rochelle N Y). 2013;9(4):215-228. doi:10.1089/ind.2013.0013

