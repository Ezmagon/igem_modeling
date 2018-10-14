%% visualization.m
%this script draws a paint4net map of (part of a cell's) metabolism
%and overlays this with the MUST sets found by the OptForce algorithm.

%% load MUST sets
%change to directory of MUST sets (found using the OptForce algorithm)
mustU = load('C:\Users\Matthijs\Desktop\igem_optforce\TestOptForce5M\InputsOptForce\mustU.mat');
mustU = mustU.mustU;
mustL = load('C:\Users\Matthijs\Desktop\igem_optforce\TestOptForce5M\InputsOptForce\mustL.mat');
mustL = mustL.mustL;

%% init cobra toolbox
addpath(genpath('cobratoolbox'));
addpath(genpath('yeast-GEM-8.1.2/yeast-GEM-8.1.2'));
initCobraToolbox();
changeCobraSolver('ibm_cplex', 'ALL');
%% loading model
fprintf('loading model..')
model=loadYeastModel(); %function defined in yeastGEM. 
fprintf('loading model done\n');
       
%% make a paint4net map of subcompartiment metabolism
%here we make a paint4net map of Phenylalanine biosynthesis and metabolism, though
%any other compartment can be specified here

rxns=findRxnsFromSubSystem(model, 'sce00400  Phenylalanine, tyrosine and tryptophan biosynthesis');
rxns2=findRxnsFromSubSystem(model, 'sce00360  Phenylalanine metabolism');
rxns=union(rxns,rxns2);
[involvedMets,deadEnds]=draw_by_rxn(model,rxns, 'true', 'struc');

%IMPORTANT! the resulsting figure must be manually saved to workspace using
%file->export to workspace->name is 'testvis'

%% colour MUST sets
backup=testvis; 
targetmet= 's_1032[c]';  %phenylanaline, change to overproduction target
fprintf("colouring mustsets....");

%first make ll edges everything blue
for i=1:length(testvis.Edges)
    testvis.Edges(i).LineColor = [0 0 0.6];
end
%colour reaction nodes light blue and mets offwhite
for i=1:length(testvis.Nodes)
    if (isempty(strfind(testvis.nodes(i).ID, 'r_'))==0)
            testvis.Nodes(i).Color = [0.8 0.9 1];
    else
        testvis.Nodes(i).Color = [1 1 0.9];
    end
end

%then colour mustL reactions red
for i=1:length(testvis.Nodes)
    for j=1:length(mustL)
        if (isempty(strfind(testvis.nodes(i).ID, mustL(j)))==0)
            testvis.Nodes(i).Color = [1 0 0];
        end
    end
end

%and the corresponding edges darker red 
for i=1:length(testvis.Edges)
    for j=1:length(mustL)
        if (isempty(strfind(testvis.Edges(i).ID, mustL(j)))==0)
            testvis.Edges(i).LineColor = [0.6 0 0];
        end
    end
end

%then colour mustU green
for i=1:length(testvis.Nodes)
    for j=1:length(mustU)
        if (isempty(strfind(testvis.nodes(i).ID, mustU(j)))==0)
            testvis.Nodes(i).Color = [0 0.8 0];
        end
    end
end
%and the corresponding edges darker green
for i=1:length(testvis.Edges)
    for j=1:length(mustU)
        if (isempty(strfind(testvis.Edges(i).ID, mustU(j)))==0)
            testvis.Edges(i).LineColor = [0 0.6 0];
        end
    end
end

%finally colour target met bright orange
for i=1:length(testvis.Nodes)
    if (isempty(strfind(testvis.nodes(i).ID, targetmet))==0)
            testvis.Nodes(i).Color = [1 0.5 0];
    end
end

fprintf(" colouring done \n");
view(testvis)
%% tidy up names
for i=1:length(testvis.Nodes)
    testvis.Nodes(i).Label=erase(testvis.Nodes(i).Label, 'Name:');
    testvis.Nodes(i).Label=erase(testvis.Nodes(i).Label, ' [cytoplasm]');
end

%% make sure reaction labels show genes and met nodes are wider
%open the testvis biograph object and select show "Label" to see this!

for i=1:length(testvis.Nodes)
    testvis.Nodes(i).ID=erase(testvis.Nodes(i).ID, ' (x)');
     if (contains(testvis.nodes(i).ID, 'r_'))
            gen=findGenesFromRxns(model, testvis.nodes(i).ID);
            gen=gen{1};
            genN=model.geneNames(find(contains(model.genes,gen{1})));
            if (length(gen)>1)
                 for j=2:length(gen)
                      genT=model.geneNames(find(contains(model.genes,gen{j})));
                      genN=strcat(genN, {'  '}, genT);
                 end
            end
            testvis.nodes(i).Label=genN{:};
     else
         testvis.nodes(i).size(1)=testvis.nodes(i).size(1)*2;
     end
end
%% resize function
scale=2;
for i=1:length(testvis.Nodes)
    testvis.nodes(i).size(1)=testvis.nodes(i).size(1)*scale;
end
    view(testvis)

