AKARNAcLoadingCell =...
    {'\\nasquatch\data\2p\stephen\SZ309\FLIM\200128_SZ309', '200128_SZ309_slice1_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ309\FLIM\200128_SZ309', '200128_SZ309_slice2_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ310\FLIM\200131_SZ310', '200131_SZ310_slice1_output.mat'};
AKARNAcLoadingCell = AKARNAcLoadingCell([1,3],:);

AKARMPOALoadingCell = ...
    {'\\nasquatch\data\2p\stephen\SZ426\200627_SZ426\FLIM\Output\', '200627_SZ426_slice1_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ426\200627_SZ426\FLIM\Output\', '200627_SZ426_slice2_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ427\200627_SZ427\FLIM\Output\', '200627_SZ427_slice1_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ427\200627_SZ427\FLIM\Output\', '200627_SZ427_slice2_output.mat'};

cADDisNAcLoadingCell =...
    {'\\nasquatch\data\2p\stephen\SZ311\FLIM\200129_SZ311', '200129_SZ311_slice1_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ311\FLIM\200129_SZ311', '200129_SZ311_slice2_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ424\200626_SZ424\FLIM\Output', '200626_SZ424_slice1_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ424\200626_SZ424\FLIM\Output', '200626_SZ424_slice2_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ424\200626_SZ425\FLIM\Output', '200626_SZ425_slice1_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ424\200626_SZ425\FLIM\Output', '200626_SZ425_slice2_output.mat'};
cADDisNAcLoadingCell = cADDisNAcLoadingCell(1,:);

cADDisMPOALoadingCell =...
    {'\\nasquatch\data\2p\stephen\SZ317\FLIM\200130_SZ317', '200130_SZ317_slice1_output.mat';...
    '\\nasquatch\data\2p\stephen\SZ312\200219_SZ312\FLIM\Output', '200219_SZ312_slice1_output.mat'};

%% Initialize data
nbrightcells = 70;

nAKARNAcExpts = size(AKARNAcLoadingCell, 1);
ncADDisNAcExpts = size(cADDisNAcLoadingCell, 1);
ncADDisMPOAExpts = size(cADDisMPOALoadingCell, 1);

AKARNAcData = zeros(nbrightcells * nAKARNAcExpts, 2);
cADDisNAcData = zeros(nbrightcells * ncADDisNAcExpts, 2);
cADDisMPOAData = zeros(nbrightcells * ncADDisMPOAExpts, 2);

%% Extract data AKAR NAc Expts
for i = 1 : nAKARNAcExpts
    loaded = load(fullfile(AKARNAcLoadingCell{i,1}, AKARNAcLoadingCell{i,2}));
    [~ ,brightnessrank] = sort(loaded.OutputStruct.Photons_start, 'descend');
    
    AKARNAcData((i-1) * nbrightcells + 1: i * nbrightcells, 1) = ...
        loaded.OutputStruct.Photonsdff_values(brightnessrank(1:nbrightcells)) - 1;
    AKARNAcData((i-1) * nbrightcells + 1: i * nbrightcells, 2) = ...
        loaded.OutputStruct.Tmdff_values(brightnessrank(1:nbrightcells));
end

%% Extract data cADDis NAc Expts
for i = 1 : ncADDisNAcExpts
    loaded = load(fullfile(cADDisNAcLoadingCell{i,1}, cADDisNAcLoadingCell{i,2}));
    [~ ,brightnessrank] = sort(loaded.OutputStruct.Photons_start, 'descend');
    
    cADDisNAcData((i-1) * nbrightcells + 1: i * nbrightcells, 1) = ...
        loaded.OutputStruct.Photonsdff_values(brightnessrank(1:nbrightcells)) - 1;
    cADDisNAcData((i-1) * nbrightcells + 1: i * nbrightcells, 2) = ...
        loaded.OutputStruct.Tmdff_values(brightnessrank(1:nbrightcells));
end

%% Extract data cADDis MPOA Expts
for i = 1 : ncADDisMPOAExpts
    loaded = load(fullfile(cADDisMPOALoadingCell{i,1}, cADDisMPOALoadingCell{i,2}));
    [~ ,brightnessrank] = sort(loaded.OutputStruct.Photons_start, 'descend');
    
    cADDisMPOAData((i-1) * nbrightcells + 1: i * nbrightcells, 1) = ...
        loaded.OutputStruct.Photonsdff_values(brightnessrank(1:nbrightcells)) - 1;
    cADDisMPOAData((i-1) * nbrightcells + 1: i * nbrightcells, 2) = ...
        loaded.OutputStruct.Tmdff_values(brightnessrank(1:nbrightcells));
end
%% Plots
figure('Position', [50,200,1600, 400])
xlims = [-0.6 0.1];
ylims = [-300 50]; 

subplot(1,3,1)
hold on
plot(AKARNAcData(:,1),AKARNAcData(:,2), '.')
plot([median(AKARNAcData(:,1)), median(AKARNAcData(:,1))], ylims, 'r-');
plot(xlims, [median(AKARNAcData(:,2)), median(AKARNAcData(:,2))], 'r-');
hold off
xlabel('dF/F0')
ylabel('Lifetime change (ps, end - start)')
title('AKAR @NAc Brightest cells')
% xlim([-0.3 0.1]);
ylim(ylims);

subplot(1,3,2)
hold on
plot(cADDisNAcData(:,1),cADDisNAcData(:,2), '.')
plot([median(cADDisNAcData(:,1)), median(cADDisNAcData(:,1))], ylims, 'r-');
plot(xlims, [median(cADDisNAcData(:,2)), median(cADDisNAcData(:,2))], 'r-');
hold off
xlabel('dF/F0')
ylabel('Lifetime change (ps, end - start)')
title('cADDis @NAc Brightest cells')
xlim(xlims);
% ylim([-100 50]);

subplot(1,3,3)
hold on
plot(cADDisMPOAData(:,1),cADDisMPOAData(:,2), '.')
plot([median(cADDisMPOAData(:,1)), median(cADDisMPOAData(:,1))], ylims, 'r-');
plot(xlims, [median(cADDisMPOAData(:,2)), median(cADDisMPOAData(:,2))], 'r-');
hold off
xlabel('dF/F0')
ylabel('Lifetime change (ps, end - start)')
title('cADDis @MOPA Brightest cells')
xlim(xlims);
ylim(ylims);