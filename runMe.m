
close all; clear all;
clear vars;

DATA_OPTION = {'Isometric-paper-generator', 'Pink-paper-mobile'};
data = DATA_OPTION{1};

addpath('Data\');
addpath('Utils\');
addpath('Utils\St-DR_v1p0c\');
addpath('D:\code\cvx-w64\cvx\'); % path to CVX installation

nPx = 3;
nPy = 3;
maxInnerIterationForTesting = 20; % recommended value = 100
maxOuterIterationForTesting = 10; % recommended value = 100
numImages = 6; % recommended >7
noiseOnPixels = 0.0;
visibilityFraction = 1.0; % in (0, 1], the lower the value, the less points are actually visible
shouldVisualizeSimulationImages = true;
showDebugPlotsWhileSimulating = false;
extensibleNoise = 0.0;
fractionSpecular = 0.9;

OP_PATH = '/tmp/';

if(strcmp(data, 'Isometric-paper-generator'))
    %%%%%%%%%%%%%%% OPTION 1 %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Isometric paper model - initialisation
        nK = 6;
        [Data, nng, K] = SpecularDataGenerator(numImages, [nPx nPy], nK, noiseOnPixels, visibilityFraction, showDebugPlotsWhileSimulating, extensibleNoise, fractionSpecular);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(strcmp(data, 'Pink-paper-mobile'))
    %%%%%%%%%%%%%%% OPTION 2 %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Data = load('Data\PinkPaper.mat');
    nK = 15;
    intrinsics = load('Data\intrinsicsPinkPaper.txt');
    K = intrinsics;
    [nng] = getNeighborhoodDuplicated(Data,nK);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

[reconstructionMDHSpecular] = NRSfM(Data, nng, K, 'SpecularNormals', {Data.specularity.normals},...
    'SpecularPoints', {Data.specularity.pts2d}, 'SpecularNNG', ...
    Data.specularNNG, 'EPSILON', 10^(-2), 'Solver', 'sdpt3');


f = figure;
f.Position = [50 50 600 600];
for ii = 1:size(Data.p,2)
    P = reconstructionMDHSpecular.reconstruction{ii};
    pcshow(P,'r', 'MarkerSize',90);
    hold on; xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on; axis on; axis equal;
    hold off; pause(0.1);
    pause(0.333);
end

