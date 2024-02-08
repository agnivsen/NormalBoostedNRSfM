function [Data, nng, K, template] = SpecularDataGenerator(nFiles, gridSize, nNeighbors, reprojectionNoise, fractionVisible, showDebugPlots, randomNoise, fractionSpecular)

    %%% This script follows the code from [Perriollat et. al., 2013] for generating paper model. 
    %
    %  [Perriollat et. al., 2013]: Perriollat, Mathieu, and Adrien Bartoli. "A computational model of bounded developable surfaces with application to image?based three?dimensional reconstruction." Computer Animation and Virtual Worlds 24.5 (2013): 459-476.
    %                                             Code at: http://igt.ip.uca.fr/~ab/code_and_datasets/index.php

    close all;
    
    pathToStDRv1p0c = './St-DR_v1p0c/'; %%% On Windows
    
    fprintf('***\nThe code from [Perriollat et. al., 2013] is expected to be at:\n %s \n***\n', pathToStDRv1p0c);
    
    modelPath = strcat(pathToStDRv1p0c, '/model/');
    miscPath =  strcat(pathToStDRv1p0c, '/misc/');
    otherPath =  strcat(pathToStDRv1p0c, '/other/');
    
    addpath(pathToStDRv1p0c);
    addpath(modelPath);
    addpath(miscPath);
    addpath(otherPath);
    
    if ~exist('fractionSpecular', 'var')
        fractionSpecular = 0.05;
    end
    
    
    assert(numel(gridSize) == 2,'SpecularDataGenerator: for <gridSize>, expecting a [1 x 2] array of dimensions');
    assert(isequal(size(gridSize,1),1),'SpecularDataGenerator: for <gridSize>, expecting a [1 x 2] array of dimensions');
    assert(gridSize(1) >= 1 && isinteger(int8(gridSize(1))),'SpecularDataGenerator: <gridSize> dimension should be integer and >= 1');
    assert(gridSize(2) >= 1 && isinteger(int8(gridSize(2))),'SpecularDataGenerator: <gridSize> dimension should be integer and >= 1');
    assert(fractionVisible>0 && fractionVisible<=1, 'SpecularDataGenerator: the fraction of visible points need to be between (0,1]');    
    assert(fractionSpecular>0, 'SpecularDataGenerator: the fraction of specular points need to be between (0,1]');    
    assert(nFiles > 1 && isinteger(int8(nFiles)), 'SpecularDataGenerator: <nFiles> should be an integer and greater than 1');
    
    STRETCH_INCREMENT = 0.0; % each data will be stretched by this much additional amount compared to the previous data
    
    SPECULAR_NNG_CUTOFF = 1;
% % %     figure;
    
    if exist('nNeighbors','var')
        nK = nNeighbors;
    else
        nK = 3;
    end

    showDbg = false;
    
    if exist('showDebugPlots','var')
        if(showDebugPlots == true)
            showDbg = true;
        end    
    end
    
    if~exist('reprojectionNoise','var')
        reprojectionNoise = 0;
    end
    
    if exist('randomNoise','var')
        randomNoiseMult = randomNoise;
    else
        randomNoiseMult = 0;
    end
    
    fx = round(random_number_within_range(700,900,1));
    fy = round(random_number_within_range(700,900,1));
    
    IMG_HEIGHT = 1080*2;
    IMG_WIDTH = 1920*2;
    
    cx = IMG_WIDTH/2;
    cy = IMG_HEIGHT/2;
    
    K = [fx 0 cx; 0 fy cy; 0 0 1];
    
    
    % generate a paper parameter structure
    pp1.type = 'al_be'; % the id of the parameterization
    pp1.r = 1.4;%random_number_within_range(1,2,1); % width length ratio
%     arc1 = random_number_within_range(0, 1, 3);
%     arc2 = random_number_within_range(1, 2, 3);
%     pp.al = [arc1; arc2];
    pp1.al = [0.57 0.6374 0.9118 ; 1.1469 0.8737 1.111]; % arc lengths of the rulings
    
    
    pp1.be = [0 0 0]; % normalized bending angles

    % interpolate the rulings
    pp1 = paperParameterisationConversion(pp1,'al_be',5);

    % compute the paper
    paper = newPaperFast(pp1);

%     % display the flat paper
%     plotFlatPaper(paper);

    M_ = gridSize(1)*gridSize(2);
    expectedSpecularPoints = round(fractionSpecular*M_);
    neededAdditionalPoints = expectedSpecularPoints*3;
    toAddPerSide = ceil(sqrt(neededAdditionalPoints));
    gridSize = gridSize + toAddPerSide;
    
    M = gridSize(1)*gridSize(2);
    actualAdditionalPoints = M - M_;
    rejectIndex = randperm(M);
    rejectIndex = sort(rejectIndex(1:actualAdditionalPoints));

    % evaluate the 3D mesh
    templatePaper = paperMesh(paper,gridSize(1),gridSize(2));
    
    
    [template3D] = meshToPointcloud(templatePaper);
    template = template3D;
    
    
    
    tempdata.p(1).p = template3D.';
    tempdata.p(1).p(:, rejectIndex) = [];
    [nng] = getNeighborhoodDuplicated(tempdata,nK);
    
    
    
        % generate a paper parameter structure
    pp.type = 'al_be'; % the id of the parameterization
    pp.r = 1.4;%random_number_within_range(1,2,1); % width length ratio
%     arc1 = random_number_within_range(0, 1, 3);
%     arc2 = random_number_within_range(1, 2, 3);
%     pp.al = [arc1; arc2];
    pp.al = [0.57 0.6374 0.9118 ; 1.1469 0.8737 1.111]; % arc lengths of the rulings
   
    
    stretch = 1;
    
    for iFile = 1:nFiles
    
        bend = random_number_within_range(-10, 10, 3);
        pp.be = bend; % normalized bending angles

        % interpolate the rulings
        ppn = paperParameterisationConversion(pp,'al_be',5);

        % compute the paper
        paper = newPaperFast(ppn);

%         % display the flat paper
%         plotFlatPaper(paper);

        % evaluate the 3D mesh
        mesh3D = paperMesh(paper,gridSize(1),gridSize(2));
        
        [pointcloud] = meshToPointcloud(mesh3D);
        
        C = mean(pointcloud);
        
        pointcloud = pointcloud.' - C.';
        pointcloud = pointcloud.';
%         pointcloud = pointcloud + 0.01.*rand(size(pointcloud));
        
        rot1 = random_number_within_range(-1,1,1);
        rot2 = random_number_within_range(-1,1,1);
        rot3 = random_number_within_range(-90,90,1);
        
        rX = rotx(rot1); rY = roty(rot2); rZ = rotz(rot3);
        
        Rs = rX*rY*rZ;
        
        tx = random_number_within_range(-0.05,0.05,1);
        ty = random_number_within_range(-0.05,0.05,1);
        tz = random_number_within_range(0.8,1.0,1);
        
        t = [tx; ty; tz];
        
        pointcloud = Rs*pointcloud.' + t;
        
        pointcloud = pointcloud.';    
        
        [normalsPoly, normalsDiscrete] = getNormals(pointcloud, template3D(:, 1:2));
        normals = normalsPoly; % using the polynomial approximation of normals as the normal estimated from specularity
        % Normals estimated discretely as 'normalsDiscrete' remains a
        % secondary/unused variable, with possible future use
        
        pointcloud = pointcloud + (randomNoiseMult*rand(size(pointcloud)));
        
        pointcloud = (stretch .* pointcloud);
        stretch = stretch + STRETCH_INCREMENT;
        
        imgPoints = [];
        
        for ii=1:size(pointcloud,1)
            P_ = pointcloud(ii,:);
            px = P_(1)/P_(3)*fx + cx + reprojectionNoise*rand;
            py = P_(2)/P_(3)*fy + cy + reprojectionNoise*rand;
            p = [px py 1];
            imgPoints = [imgPoints;p];
        end
        
        if(showDbg)
            subplot(1,2,1);
            scatter(imgPoints(:,1), imgPoints(:,2), 'r*');
            xlim([0 IMG_WIDTH]);
            ylim([0 IMG_HEIGHT]);
            titStr = strcat('Frame no.', num2str(iFile));
            title(titStr);
        end

        Data.p(iFile).p = imgPoints.';
        Data.Pgth(iFile).P = pointcloud.';
        Data.Ngth(iFile).N = normals;
        Data.Ngth(iFile).N_discrete = normalsDiscrete;
        
        specular = normals(:, rejectIndex);
        specularPts = imgPoints(rejectIndex, :).';
        specularGT = pointcloud(rejectIndex, :).';
        sSpec = size(specular,2);
        acceptIndex = randperm(sSpec);
        acceptIndex = acceptIndex(1:expectedSpecularPoints);
        specular = specular(:, acceptIndex);
        specularPts = specularPts(:, acceptIndex);
        specularGT = specularGT(:, acceptIndex);
        
        Data.specularity(iFile).normals = specular;
        Data.specularity(iFile).pts2d = specularPts;
        Data.specularity(iFile).pts3dGT = specularGT;
        
        Data.p(iFile).p(:, rejectIndex) = [];
        Data.Pgth(iFile).P(:, rejectIndex) = [];
        Data.Ngth(iFile).N(:, rejectIndex) = [];
        Data.Ngth(iFile).N_discrete(:, rejectIndex) = [];
        
        if(nNeighbors > 2)
            snC = 1;%max(3, round(nNeighbors/3));
            arbitraryNeighborDist = EuclideanDistance(Data.p(iFile).p(:, 1),Data.p(iFile).p(:, nng(1,snC)));
        else
            arbitraryNeighborDist = EuclideanDistance(Data.p(iFile).p(:, 1),Data.p(iFile).p(:, nng(1,1)));
        end
        
        for iP = 1:M_
            p = Data.p(iFile).p(:,iP);
            specIndex = 1;
            for iS = 1:size(specularPts,2)
                s = specularPts(:,iS);
                d = EuclideanDistance(p,s);
                if(d < arbitraryNeighborDist)
                    nngSpecular{iP, specIndex} = iS;
                    specIndex = specIndex + 1;
                end
            end
        end
        
        nngSpecular{iP, specIndex} = [];
        
        if exist('nngSpecular', 'var')
            tf = cellfun('isempty',nngSpecular); % true for empty cells
            nngSpecular(tf) = {0};

            nngSpecular = cell2mat( nngSpecular );
            
            if(size(nngSpecular,2)>SPECULAR_NNG_CUTOFF)
                nngSpecular(:, SPECULAR_NNG_CUTOFF+1:end) = [];
            end
        else
            nngSpecular = 0;
        end
        
        Data.specularNNG(iFile).S = nngSpecular;
        
        clear nngSpecular;
        
        if(showDbg)
            subplot(1,2,2);
            pcshow(pointcloud,'MarkerSize', 70); hold on;
            set(gcf,'color','w');
            set(gca,'color','w');
            pause(1);
        end
    
    end
    
    nFeats = M_;
    
    visibilityMatrix = ones(nFiles, nFeats);
    
    if(fractionVisible < 1)
        visNum = round((1 - fractionVisible) * nFeats);
        for ii = 2:nFiles % because the 1st image always need to have all keypoints visible
            visList = randperm(nFeats,visNum);
            visibilityMatrix(ii,visList) = 0;
        end
    end
    
    Data.v = visibilityMatrix;
    
    if(showDbg)
        close all;
    end

end

function [normals, normalsDiscreteApproximation] = getNormals(pointcloud, template)
    assert(size(pointcloud,2) == 3, 'Pointcloud must be a [N x 3] array!');
    assert(size(template,2) == 2, 'Template must be a [N x 2] array!');
    assert(size(template,1) == size(pointcloud,1), 'Template and pointcloud must be pf same length');
    
    EPS = 10^(-12);
    
    X= template(:,1).' + EPS; Y= template(:,2).' + EPS;
    
    O1 = ones(1, size(template,1));
    Z0 = zeros(1, size(template,1));
    
    monomialBasis = [X.^4; Y.^4; Y.*X.^3; X.*Y.^3; X.^2.*Y.^2; X.^3; Y.^3; X.^2.*Y; X.*Y.^2; X.^2; Y.^2; X.*Y; X; Y; O1];
    
    monomialCoeff = pointcloud.'*pinv(monomialBasis);
    
    xDeriv = [4*X.^3; Z0; 3.*X.^2.*Y; Y.^3; 2.*X.*Y.^2; 3.*X.^2; Z0; 2.*X.*Y; Y.^2; 2.*X; Z0; Y; O1; Z0; Z0];
    yDeriv = [Z0; 4*Y.^3; X.^3; 3.*Y.^2.*X; 2.*Y.*X.^2; Z0; 3.*Y.^2; X.^2; 2.*X.*Y; Z0; 2.*Y; X; Z0; O1; Z0];
    
    Px = monomialCoeff*xDeriv;
    Py = monomialCoeff*yDeriv;
    
    normals = zeros(3, size(template,1));
    
    normalsDiscreteApproximation = pcnormals(pointCloud(pointcloud));
    normalsDiscreteApproximation = normalsDiscreteApproximation.';
    
    for iP = 1:size(template,1)
        P1 = Px(:, iP);
        P2 = Py(:, iP);
        
        P1 = P1./norm(P1);
        P2 = P2./norm(P2);
        
        N_ = cross(P1, P2);
        
        normals(:, iP) = N_./norm(N_);
        
        Pg = pointcloud(iP,:).';
        Ndiscrete = normalsDiscreteApproximation(:, iP);
        
        if(dot(Pg, Ndiscrete) < 0)
            normalsDiscreteApproximation(:, iP) = -normalsDiscreteApproximation(:, iP);
        end
        
    end
end

function [pointcloud] = meshToPointcloud(mesh3D)

    numPts = size(mesh3D,1)*size(mesh3D,2);
    
    pointcloud = zeros(numPts, 3);
    
    index = 1;
    
    for ii = 1: size(mesh3D,1)
        for jj = 1:size(mesh3D,2)
            pointcloud(index,:) = [mesh3D(ii,jj, 1) mesh3D(ii,jj, 2) mesh3D(ii,jj, 3)];
            index = index + 1;
        end
    end
    
    
end
