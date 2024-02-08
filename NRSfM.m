% NRSfM is the method that reconstructs deforming shapes, with or without
% the presence of specular highlights in the scene, following the method
% described in equation (2) of our IEEE ICRA-2024 paper, "Using Specularities to Boost 
% Non-Rigid Structure-from-Motion"
%
function [resultNRSfM] = NRSfM(Data, nng, K, varargin)

    p = inputParser;
    expectedSolvers = {'mosek', 'sdpt3', 'yalmip', 'sedumi'};
    defaultEpsilon = 0.00;
    defaultSolver = expectedSolvers{1};
    defaultScaleCompensation = true;

    addRequired(p,'Data');
    addRequired(p,'nng');
    addRequired(p,'K');
    addParameter(p,'Solver',defaultSolver,@(x) any(validatestring(x,expectedSolvers)));
    addParameter(p,'EPSILON',defaultEpsilon);
    addParameter(p,'RescaleData',defaultScaleCompensation);
    addParameter(p,'SpecularNormals',[]);
    addParameter(p,'SpecularPoints',[]);
    addParameter(p,'SpecularNNG',[]);

    parse(p,Data, nng, K, varargin{:});     
    edI = EdgeIndexer(nng);    
    params = p.Results;

    [nFiles, nPts, kN, hasSpecularity] = checkDataSanity(params);
    [directionVector] = getDirectionVectors(Data, K);
    
    
    cvx_begin
    cvx_precision low
    if(strcmp(params.Solver,expectedSolvers{1}))
        cvx_solver mosek
    elseif(strcmp(params.Solver,expectedSolvers{2}))
        cvx_solver sdpt3
    elseif(strcmp(params.Solver,expectedSolvers{3}))
        cvx_solver yalmip
    elseif(strcmp(params.Solver,expectedSolvers{4}))
        cvx_solver sedumi
    end

   
    variable geod(edI.numberOfEdges);
    
    variable deltaD(nFiles,nPts);
    
    resSpecular = [];

    for iF = 1:nFiles
        fprintf('Modelling constraints for <strong>image #%d</strong>\n', iF);
        textprogressbar('Progress:   ');
        for iP = 1:nPts
            if(Data.v(iF,iP)~=0)
                vecJ = directionVector{iF}(:,iP);
                delJ = deltaD(iF, iP);
                delJ >= 0;
                
                P = delJ.*vecJ;
                
                textprogressbar((round((iP/nPts)*100)));
                
                if(hasSpecularity) % Specular NRSfM
                    specularNeighbors_j = params.SpecularNNG(iF).S(iP,:);
                end
                
                for iN = 1:(kN)
                    nN1 = nng(iP, iN);
                    if(Data.v(iF,nN1)~=0)
                        vecQ = directionVector{iF}(:,nN1);                
                        delQ = deltaD(iF, nN1);
                        edgeIndex = edI.getIndex(iP, nN1);
                        Q = delQ.*vecQ;
                        
                        if(hasSpecularity) % Specular NRSfM
                            for iS = 1:numel(specularNeighbors_j)
                                iSs = specularNeighbors_j(iS);
                                if(iSs>0)
                                    specNorm =  params.SpecularNormals{iF}(:,iSs);
                                    
                                    dotRes = dot(P-Q, specNorm);
                                    resSpecular = [resSpecular abs(dotRes)];
                                end
                            end
                        end  
                        
                        eucl_PQ = norm(P - Q, 2);
                        
                        eucl_PQ <= geod(edgeIndex); 
                        
                        if(iF == 1)
                            geod(edgeIndex) >=0;
                        end
                    end
                end
            end
        end
        textprogressbar('done.');
        
    end
    
    sum(sum(geod)) == 1;
    
    if(hasSpecularity)
        minimize(-sum(sum(deltaD)) + (sum(resSpecular)));
    else
        minimize(-sum(sum(deltaD)));
    end


    cvx_end
    
    
    
    
    fprintf('\n---------------------------------------------------------------------------\n');

    depthOfPoints = zeros(nFiles, nPts);

    for iF = 1:nFiles
        depthOfPoints(iF,:) = deltaD(iF,:);
        depthOfPointsFormat2{iF} = deltaD(iF,:);
    end

    geodesicMatrix = geod;

    
    [reconstruction, scaleList] = getReconstruction(Data, depthOfPoints, directionVector, nFiles);
    [reconstructionUnscaled, ~] = getReconstructionUnscaled(Data, depthOfPoints, directionVector, nFiles); 
     

    resultNRSfM.status = cvx_status;
    resultNRSfM.reconstruction = reconstruction;
    resultNRSfM.reconstructionUnscaled = reconstructionUnscaled;
    resultNRSfM.scaleList = scaleList;
    resultNRSfM.geodesicMatrix = geodesicMatrix;
    resultNRSfM.depthOfPoints = depthOfPointsFormat2;
    resultNRSfM.delta = deltaD;
    

end

function [nFiles, nPts, kN, hasSpecularity] = checkDataSanity(params)
    nFiles = size(params.Data.Pgth,2);
    nPts = size(params.Data.Pgth(1).P,2);
    kN = size(params.nng,2);
    
    hasSpecularity = false;
    if(size(params.SpecularNNG) >= 1)
        hasSpecularity = true;
    end
    
    fprintf('\n\n************************************************************************\n');
    timestmp = strrep(datestr(datetime('now')),' ','_');
    fprintf('Running at (timestamp) = [%s]\n', timestmp);
    fprintf('*   ==: Solver = <strong>%s</strong>\n',params.Solver);
    fprintf('*   ==: Epsilon (lower bounds on +ve values) = <strong>%d</strong>\n',params.EPSILON);                                                                                                                                                        
    fprintf('*   ==: No. of neighbors = <strong>%d</strong>\n',kN);
    fprintf('*   ==: No. of images = <strong>%d</strong>\n',nFiles);
    fprintf('*   ==: No. of points/image = <strong>%d</strong>\n',nPts);
    if(hasSpecularity)
        cprintf('red', '*   ==: Has specularity information\n');
    end
    fprintf('************************************************************************\n\n');
    
end
