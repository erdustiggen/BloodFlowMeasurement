% faster implementation of vectorDoppler using precalculated aliasing
% residuals

function [vectorEst, aliasPat, aliasMat, minRes] = vectorDoppler(R1,p, aliasMat)

if(isfield(p,'useGPU'))
    useGPU = p.useGPU; 
else
    useGPU = 0;
end

if isfield(p,'doAliasingCorrection')
    if p.doAliasingCorrection
        
        switch p.application
            
            case 'carotid'
             aliasLimit = 3;
        
            case 'neo'
             aliasLimit = 3;
        end
    else
        aliasLimit = 1;
    end
else
    aliasLimit = 3;
end

sizeR1 = size( R1 );
% [txangles,rxangles] = meshgrid(p.angles,p.rxangles);
% txangles = txangles.';
% rxangles = rxangles.';
% txangles = txangles(:);
% rxangles = rxangles(:);
rxangles = p.rxVals';
txangles = p.txVals';


%   vectorEst = single(zeros(size(R1,1)*size(R1,2),2));
fMat = (reshape(angle(R1),[size(R1,1)*size(R1,2),size(R1,3)])*p.PRF/(2*pi) ).';
%aMat = [-sin(p.angles) (1+cos(p.angles))]./2; % When rx steering angle is always zero 
aMat = [-sin(txangles)-sin(rxangles) cos(txangles)+cos(rxangles)]./2; % Valid also when using multiple steering angles on tx and rx

pseudoInv = pinv(aMat);
nrAngles = sizeR1(3);

if nargin < 3,
    nrCases = aliasLimit^nrAngles; %including no aliasing, but not aliasing on all angles
    sprintf('Number of aliasing cases/patterns is %i',nrCases);
    aliasMat = zeros(nrAngles, nrCases );
    for patternCode = 0:nrCases-1,
    %     str = char( dec2base( patternCode, aliasLimit) );
    %     tail = str2num( str(:) );
        pattern = zeros(nrAngles, 1);
        tail = dec2baseJorg( patternCode, aliasLimit);
        pattern(end-length(tail)+1:end) = tail;
        aliasMat(:, patternCode+1) = pattern;
    end
end
validinds = find( any( ~aliasMat, 1) );
aliasMat = aliasMat(:, validinds);
nrCases = size( aliasMat,2);
nrPoints = size(fMat,2);

gRes = ( ( aMat*pseudoInv-eye(nrAngles) )*aliasMat*p.PRF ).';
fRes = ( aMat*pseudoInv-eye(nrAngles) )*fMat;

temp = zeros(nrCases, size(R1,1)*size(R1,2));                                
% GPU        
if(useGPU)
    temp_g = gpuArray(temp);        
    gRes_g = gpuArray(gRes);
    fRes_g = gpuArray(fRes);        
    for i = 1:nrAngles,                
        temp_g = temp_g + bsxfun(@plus, fRes_g(i,:), gRes_g(:,i)).^2;
    end        
    temp = gather(temp_g);
else
    for i = 1:nrAngles,
       %temp = temp + ( repmat( fRes(i,:), [size(gRes,1) 1])+repmat( gRes(:,i), [1 size(fRes,2)]) ).^2;
       temp = temp + bsxfun(@plus, fRes(i,:), gRes(:,i)).^2;
    end            
end        
[minRes, aliasPat] = min( temp, [], 1);
vectorEst = ( pseudoInv*( fMat+aliasMat(:,aliasPat)*p.PRF ) ).';

%         temp2 = double( zeros(2, 2059, nrPoints) );
%         parpool(4);
%         parfor kk = 1:nrPoints
%             bigfMat = repmat( fMat(:,kk), [1 size( aliasMat, 2)] ) + aliasMat*p.PRF;
%             temp2(:,:,kk) = pseudoInv*bigfMat;
%             res2(kk,:) = sum( ( aMat*temp2(:,:,kk) - bigfMat ).^2, 1 );
%         end
%         delete(gcp);
%         toc
%         
% %       Debug part
% %         [minres, minind] = min( res(4,:) );
% %         for i = 1:size(temp,3),
% %             col = 'b';
% %             if i == minind,
% %                 col = 'r';
% %                 atand(temp(1,2,i)./temp(1,1,i) )
% %             end
% %             ang(i) = atand(temp(1,2,i)./temp(1,1,i) );
% %             figure(100), hold on, plot([0 temp(1,1,i)], [0 temp(1,2,i)], col );
% %         end
%  
%         for kk = 1:size(fMat,2),
%             [minres2, minind2] = min( res2(kk,:) );
%             aliasPat2(kk) = minind2;
%             minRes2tab(kk) = minres2;
%             vectorEst2(kk,:) = temp2(:,minind2,kk);
%         end

% figure(); plot(fMat,'o'); title('All frequency estimates from fMat');





end