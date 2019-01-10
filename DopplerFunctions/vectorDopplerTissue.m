% no aliasing, and input is velocities, not phase

% faster implementation of vectorDoppler using precalculated aliasing
% residuals


function [vectorEst, residuals] = vectorDopplerTissue(vels,p)

if(isfield(p,'useGPU'))
    useGPU = p.useGPU; 
else
    useGPU = 0;
end

aliasLimit = 3;

sizeV = size( vels );
% [txangles,rxangles] = meshgrid(p.angles,p.rxangles);
% txangles = txangles.';
% rxangles = rxangles.';
% txangles = txangles(:);
% rxangles = rxangles(:);
rxangles = p.rxVals';
txangles = p.txVals';


%   vectorEst = single(zeros(size(R1,1)*size(R1,2),2));
vMat = (reshape( vels,[size(vels,1)*size(vels,2),size(vels,3)]) ).';
%aMat = [-sin(p.angles) (1+cos(p.angles))]./2; % When rx steering angle is always zero 
aMat = [-sin(txangles)-sin(rxangles) cos(txangles)+cos(rxangles)]/2; % Valid also when using multiple steering angles on tx and rx

pseudoInv = pinv(aMat);

nrAngles = sizeV(3);

vRes = ( aMat*pseudoInv-eye(nrAngles) )*vMat;       
sumvRes = sum( vRes.^2 );
residuals = min( sumvRes);
vectorEst = ( pseudoInv*vMat ).';

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