% faster implementation of vectorDoppler using precalculated aliasing
% residuals

function [vectorEst, aliasPat, aliasMat, minRes] = vectorDopplerWeighted(R1,p,S)



aliasLimit = 3;

sizeR1 = size( R1 );
% [txangles,rxangles] = meshgrid(p.angles,p.rxangles);
% txangles = txangles.';
% rxangles = rxangles.';
% txangles = txangles(:);
% rxangles = rxangles(:);
rxangles = p.rxVals';
txangles = p.txVals';

if isreal(R1)
    
        vectorEst = single(zeros(size(R1,1), size(R1,2),2));
        temp = vectorEst;
        temp = reshape(temp, [size(temp,1)*size(temp,2),size(temp,3)]);
        fMat = (reshape(R1,[size(R1,1)*size(R1,2),size(R1,3)])*p.PRF/(2*pi) ).';

        aMat = [-sin(p.angles) (1+cos(p.angles))]./2;


        pseudoInv = pinv(aMat);
        for kk = 1:size(fMat,2)
           temp(kk,:) = pseudoInv*fMat(:,kk);
        end


        vectorEst = reshape(temp,[size(vectorEst,1), size(vectorEst,2), size(vectorEst,3)]);

else

    if sizeR1(2) ~= 1
%         vectorEst = single(zeros(size(R1,1)*size(R1,2),2));
        fMat = (reshape(angle(R1),[size(R1,1)*size(R1,2),size(R1,3)])*p.PRF/(2*pi) ).';
        fMat = fMat.*S;
        %aMat = [-sin(p.angles) (1+cos(p.angles))]./2; % When rx steering angle is always zero 
        aMat = [-sin(txangles)-sin(rxangles) cos(txangles)+cos(rxangles)]./2; % Valid also when using multiple steering angles on tx and rx

%         pseudoInv = pinv(aMat);
        pseudoInv = zeros([size( aMat,2) size( aMat,1) size(S,2)] );
        for i = 1:size(S,2)
            pseudoInv(:,:,i) = pinv( repmat( S(:,i),[1 2]).*aMat );
        end
        rPinv = reshape( pseudoInv, [size( pseudoInv,1) size( pseudoInv,2)*size( pseudoInv,3)] );

        nrAngles = sizeR1(3);
        nrCases = aliasLimit^nrAngles; %including no aliasing, but not aliasing on all angles
        sprintf('Number of aliasing cases/patterns is %i',nrCases);
        aliasMat = zeros(nrAngles, nrCases );
        for patternCode = 0:nrCases-1,
            str = char( dec2base( patternCode, aliasLimit) );
            tail = str2num( str(:) );
            pattern = zeros(nrAngles, 1);
            pattern(end-length(tail)+1:end) = tail;
            aliasMat(:, patternCode+1) = pattern;
        end
        
        validinds = find( any( ~aliasMat, 1) );
        aliasMat = aliasMat(:, validinds);
        nrCases = size( aliasMat,2);
        nrPoints = size(fMat,2);

        temp = reshape(aMat*rPinv, [size( pseudoInv,2) size( pseudoInv,2) size( pseudoInv,3)] ) ...
            - repmat( eye( nrAngles), [1 1 size( pseudoInv,3)] );
        temp2 = reshape( permute( temp, [1 3 2] ), [size( temp,1)*size(temp,3) size( temp,2)] );
%         temp2 = temp2*aliasMat*p.PRF;
%         temp2 = permute( reshape( temp2, [size( temp,1) size(temp,3) size( aliasMat,2)] ), [1 3 2] );
        gRes = permute( reshape(temp2*aliasMat*p.PRF, [size( temp,1) size(temp,3) size( aliasMat,2)] ), [1 3 2] );
%         fRes = temp2*fMat;
        fRes = sum( temp.*repmat( permute( fMat, [3 1 2] ), [7 1 1]), 2);
        temp = zeros(nrCases, size(R1,1)*size(R1,2));
        for i = 1:nrAngles,
            temp = temp + squeeze( ( repmat( fRes(i,:,:), [1 size(gRes,2) 1])+gRes(i,:,:) ).^2 );
        end
%         
%         gRes = ( ( aMat*pseudoInv-eye(nrAngles) )*aliasMat*p.PRF ).';
%         fRes = ( aMat*pseudoInv-eye(nrAngles) )*fMat;
%         
%         temp = zeros(nrCases, size(R1,1)*size(R1,2));
%         
%         for i = 1:nrAngles,
%             temp = temp + ( repmat( fRes(i,:), [size(gRes,1) 1])+repmat( gRes(:,i), [1 size(fRes,2)]) ).^2;
%         end
        
        
        [minRes, aliasPat] = min( temp );
        vectorEst = squeeze( sum( repmat( permute( fMat+aliasMat(:,aliasPat)*p.PRF, ...
            [3 1 2] ), [2 1 1] ).*pseudoInv, 2) );
%         vectorEst = ( pseudoInv*( fMat+aliasMat(:,aliasPat)*p.PRF ) ).';

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
        
    else
%         vectorEst = single(zeros(size(R1,1),2));

        % three types of aliasing
        % 1: largest component on leftmost angle
        % 2: largest component on rightmost angle
        % 3: components are first increasing, then decreasing, from left to
        % right

        fMat = ( angle(R1)*p.PRF/(2*pi) ).';

        aMat = [-sin(p.angles).'; (1+cos(p.angles)).']./2;
    %     aMat = [-sin((p.angles/2).'); cos((p.angles/2).')];
        pseudoInv = pinv(aMat).';
        vectorEst = pseudoInv*fMat';
        
        vectorEst = reshape( vectorEst, [size( R1,1), size(R1, 2), 2] );


    end
end

end