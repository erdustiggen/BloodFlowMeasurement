function [finalEsts, corrValtot, chosenAliasPat] = blockMatching( iqhpComp, vCand, ...
            aliasMat, aliasPat, bfPars, npx)

    lowResPacketSize = size( iqhpComp, 4);

    %% Find correlation for all candidate velocity vectors using compounded data (no GPU)   

        kernelsize = [npx npx];

        dx = bfPars.x_axis(2)-bfPars.x_axis(1);
        dz = bfPars.z_axis(2)-bfPars.z_axis(1);
        nrpoints = size( vCand,2);
        corrVal = zeros( size( vCand, 3), nrpoints, lowResPacketSize-1);
        gData = iqhpComp;
        gAxisX = bfPars.x_axis;
        gAxisZ = bfPars.z_axis;

        coords_origx = -(kernelsize(1)-1)/2*dx:dx:(kernelsize(1)-1)/2*dx;
        coords_origz = -(kernelsize(2)-1)/2*dz:dz:(kernelsize(2)-1)/2*dz;
        [oXb, oZb] = meshgrid(coords_origx, coords_origz);
        oX = repmat( oXb, [1 nrpoints] )+ reshape( repmat( bfPars.mymaskX(1:nrpoints).', [kernelsize(1)*kernelsize(2) 1] ), [kernelsize(1) kernelsize(2)*nrpoints] );
        oZ = repmat( oZb, [1 nrpoints] )+ reshape( repmat( bfPars.mymaskZ(1:nrpoints).', [kernelsize(1)*kernelsize(2) 1] ), [kernelsize(1) kernelsize(2)*nrpoints] );

        for ictr = 1:size( vCand, 3),

            sX = oX + reshape( repmat( vCand(1,1:nrpoints,ictr)*bfPars.c/(2*bfPars.f_demod)/bfPars.PRF, ...
                [kernelsize(1)*kernelsize(2) 1] ), [kernelsize(1) kernelsize(2)*nrpoints] );
            sZ = oZ - reshape( repmat( vCand(2,1:nrpoints,ictr)*bfPars.c/(2*bfPars.f_demod)/bfPars.PRF, ...
                [kernelsize(1)*kernelsize(2) 1] ), [kernelsize(1) kernelsize(2)*nrpoints] );

            for timestep = 1:lowResPacketSize-1,
                data_orig = interp2(gAxisX ,gAxisZ, abs(gData(:,:,1,timestep) ), oX, oZ);
                data_orig = reshape(data_orig, [kernelsize(1)*kernelsize(2) nrpoints]);
                data_orig = data_orig-repmat( mean(data_orig,1), [kernelsize(1)*kernelsize(2) 1] );
                orignorm = sqrt( sum( data_orig.^2, 1) );
                
                data_shift = interp2(gAxisX,gAxisZ, abs(gData(:,:,1,timestep+1) ), sX, sZ);
                data_shift = reshape(data_shift, [kernelsize(1)*kernelsize(2) nrpoints]);
                data_shift= data_shift-repmat( mean(data_shift,1), [kernelsize(1)*kernelsize(2) 1] );
                shiftnorm = sqrt( sum( data_shift.^2, 1) );
                
                corrVal(ictr, :, timestep) = sum( conj(data_orig).*data_shift, 1 )./(orignorm.*shiftnorm);
            end

        end
        
        corrValtot = mean( corrVal, 3);


        spanVals = max( aliasMat,[],1);
        spanVec = spanVals( aliasPat);
        for i = 1:max( spanVec),
            corrValtot(end-i+1:end, spanVec >= i) = 0;
        end


        [amps, chosenAliasPat] = max( corrValtot, [], 1);

%         chosenAliasPat(aliasPat == 1) = find( aliasfactTab == 0);

        inds = [(chosenAliasPat-1).*ones(1,nrpoints)*2*size(vCand,2)+(0:nrpoints-1)*2+1; ...
            (chosenAliasPat-1).*ones(1,nrpoints)*2*size(vCand,2)+(0:nrpoints-1)*2+2];

        finalEsts = vCand(inds)*bfPars.c/(2*bfPars.f_demod);

        %disp(sprintf('Number of origopoints for tracking is %i, kernelsize is %i x %i,averaged timesteps is %i',nrpoints,kernelsize(1),kernelsize(2),size(corrVal,3)));


%         corrValtot = gather(corrValtot);
%         chosenAliasPat = gather(chosenAliasPat);
        
   