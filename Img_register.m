fprintf('computing session transform...\n');
    load(sriDir);
    % make sure the size of sessregim matches trialdata
    sessregim = imresize(sessregim,[size(trialdata,1) size(trialdata,2)]);
    % initialize the transform...
    par.transform = 'homography';
    par.levels = 3;
    par.iterations = 30;
    % ... and then fit it
    X = nanmean(trialdata,3);
    % z-score the images according to the stats of their CENTERS
    % (the boundaries are useful for registration but often have
    % strange limits -- so we don't use them for this part) -mjm
    mask = func_inscrsquareROI(trialdata,false);
    sz = size(trialdata);
    centsri = sessregim(mask(:,round(sz(2)/2)),mask(round(sz(1)/2),:));
    u_sri = nanmean(centsri(:)); s_sri = nanstd(centsri(:));
    centX = X(mask(:,round(sz(2)/2)),mask(round(sz(1)/2),:));
    u_X = nanmean(centX(:)); s_X = nanstd(centX(:));
    
    %fitter = (sessregim-u_sri)./s_sri;
    fitter = sessregim/prctile(sessregim(:),99);
    fittee = X/max(X(:));%(X-u_X)./s_X;
    fprintf('registering session...\n');
    [moving_out,fixed_out] = cpselect(fittee,fitter,'wait',true);
    moving_out = cpcorr(moving_out,fixed_out,fittee,fitter);
    sessregT = fitgeotrans(moving_out,fixed_out,'projective');
    % then perform the registration on all frames
    [M,N] = size(sessregim);
    newtrialdata = nan(M,N,size(trialdata,3));
    R = imref2d([M N]);
    for m = 1:size(trialdata,3)
        newtrialdata(:,:,m) = imwarp(trialdata(:,:,m),sessregT,'OutputView',R);
    end
    trialdata = newtrialdata; clear newtrialdata
    
    %% 04 // validation and saving
    % plot to confirm effectiveness of the mapping
    fgH=figure(101);
    set(fgH,'color','w','position',[50,580,1080,370]);
    subplot(1,3,1);
    imagesc((sessregim-u_sri)./s_sri); axis image; set(gca,'clim',[-3 3])
    title('retino','fontsize',15)
    
    subplot(1,3,2);
    imagesc((X-u_X)./s_X); axis image; set(gca,'clim',[-3 3])
    title('this sess','fontsize',15)
    
    tmp1 = nanmean(trialdata,3); tmp2 = tmp1(mask(:,round(sz(2)/2)),mask(round(sz(1)/2),:));
    tmp1 = (tmp1 - nanmean(tmp2(:)))./nanstd(tmp2(:));
    
    subplot(1,3,3);
    imagesc(tmp1); axis image; set(gca,'clim',[-3 3])
    title('this sess reg to retino','fontsize',15)
    
    fgH=figure(102);
    set(fgH,'color','w','position',[1135,580,670,370]);
    subplot(1,2,1);
    imshow(cat(3,fitter,fittee,zeros(size(tmp1))));
    title('before reg','fontsize',15)
    subplot(1,2,2);
    imshow(cat(3,fitter,tmp1,zeros(size(tmp1))));
    title('after reg','fontsize',15)
    drawnow;