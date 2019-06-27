function out = loadWFtrials(fpath,trialfiles, dat, varargin)
% generalizized, multi-purpose function for loading BVS/behavioral
% widefield imaging data. functionality determined by the inputParser
% object (see below), which sets available parameter-value pairs for the
% function. 
%   author's note: in principle, this should work for retinotopy data also,
%   i just need to append the appropriate variables
%   secondary note: if you edit this file, check up on "loadWFstruct_fast",
%   to make sure that everything is consistent.

% -basic code from mohammad (mainNeural, mainWidefield, ...) 2015
% -mjm, last update: 20160303
% -mohammad, major update 20170110: changed the order of processing steps, added plotting, removed redundant parts; 

if isempty(dat), fpsPCO =5; else if isfield(dat,'fps'), fpsPCO = dat.fps(1); else fpsPCO=30; end; end   
if fpsPCO ==50, HiCutFq = 10; elseif fpsPCO ==30, HiCutFq = 10; elseif fpsPCO ==5, HiCutFq = 5;else HiCutFq = fpsPCO/2; end

%% 001 // input parsing
    ip = inputParser;
    % which trials of the set will you load?
    addParameter(ip,'tset',1:numel(trialfiles),@isnumeric);
    addParameter(ip,'trimendframe',2,@isnumeric);
    % will you perform motion correction?
    addParameter(ip,'moco_flag',true,@islogical);
    addParameter(ip,'mimg',[],@isnumeric);
    % will you fit/apply a projective transform to these data?
    addParameter(ip,'sessregmethod','LK'); % manual (CP) or auto (LK)
    addParameter(ip,'sessregim',[],@isnumeric); % reference image (fit)
    addParameter(ip,'sessregT',[]);  % transform (apply)
    % dF/F // spatial and temporal filters
    addParameter(ip,'filter_flag',[],@islogical);
    addParameter(ip,'F0',[],@isnumeric);         % dF/F normalization
    addParameter(ip,'LoCutFreq',0.1,@isnumeric);%0.1
    addParameter(ip,'HiCutFreq', HiCutFq,@isnumeric); % temporal band-pass freqs
    addParameter(ip,'spatwid',3,@isnumeric);     % spatial gaussian width
    addParameter(ip,'PreFrm',[],@isnumeric);   % frame indices for frame0
    addParameter(ip,'fpsPCO',fpsPCO,@isnumeric);     % camera framerate
    % will you resize the image? (keeping pixels is too costly)
    addParameter(ip,'downsamp',[],@isnumeric); % give the major dim size
    % ROI mask type ('none','auto','manual')
    addParameter(ip,'masktype','none',@ischar);
    % timing data (use special frame zero if nonempty)
    addParameter(ip,'eventdata',[]);
    % tiling analysis (in place of actual visual segmentation)
    addParameter(ip,'tilesz',[],@isnumeric);
    % BEHAVIORAL DATA SPECS
        % will you use retinotopy to extracted segmented time sequences?
        addParameter(ip,'visareas',[],@isnumeric);
        addParameter(ip,'areatavg_flag',false,@islogical);
        % if so, will you average over pixels (0) or keep pixel info (1)?
        addParameter(ip,'keeppixels',false,@islogical);
    % RETINOTOPY DATA SPECS
        % will you trial-average? will you keep all trials?
        addParameter(ip,'trialavg_flag',false,@islogical);
        addParameter(ip,'alltrials_flag',false,@islogical);
        addParameter(ip,'plot_flag',[],@islogical);
    parse(ip,varargin{:}); struct2vars(ip.Results);
    
%% 002 // loading trials (main loop)
out = struct('xavg',[],'xall',[],'x',[],'xpix',[],'pixref',[],'datamask',[],'visareas',[],'bg',[]);
countr = 0;
for t =tset
    %% load the data and rotate it
    countr = countr + 1;
    fprintf('trial %d of %d: %s\n',countr,numel(tset), trialfiles{t});
    filename   = fullfile(fpath, trialfiles{t});
    if exist(filename,'file')
        [~, ~, ~, trialdata, ~] = loadPCOFile(filename); %[nRows, nCols, timeStamps, imgData, startTime]= loadPCOFile(filename);
    else
        out = {'no file'};
        continue
    end
    out.bg = cat(1,out.bg,mean(trialdata(:)));
    trialdata  = (fliplr(imrotate(single(trialdata),-90))); % rotate PCO camera images
    % this (fliplr(imrotate(single(trialdata),-90))) the same a transpose done on every frame MA 20190410
    
    % resize sessregim (because binning changes the size of trialdata)
    if ~isempty(sessregim)
        sessregim = imresize(sessregim,[size(trialdata,1) size(trialdata,2)]);
    end
    
    %% check if we 
    fncrop = fullfile(fpath(1:(end-4)),'crop.mat');
    if exist(fncrop,'file')
        disp('CROPPING THIS SESSION ADDITIONALLY')
        load(fncrop)
        pos = round(pos);
        %trialdata = trialdata((pos(2)+1):(pos(2)+pos(4)),(pos(1)+1):(pos(1)+pos(3)),:);
        trialdata = trialdata((pos(2)):(pos(2)+pos(4)-1),(pos(1)):(pos(1)+pos(3)-1),:);
    end

    %% 003 // correct motion
    ds = [];
    if moco_flag
        niter = 10;
        if isempty(mimg) % first trial only (or sent-in variable)
            fprintf('iterating for mean image ...\n');
            [~, mimg, ~, ~, ~] = align_iterative_sa(trialdata, niter);
        end
        fprintf('correting for motion ...\n')
        [ds, ~]   = registration_offsets_sa(trialdata, mimg, 0, Inf);
        trialdata = register_movie_sa(trialdata, ds);
    end
    
    %% Apply ROI in case the size of trialdata and sessregim do not match (trialdat>sessregim)
    if ~isempty(dat)
        if exist(['\\Labserver5\data\MOUSE\Segmentation\animalseg\Roi_' dat.id(2:end) '.mat'],'file')
            [NR,NC,~]=size(trialdata); 
            [nr,nc] = size(sessregim);
            if nr < NR && nc<NC
                load('\\Labserver5\data\MOUSE\Segmentation\animalseg\Roi_15309.mat','Roi');
                trialdata = func_TensorRoi(trialdata,Roi);
            end
        end
        if plot_flag && t==1
            drawnow
            fgh = figure('color','w'); set(fgh,'position',[680,560,960,420]);
            %subplot(1,2,1);imshow(mean(sessregim,3),[]); title([dat.animal ', refrence img'],'interpreter','none')
            %subplot(1,3,2);imshow(mean(trialdata,3),[]); title('bvs img before reg')
            fitter = sessregim/prctile(sessregim(:),99);
            thisImg = mean(trialdata,3);
            fittee = thisImg/max(thisImg(:));%(X-u_X)./s_X;
            subplot(1,2,1); imshow(cat(3,fitter,fittee,zeros(size(sessregim))),[]); title('before reg')
        end
    end
    
    %% 004 // register with other session, if available
    if ~isempty(sessregT) % this will run from trials 2-end
        fprintf('registering with retinotopy ...\n');
        [NR,NC] = size(sessregim);
        newtrialdata = nan(NR,NC,size(trialdata,3));
        %newtrialdata = nan(size(trialdata));

        switch sessregmethod
            case 'LK'
                for k = 1:size(trialdata,3)
                    [newtrialdata(:,:,k), supportLK] = iat_inverse_warping...
                        (trialdata(:,:,k), sessregT, par.transform, 1:NC, 1:NR);  %#ok<ASGLU>
                end
            case 'CP'
                %R = imref2d(size(sessregim)); 
                R = imref2d([NR,NC]); 
                out.datamask = imwarp(ones(size(trialdata(:,:,1))),sessregT,'OutputView',R,'Interp','Nearest');                
                for k = 1:size(trialdata,3)
                    newtrialdata(:,:,k) = imwarp(trialdata(:,:,k),sessregT,'OutputView',R,'Interp','Nearest');
                end
        end
        trialdata = newtrialdata; clear newtrialdata
    else
        cprintf('m','NO registration found. run "cpalignsession". this will continue without registration!\n');
    end
    
    %% 005 // ROI mask
    % df/f, pt 1
    if isempty(F0) && filter_flag
        sz = size(trialdata);
        F0mask = func_inscrsquareROI(trialdata,false);
        F0 = mean(reshape(trialdata(F0mask(:,round(sz(2)/2)),F0mask(round(sz(1)/2),:),:),[],1));
    end
    switch masktype
        case 'auto'
            % use the square inscribed in the circle inscribed in the rect. ROI
            mask = func_inscrsquareROI(trialdata,false);
            
            if strcmp(dat.id,'A15309') % MA added for mouse 15309 offcentered retino
                npix = 60; % num pix for translating the mask to capture later visual areas
                mask(:,end-npix+1:end) = [];
                mask = [false(size(mask,1),npix), mask];    
            end
            
            mask = imdilate(mask,strel('square',60)); % hack 160409%mjm
            
            sz = size(trialdata);
            if ~isempty(visareas) && countr <= 1 % only do it once
                visareas = bwlabel(imresize(visareas>eps,size(trialdata(:,:,1))),4);
                visareas = visareas(mask(:,round(sz(2)/2)),mask(round(sz(1)/2),:));
            end
            trialdata = trialdata(mask(:,round(sz(2)/2)),mask(round(sz(1)/2),:),:);
        case 'edge'
            % use the square inscribed in the circle inscribed in the rect.
            % ROI, including also the top block of the boundary
            mask = func_inscrsquareROI(trialdata,true);
            mask = imdilate(mask,strel('square',60)); % hack 160409%mjm
            sz = size(trialdata);
            if ~isempty(visareas) && countr <= 1 % only do it once
                visareas = bwlabel(imresize(visareas>eps,size(trialdata(:,:,1))),4);
                visareas = visareas(mask(:,round(sz(2)/2)),mask(round(sz(1)/2),:));
            end
            trialdata = trialdata(mask(:,round(sz(2)/2)),mask(round(sz(1)/2),:),:);
        case 'none'
            mask = [];
        case 'manual'
            mask = []; % not implemented/supported yet... (why would you want this?)
    end
%     fprintf('4 - data roi **BEFORE DF/F AND FILTERS\n'); keyboard;
    
    if ~isempty(dat)
        if plot_flag && t==1
            if isempty(mask)
                mask = func_inscrsquareROI(trialdata,false);
            end
            tmp1 = nanmean(trialdata,3); tmp2 = tmp1(mask(:,round(sz(2)/2)),mask(round(sz(1)/2),:));
            tmp1 = (tmp1 - nanmean(tmp2(:)))./nanstd(tmp2(:)); 
            subplot(1,2,2); imshow(cat(3,fitter,tmp1,zeros(size(sessregim))),[]); title('after reg')            
            clear tmp1 tmp2
            mask = [];
            drawnow
            pause(1)
        end
    end
    
    %% 006 // deltaF/F & spatial-temporal filtering       
    if filter_flag
        fprintf('deltaF/F normalizing...\n'); % df/f normalization should be done before tf and sf
        trialdata = (trialdata-F0)./F0;
        
        fprintf('spatial filtering [%.2f]...\n',spatwid);
        tmp = func_SpatFiltering(spatwid, trialdata); % Spatial Filtering
        trialdata = tmp{1};
        
        fprintf('temporal filtering [%.2f %.2f]...\n',LoCutFreq,HiCutFreq);
        tmp = func_TempFiltering( LoCutFreq, HiCutFreq, trialdata, fpsPCO ); % Temporal Filtering
        trialdata = tmp{1};
        
        if isempty(eventdata)
            fprintf('frame-0 correcting...\n');
            frame0 = nanmean(trialdata(:,:,PreFrm),3);
            if all(isnan(frame0(:))), frame0 = 0; end;
        else % 3 frames immediately before stimulus onset (edit 160419@mjm per ab suggestion)
            fprintf('frame-0 correcting pre-stimonset...\n');
            endsamp = floor((eventdata(t).stimulusCueStartedTime - eventdata(t).stimulusBackgroundStartedTime)*fpsPCO);
            frame0 = nanmean(trialdata(:,:,endsamp-4:endsamp),3);
        end
        trialdata = bsxfun(@minus,trialdata,frame0);
    end
    if ~isempty(downsamp) && abs(downsamp-1) > eps
            fprintf('downsampling... ');
        newtrialdata = [];
        for n = 1:size(trialdata,3)
            newtrialdata = cat(3,newtrialdata,imresize(trialdata(:,:,n),downsamp));
        end
        trialdata = newtrialdata;
    end
    
    %% 007 // apply visual area signal extraction, if available
    if ~isempty(visareas) && isempty(tilesz) % regular area selection
        visareas = bwlabel(imopen(imresize(visareas>eps,size(trialdata(:,:,1))),strel('disk',1)),4);
        xpix = cell(1,max(visareas(:)));
        pixref = xpix;
        for iSeg = 1:numel(xpix)
            fprintf('calculating timeseries for visual area %d...\n',iSeg);
            Ind = (visareas==iSeg);
            for iFrm = 1:size(trialdata,3)
                dataframe = trialdata(:,:,iFrm);
                xpix{iSeg} = cat(3,xpix{iSeg},dataframe(Ind));
            end
            pixref{iSeg} = int16(Ind.*reshape(cumsum(Ind(:)),size(Ind)));
        end
        x = cell2mat(cellfun(@nanmean,xpix,repmat({1},size(xpix)),'uni',0));
        % store one of the above (it's redundant to store both x and xpix)
        if keeppixels
            out.xpix = cat(1,out.xpix,xpix);
            if isempty(out.pixref) % don't duplicate this for each trial
                out.pixref = pixref;
            end
        end
        out.x = cat(1,out.x,{x});
            % note: we keep this in cell form because the timings in the
            % behavioral data are not necessarily consistent, and so the
            % time points do not necessarily line up -mjm 160310
        out.visareas = visareas;
    elseif ~isempty(visareas) && ~isempty(tilesz) % tiling analysis
        % note: this seems to rely on the 'edge' masking filter being implemented
        fprintf('sampling tiles...\n');
        % make the filter and convolve with the data
        tilesz = tilesz + mod(tilesz,2); % let's make tilesz an even number
        ff = fspecial('gaussian',tilesz,tilesz/6);
        tdat = convn(trialdata,ff,'valid');
        sz = size(tdat);
        tdat = tdat(1:(tilesz/2):sz(1),1:(tilesz/2):sz(2),:);
        tdat = reshape(tdat,1,size(tdat,1)*size(tdat,2),[]);
        % make a legend showing what pixel corresponds to where
        innergrid = zeros(size(trialdata,1),size(trialdata,2));
        innergrid((tilesz/2):(tilesz/2):(size(innergrid,1)-(tilesz/2)),(tilesz/2):(tilesz/2):(size(innergrid,2)-(tilesz/2))) = 1;
        if ~isempty(mask)
            fullgrid = zeros(size(mask));
            fullgrid(mask) = bwlabel(innergrid);
        else
            fullgrid = bwlabel(innergrid);
        end
        out.pixref = fullgrid;
        out.visareas = visareas;
        out.x = cat(1,out.x,{tdat});
    elseif isempty(visareas) && isempty(tilesz) && ~trialavg_flag % stack raw files into output
        out.x = cat(1,out.x,{trialdata});
    end

    %% 008 // package into trial average/single trial collection
    fprintf('packaging...\n\n');
    % add current trial to the trial average
    if trialavg_flag
        % make every trial as large as the longest trial (CHANGE?)
        if isempty(out.xavg)
            out.xavg = trialdata./numel(tset);
        else
            if size(out.xavg,3) > size(trialdata,3)
                sz = size(trialdata);
                trialdata = cat(3,trialdata,nan(sz(1),sz(2),size(out.xavg,3)-sz(3)));    
            elseif size(out.xavg,3) < size(trialdata,3)
                sz = size(out.xavg);
                out.xavg = cat(3,out.xavg,nan(sz(1),sz(2),size(trialdata,3)-sz(3)));
            end 
            out.xavg = nansum(cat(4,out.xavg,trialdata./numel(tset)),4);
        end
    end
    % append current trial to collection of (unaveraged) other trials
    if alltrials_flag
        % note that we don't have the hard disk space to 1) save or 2)
        % store all of the trials at once -- instead, we're going to store
        % a structure to rebuild the preprocessed trials as quickly as
        % possible. send these structure array elements one-by-one to
        % 'loadWFstruct_fast' -mjm
        out.xall(t,1).fnm = filename;
        out.xall(t,1).trimendframe = trimendframe;
        out.xall(t,1).regshift = ds;
        out.xall(t,1).ROImask = mask;
        if ~isempty(sessregT) % ==> this trial was session-registered
            out.xall(t,1).sessregim   = sessregim;
            out.xall(t,1).sessregT    = sessregT;
        end
        if filter_flag        % ==> this trial was filtered
            out.xall(t,1).F0          = F0;
            out.xall(t,1).loCutFreq   = LoCutFreq;
            out.xall(t,1).hiCutFreq   = HiCutFreq;
            out.xall(t,1).spatwid     = spatwid;
            out.xall(t,1).PreFrm      = PreFrm;
        end

    end
end
