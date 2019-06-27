function tensor = func_SpatFiltering(gSigma, tensor)
%% Temporal Filtering in Hz (in Freq domain)
global VDAQ

if nargin < 2
    % build filter
    [nr, nc, nt] = size(VDAQ.tensor{1});
    fFG = buildGfilt(gSigma,nr,nc);
    % perform filtering
    SpatFlt = cell(VDAQ.tensor);
    for istim   = 1: VDAQ.nstim;
        for frm = 1:nt
            %pad the original image
            if nr > nc
                padImg = repmat(VDAQ.tensor{istim}(:,1,frm),[1 abs(nr-nc)]);
                Img = [padImg VDAQ.tensor{istim}(:,:,frm)];
                fImg = fft2(Img);
                ImgFlt = real(ifftshift(ifft2(fImg .* fFG)));
                ImgFlt = ImgFlt(:,abs(nr-nc)+1:end);
            elseif nr < nc
                padImg = repmat(VDAQ.tensor{istim}(1,:,frm),[abs(nr-nc) 1]);
                Img = [padImg; VDAQ.tensor{istim}(:,:,frm)];
                fImg = fft2(Img);
                ImgFlt = real(ifftshift(ifft2(fImg .* fFG)));
                ImgFlt = ImgFlt(abs(nr-nc)+1:end,:);
            else
                Img = VDAQ.tensor{istim}(:,:,frm);
                fImg = fft2(Img);
                ImgFlt = real(ifftshift(ifft2(fImg .* fFG)));
            end
            SpatFlt{istim}(:,:,frm) = ImgFlt;
        end
    end
    VDAQ.tensor = SpatFlt;
else
    if ~iscell(tensor), tensor = {tensor}; end;
    % build filter
    [nr, nc, nt] = size(tensor{1});
    fFG = buildGfilt(gSigma,nr,nc);
    % perform filtering
    SpatFlt= cell(tensor);
    for istim   = 1: numel(tensor);
        for frm = 1:nt
            %pad the original image
            if nr > nc
                padImg = repmat(tensor{istim}(:,1,frm),[1 abs(nr-nc)]);
                Img = [padImg tensor{istim}(:,:,frm)];
                fImg = fft2(Img);
                ImgFlt = real(ifftshift(ifft2(fImg .* fFG)));
                ImgFlt = ImgFlt(:,abs(nr-nc)+1:end);
            elseif nr < nc
                padImg = repmat(tensor{istim}(1,:,frm),[abs(nr-nc) 1]);
                Img = [padImg; tensor{istim}(:,:,frm)];
                fImg = fft2(Img);
                ImgFlt = real(ifftshift(ifft2(fImg .* fFG)));
                ImgFlt = ImgFlt(abs(nr-nc)+1:end,:);
            else
                Img = tensor{istim}(:,:,frm);
                fImg = fft2(Img);
                ImgFlt = real(ifftshift(ifft2(fImg .* fFG)));
            end
            SpatFlt{istim}(:,:,frm) = ImgFlt;
        end
    end
    tensor = SpatFlt;
end

%% auxiliary functions
function filtout = buildGfilt(sig,r,c)
% Create the Filter
disp(['Gaussian Sigma: ' num2str(sig) ', Spatial filtering ...']);
if r >= c
    GF = fspecial('gaussian',[r r], sig); % when the filter and image are squre, the fft works better
else
    GF = fspecial('gaussian',[c c], sig);
end
GF = GF/sum(GF(:));
filtout = fft2(GF);





%% Log
%{
- Mohammad 2015; wrote the file
%}
