function fts = TensorFreq( tensor, ff, dur )
% TensorFreq computes the FT at one frequency per stimulus
% 
% fts = TensorFreq( tensor, ff, dur ) computes the FT of the data in
% tensor, which has dimensions tensor{nstim}(nx,ny,nt). The frequency is
% indicated in ff, which is either a number or a vector (one frequency per
% stimulus).
%
% 2004-05 MC + AB
% 2005-12 MC removed option to not indicate "dur" (was buggy)
% 2005-12 MC changed scaling for case f=0 (was twice the correct one)
% 2005-12 MC changed the operation to 3D to save time
% 2006-01 MC made sure the fts are double (not single)
% 2008-10 LB added check for tensor dimensions and made matlab happy by
%   replacing '[1:nt]' by (1:nt), and | by ||, initializing fts{}

if nargin < 3
    error('You need to specify the duration dur');
end

% check the dimensions of the tensor
nstim = length(tensor);
for istim = 1 : nstim
    if ndims(tensor{istim}) ~= 3
        error('<TensorFreq> tensor must have format tensor{nstim}(nx,ny,nt)');
    end
end

% if the ff are all the same, make them become one number
if length(unique(ff))==1, ff = ff(1); end

[nx, ny, nt] = size( tensor{1} );

tt = (1:nt)'/nt *dur;

yy = zeros(1,1,nt);
aaa = zeros(nx,ny,nt);
fts = cell(1,nstim);
for istim = 1:nstim
    if istim==1 || length(ff)>1
        if ff(istim)==0
            yy(:) = ones(size(tt));
        else
            yy(:) = 2*exp(- tt*2*pi*1i*ff(istim));
        end
        aaa = repmat(yy,[nx,ny,1]);
    end
    fts{istim} = double(mean(tensor{istim}.*aaa, 3));
end


return

