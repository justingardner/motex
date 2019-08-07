% motexGetEdges.m
%
%      usage: motexGetEdges(timeseries,cutoff)
%         by: justin gardner
%       date: 12/08/03
%       e.g.: motexGetEdges(timeseries,4.5)
%    purpose: returns edge fall and rise times
%
function edges = motexGetEdges(timeseries,cutoff,dilation)

% check command line arguments
if (nargin == 2)
  dilation = 1;
elseif (nargin ~=3)
  help motexGetEdges;
  return
end

timeseries = timeseries(:)';

% find when timecourse exceeds cutoff value (if cutoff is negative
% assume that means less than cutoff. If cutoff is positive assume
% we are looking for values larger than cutoff).
if (cutoff < 0)
  cutofftimes = timeseries < cutoff;
else
  cutofftimes = timeseries > cutoff;
end

% no events, give up
if (isempty(find(cutofftimes)))
  edges.cutofftimes = cutofftimes;
  edges.rising = [];
  edges.falling = [];
  edges.n = 0;
  edges.dilated = edges.cutofftimes;
  return
end

% find rising edges
rising = [0 find(diff(find(cutofftimes)) > 1)]+1;
% make sure that last one is not problematic
if (rising(length(rising)) > length(cutofftimes))
  rising(length(rising)) = risingedgretimes(length(rising))-1;
end

% find falling edges
falling = [find(diff(find(cutofftimes)) > 1) length(find(cutofftimes))];

% match length
if (length(rising) ~= length(falling))
  falling = falling(1:length(rising));
end

% get times where the signal is over cutoff
findcutofftimes = find(cutofftimes);

% pack return structure
edges.cutofftimes = cutofftimes;
edges.rising = findcutofftimes(rising);
edges.falling = findcutofftimes(falling);

% dilate edges 
dilatedrise = edges.rising-dilation;
dilatedrise(dilatedrise <= 0) = 1;
dilatedfall = edges.falling+dilation;
dilatedfall(dilatedfall > length(timeseries)) = length(timeseries);

% set dilated edges to true
edges.dilated = cutofftimes;
for i = 1:length(dilatedrise);
  edges.dilated(dilatedrise(i):dilatedfall(i)) = 1;
end

edges.n = length(edges.rising);

% get edge length
edges.completeN = min(length(edges.rising),length(edges.falling));
edges.len = edges.falling(1:edges.completeN) - edges.rising(1:edges.completeN);