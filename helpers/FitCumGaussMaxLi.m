function [FitPa,fval,exitflag,fitx,fity] = FitCumGaussMaxLi(x,y, n_resp,...
    varargin)
 % hn 25/08/03
 % fits cumulative gaussian to data using maximum likelihood
 % x: values (e.g.) signal strength
 % y: number (not percentage) of 'detect' (preferred) responses
 % n_resp: number of samples contributing to each x,y datapoint
 % Options:
 %  'fixbase': base
 %  'fixscale': makes it go from yoffset to (yoffset+amp*1), e.g. for
 %  psychophysics: teh Cum Gaussian should go from 0.5 to 1 (i.e. yoff=0.5,
 %  scale=0.5
 
 % returns :
 % fitparams(1) = Mean
 %  "       (2) = SD
 %  "       (3) = y-offset (base)
 %  "       (4) = amplitude
 % fval      residual (sq) at fitparams
 % exitflag
 % fitx     (optional) dadapoints to plot fit
 % fity
basefix = 'off';
ampfix = 'off';
meanfix = 'off';
ploton = 0    ;   
MEAN = 1;
SD = 2;
BASE = 3;
AMP = 4;

guess (1:4) =NaN;
nvar = nargin - 3;
j = 1;
while j  <= nvar
  str = varargin{j};
  if strcmpi('log',str)
    type = varargin{j};    
  elseif strcmpi('lin',str)
    type = varargin{j};
  elseif strcmpi('plot',str)
    ploton = 1;
  elseif strcmpi ('amp', str)
    j = j+1;
    guess (AMP) = varargin{j};    
    elseif strcmpi ('sd', str)
        j = j+1;
        guess (SD) = varargin{j};
    elseif strcmpi ('mean', str)
        j = j+1;
        guess (MEAN) = varargin{j};
        

    elseif strcmpi ('meanfix', str)
        j = j+1;
        fixedmean = varargin{j};
        meanfix = 'on';
    elseif strcmpi ('base',str)
        j = j+1;
        guess (BASE) = varargin {j};
    elseif strcmpi ('basefix',str)
        j = j+1;
        basefix = 'on';
        fixedbase = varargin{j};
    elseif strcmpi ('ampfix',str)
        j = j+1;
        ampfix = 'on';
        fixedamp = varargin{j};
    
  end
  j = j+1;
end
if isempty (x) | size(x) ~= size(y) 
    
    disp( ' to break')
    
    FitPa.amp = [];
    FitPa.sd = [];
    FitPa.mean = [];
    FitPa.base = [];
    
    exitflag = [];
    fval = [];
    type = [];

    fitx = [];
    fity = [];
    
    return;
end


if isnan (guess(AMP)) 
    guess(AMP) = (max(y) - min(y));
end

if  isnan (guess(MEAN))
    [ maximum maxpos] = max(y);
   
    guess(MEAN) = x(maxpos);         
end

if isnan (guess(SD))
    guess(SD) = std(x .* y)/mean(y);
end

if isnan (guess(BASE))
    guess(BASE) = 0.5-min(y);
end

options = optimset('MaxFunEvals',1000,'maxiter',1000);
LB = [ 0 0 0 0];
exitflag1 = 0;
fval1 = NaN;
if strcmpi(basefix,'on') & strcmpi(ampfix,'on') & strcmpi(meanfix,'on')
     LB(3) = fixedbase;
     LB(4) = fixedamp;
     LB(1) = fixedmean;
     guess(AMP) = [];
    guess(BASE) = [];
    
        [fitparams,fval,exitflag,output]= fminsearch(@cumgaussllike_meanfix,...
            guess,options,x,y,n_resp,LB); %,lower_bound,upper_bound,options)
    fitparams(BASE) = fixedbase;
    fitparams(AMP) = fixedamp;
    fitparams(MEAN) = fixedmean;


elseif strcmpi(basefix,'on') & strcmpi(ampfix,'on')
     LB(3) = fixedbase;
     LB(4) = fixedamp;
     guess(AMP) = [];
    guess(BASE) = [];
    
        [fitparams,fval,exitflag,output]= fminsearch(@cumgaussllike,...
            guess,options,x,y,n_resp,LB); %,lower_bound,upper_bound,options)
    fitparams(BASE) = fixedbase;
    fitparams(AMP) = fixedamp;
    
else 
    fitparams(BASE) =0.5;
    guess(BASE) = 0;
    guess(AMP) = 0.5;
    guess(SD) = 0.2;
    guess(MEAN) = 0;
    
        [fitparams,fval,exitflag,output]= fminsearch(@cumgaussllike2,...
            guess,options,x,y,n_resp,LB); %,lower_bound,upper_bound,options)

end

FitPa.amp = fitparams(AMP);
FitPa.sd = fitparams(SD);
FitPa.mean = fitparams(MEAN);
FitPa.base = fitparams(BASE);

fitx = [];
fity = [];
incr = (max(x)-min(x))/100;
if ploton
fitx = min(x): incr : max(x);
fity = fitparams(AMP)*normcdf(fitx,fitparams(MEAN),fitparams(SD)) + fitparams(BASE);
end
% -------------------------------------------------------------------------
% -------
function result = cumgaussllike(X0,x,y,n_resp,LB);

f=[];
% y = normcdf
p = LB(4)* normcdf(x,X0(1),X0(2)) + LB(3); 
q= 1-p;

result = sum(((n_resp - y) .* log(q)) + ((y .* log(p))));
result = -result;
% -------------------------------------------------------------------------
% -------
% -------
function result = cumgaussllike_meanfix(X0,x,y,n_resp,LB);
f=[];
% y = normcdf
p = LB(4)* normcdf(x,LB(1),X0(2)) + LB(3); 
q= 1-p;

result = sum(((n_resp - y) .* log(q)) + ((y .* log(p))));
result = -result;
% -------------------------------------------------------------------------

function result = cumgaussllike2(X0,x,y,n_resp,LB);

f=[];
% y = normcdf
p = X0(4)* normcdf(x,X0(1),X0(2)) + X0(3); 
q= 1-p;

result = sum(((n_resp - y) .* log(q)) + ((y .* log(p))));
result = -result;
