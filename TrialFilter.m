function [data, trials, varargin] = TrialFilter(data,varargin);
VALID_PARS = {'BLOCK','TYPE','SHK','REW','OMIT','OMITALL','ALL','WSLS','STTOCHO','WINSTAY','LOSESHIFT','LOSEOMIT', 'WIN', 'LOSS', 'AA'};
%for varargin
%JSFilter will assign the valid parameter equal to the argument that
%immediately succeeds that paramter.  So, for example, if you call the
%function as follows: d=JSFilter(data,BigRew_choice,'SESSION',SESSION3), it
%will set SESSION=Session3.  Parameters should be entered as strings, and
%assignment values followng the parameters should be vectors.


%initialize variables
BLOCK=[]; %BLOCK number (1,2,3,4,5)
TYPE=[]; %use 1 for force, 0 for free
SHK=[]; %use 0 for no shock, 1 for shock
REW=[]; % use 1.2 for big, 0.3 for small
OMIT=[]; %use 0 for non-omissions, 1 for free trial omissions, 2 for force trial omissions
ALL=[]; %if you want to see all trials, set ALL to 1 (or any number should work)
WSLS=[]; %1=win; 2=win+1 (trial after win); 3=loss; 4=loss+1;
WINSTAY=[]; % 1 = winstay
LOSESHIFT=[]; % 1 = lose/shift
LOSEOMIT=[]; % 1 = lose/omit
WIN=[]; % 1 = win
LOSS=[]; % 3 = loss
AA=[]; % 1 = large abort, 2 = small abort (changed 12/17/2021) (old {'large abort'}, {'small abort'})



% parse varargin
for ii = 1:2:length(varargin)
    if ~ismember(upper(varargin{ii}), VALID_PARS)
        error('%s is not a valid parameter', upper(varargin{ii}));
    end
    eval([upper(varargin{ii}) '=varargin{ii+1};']);  %sets the argument in equal to the value that follows it.  
end

% added on 1/22/2023. Some of these mice were run with an old version of
% AEBTtoTable that did not change feeder values of 0.5 to 0.3 (these are
% small reward sizes, but particular value is dependent on how well the
% feeder itself functions). This code will check if the filter variable is
% 0.3 (corresponding to small reward), and then will update the BehavData
% (data) to reflect that all 0.5s should really be 0.3s
for zz = 1:size(varargin, 2)
    if varargin{zz} == 0.3
        for pp = 1:size(data,1)
            if data.bigSmall(pp) == 0.5
                data.bigSmall(pp) = 0.3;
            end
        end
    else
end
end


for k=1:2:length(varargin)
    %Filtering by SESSION: SESSION will be set equal to the vector which
    %corresponds to the SESSION of interest.
    if strcmp(upper(varargin{k}),'BLOCK')
        trials = table2cell(data([find(data.Block==BLOCK)],1));
        data([find(data.Block~=BLOCK)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'TYPE')
        trials = table2cell(data([find(data.ForceFree==TYPE)],1));
        data([find(data.ForceFree~=TYPE)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'SHK')
        trials = table2cell(data([find(data.shock==SHK)],1));
        data([find(data.shock~=SHK)],:)=[];
    end
    
    
    if strcmp(upper(varargin{k}),'REW')
       trials = table2cell(data([find(data.bigSmall==REW)],1)); 
       data([find(data.bigSmall~=REW)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'OMIT')
        trials = table2cell(data([find(data.omission==OMIT)],1));
        data([find(data.omission~=OMIT)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'OMITALL')
        trials = table2cell(data([find(data.omissionALL==OMITALL)],1));
        data([find(data.omissionALL~=OMITALL)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'WSLS')
        trials = table2cell(data([find(data.WSLScode==WSLS)],1));
        data([find(data.WSLScode~=WSLS)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'ALL')
       trials = num2cell(data.Trial);
       disp('Using all trials')
    end
    
    if strcmp(upper(varargin{k}),'WINSTAY')
        trials = table2cell(data([find(data.win_stay==WINSTAY)],1));
        data([find(data.win_stay~=WINSTAY)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'LOSESHIFT')
        trials = table2cell(data([find(data.lose_shift==LOSESHIFT)],1));
        data([find(data.lose_shift~=LOSESHIFT)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'LOSEOMIT')
        trials = table2cell(data([find(data.lose_omit==LOSEOMIT)],1));
        data([find(data.lose_omit~=LOSEOMIT)],:)=[];
    end    
    
    if strcmp(upper(varargin{k}),'WIN')
        trials = table2cell(data([find(data.WL==WIN)],1));
        data([find(data.WL~=WIN)],:)=[];
    end    
    
    if strcmp(upper(varargin{k}),'LOSS')
        trials = table2cell(data([find(data.WL==LOSS)],1));
        data([find(data.WL~=LOSS)],:)=[];
    end   
    
    if strcmp(upper(varargin{k}),'AA')
        trials = table2cell(data([find(data.type_binary==AA)],1));
        data([find(data.type_binary~=AA)],:)=[];
        
    end   

end

end

