function [data, trials, varargin] = TrialFilter_test(data,varargin)
      VALID_PARS = {'BLOCK','TYPE','SHK','REW','OMIT','OMITALL','ALL','WSLS','STTOCHO','WINSTAY','LOSESHIFT','LOSEOMIT', 'LOSESTAY', 'WIN', 'LOSS', 'AA', 'BLANK_TOUCH'};
    
    % Initialize logical index arrays for each condition
    idx = true(height(data), 1);
    condition_idxs = cell(1, length(varargin) / 2); % Cell array to store condition indices
    condition_num = 1; 
    % Parse varargin
    for ii = 1:2:length(varargin)
        if ~ismember(upper(varargin{ii}), VALID_PARS)
            error('%s is not a valid parameter', upper(varargin{ii}));
        end
        
        % Initialize logical index array for current condition
        condition_idx = true(height(data), 1);
        
        % Filter data based on the current condition
        switch upper(varargin{ii})
            case 'BLOCK'
                condition_idx = ismember(data.Block, varargin{ii+1});
            case 'TYPE'
                condition_idx = (data.ForceFree == varargin{ii+1});
            case 'SHK'
                condition_idx = (data.shock == varargin{ii+1});
            case 'REW'
                condition_idx = (data.bigSmall == varargin{ii+1});
            case 'OMIT'
                condition_idx = (data.omission == varargin{ii+1});
            case 'OMITALL'
                condition_idx = (data.omissionALL == varargin{ii+1});
            case 'ALL'
                % No filtering needed for 'ALL'
            case 'WSLS'
                condition_idx = (data.WSLScode == varargin{ii+1});
            case 'WINSTAY'
                condition_idx = (data.win_stay == varargin{ii+1});
            case 'LOSESHIFT'
                condition_idx = (data.lose_shift == varargin{ii+1});
            case 'LOSEOMIT'
                condition_idx = (data.lose_omit == varargin{ii+1});
            case 'LOSESTAY'
                condition_idx = (data.lose_stay == varargin{ii+1});
            case 'WIN'
                condition_idx = (data.WL == varargin{ii+1});
            case 'LOSS'
                condition_idx = (data.WL == varargin{ii+1});
            case 'AA'
                condition_idx = (data.type_binary == varargin{ii+1});
            case 'BLANK_TOUCH'
                condition_idx = (data.Blank_Touch == varargin{ii+1});
        end
        
        % Store the condition index
        condition_idxs{condition_num} = condition_idx;
        condition_num = condition_num+1; 
    end
    
    % Find the intersection of indices
   if ~isempty(condition_idxs)
        idx = intersect([condition_idxs{:}], 'rows');
    end
    
    % Apply filtering
    data = data(idx, :);
    trials = table2cell(data(:, 1));
end