test = zeros(1,size(fieldnames(final),1));

for i = size(fieldnames(final),1)
    currentanimal = char(animalIDs(i));
   
    if isfield(final.(currentanimal), session_to_analyze)
        test(1,i) = 1;
    elseif ~isfield(final.(currentanimal), session_to_analyze)
        test(1,i) = 2;
    end

end