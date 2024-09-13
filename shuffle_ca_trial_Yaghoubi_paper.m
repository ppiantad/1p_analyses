ca_test = ca(1, :);
[~, num_samples] = size(ca);
shuffled_data = zeros(uv.resamples, num_samples); % Preallocate matrix for efficiency
for hh = 1:uv.resamples


    shift_val = randi(num_samples) % Generate a random shift value for each signal RUAIRI RECOMMENDED KEEPING THE SAME SHIFT VAL, rather than randomizing per neuron. this is because then you keep the overall correlation b/w the neurons, but disrupt the relationship to the event timestamps

    shuffled_data(hh,:) = circshift(ca_test, shift_val,2); % Perform the circular shuffle




end
