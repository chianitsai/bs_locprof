numObjectsValues = zeros(size(frames)); % Initialize the array
for i = 1:numel(frames) % Loop over all elements in frames
  if isfield(frames(i).cells, 'NumObjects') % Check if the structure has the field 'NumObjects'
    numObjectsValues(i) = frames(i).cells.NumObjects; % Access the NumObjects field using dot notation
  else
    warning('Structure %d does not contain a field "NumObjects".', i); % Issue a warning
    numObjectsValues(i) = NaN; % Placeholder for missing NumObjects
  end
end

totalNumObjects = sum(numObjectsValues); % Sum all values in the numObjectsValues array
disp('Total sum of NumObjects values:');
disp(totalNumObjects); % Display the total sum