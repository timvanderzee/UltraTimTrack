function parms = adjust_hough_parameters(parms)

%% Init variables
adj_parms = [parms.fas.range parms.apo.super.maxangle parms.apo.deep.maxangle]; %get the 4 parms I can change
%% initialize text var for the inputs
% Create a figure window
fig = figure('Units', 'normalized', 'Position', [0.4 0.4 0.3 0.4], ...
    'Name', 'Set TimTrack settings', 'NumberTitle', 'off');

text_prompt = {'Min. fascicle angle (deg):', 'Max. fascicle angle (deg):', 'Max. superficial aponeurosis angle (deg):', 'Max. deep aponeurosis angle (deg):'};
fields = {'fas_min', 'fas_max', 'sup_apo', 'deep_apo'};
num_fields = numel(fields);
inputs = gobjects(num_fields, 1);

%% create the object to handle
% Create a figure window
for i = 1:num_fields
    uicontrol('Style', 'text', 'String', text_prompt{i}, ...
        'Units', 'normalized', 'Position', [0.1 0.8 - (i-1) * 0.18, 0.4, 0.1], ...
        'FontSize', 12, 'HorizontalAlignment', 'left');

    inputs(i) = uicontrol('Style', 'edit', 'String', num2str(adj_parms(i)), ...
        'Units', 'normalized', 'Position', [0.5 0.8 - (i-1) * 0.18, 0.3, 0.1], ...
        'FontSize', 12, 'Callback', @(src, ~) update_param(i, src));
end

% OK Button
uicontrol('Style', 'pushbutton', 'String', 'Save', ...
    'Units', 'normalized', 'Position', [0.35 0.05 0.3 0.15], ...
    'FontSize', 12, 'Callback', @okButtonCallback);


% OK Button
uicontrol('Style', 'pushbutton', 'String', 'Save', ...
    'Units', 'normalized', 'Position', [0.35 0.05 0.3 0.15], ...
    'FontSize', 12, 'Callback', @okButtonCallback);

%% Callback function to update parameter values
    function update_param(field, src)
        value = str2double(get(src, 'String'));
        if isnan(value) || isempty(value) || (value>=90) 
              set(src, 'String', num2str(adj_parms(field)));
            errordlg('Invalid input. Please enter a numeric value <90.', 'Error', 'modal');
        else        
            adj_parms(field) = value;
        end
    end

%% OK Button Callback and update parms to the struct before returning it
    function okButtonCallback(~, ~)

        parms.fas.range = adj_parms(1:2);
        parms.apo.super.maxangle = adj_parms(3);
        parms.apo.deep.maxangle= adj_parms(4);
        close(fig);
    end

% Wait for the figure to close before returning updated parameters
uiwait(fig);

end