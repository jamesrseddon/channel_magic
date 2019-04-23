% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
%
%  Copyright © 2019 James R. Seddon
%  This file is part of project: Channel Magic
%
function [ summary_result_array,detail_result_array ] =...
    noisyChannelParamSweep(channel_type,param_names,parameters,...
                robustness_config,varargin)
% NOISYCHANNELPARAMSWEEP Sweep parameters for a noisy channel and obtain
% robustness results
%   channel_type = string: specifies channel as per constructChannel
%   function, eg. 'AMP_X' is amplitude damping followed by X-rotation.
%
%   param_names = cell aray: names of parameters required for the specified
%   channel_type. eg. for 'AMP_X' this should be {'p','theta'}.
%
%   parameters  = cell array. Each cell contains a row vector of numerical
%   values. eg.{[0.1, 0.01, 0.001], [1:100]*0.01*pi/4}. First vector can
%   have any number of entries. The other vectors should all be of the same
%   length.
%
%   robustness_config = struct: specifies which quantities to calculate
%   for each instance, as well as containing the A matrices required for
%   the linear program.
%       robustness_config.optimisations = cell array giving IDs for quantities
%       to calculate, indices for A matrices, and label for results column.
%       eg. {'R_CHOI', 1, 'Robustness of Choi state';
%            'R_CHOI', 2, 'R_CPR'}
%       
%       robustness_config.A_matrices = cell array: contains A matrices for
%       the optimisation problems specified by .optimisations, eg:
%       {A_mat_2; A_mat_CP_1};
%
%   varargin{1} can be a string. This then becomes the filename for a text
%   file where the results will be saved.

%%%% INPUT VALIDATION %%%%
if ~ischar(channel_type)
    error_struct.message = ['channel_type must be a string.'];                
    error_struct.identifier = ['quasi:noisyChannelParamSweep:'...
                'invalidChannelType'];
    error(error_struct)
end

[row_check_names,num_param_names] = size(param_names);

if (row_check_names > 1)||(~iscell(param_names))
    error_struct.message = ['param_names should be a cell array with '...
        'format: {''param_1'', ''param_2'',...}'];                
    error_struct.identifier = ['quasi:noisyChannelParamSweep:'...
                'invalidParamNames'];
    error(error_struct)
end

for kk=1:num_param_names
    if ~ischar(param_names{1,kk})
        error_struct.message = ['All cells of param_names should '...
            'contain strings.'];
        error_struct.identifier = ['quasi:noisyChannelParamSweep:'...
                'paramNameNotString'];
            error(error_struct)
    end
end

[row_check_params,num_params] = size(parameters);

if (row_check_params > 1)||(~iscell(parameters))
    error_struct.message = ['parameters should be a cell array with '...
        'format: {params_1, params_2,...}'];                
    error_struct.identifier = ['quasi:noisyChannelParamSweep:'...
                'invalidParameterFormat'];
    error(error_struct)
end

if num_params ~= num_param_names
    error_struct.message = ['There are ' num2str(num_param_names)...
        ' parameter names, but ' num2str(num_params) ' parameters.' ];                
    error_struct.identifier = ['quasi:noisyChannelParamSweep:'...
                'parameterCountMismatch'];
    error(error_struct)
end

for kk=1:num_params
    if ~isnumeric(parameters{1,kk})
        error_struct.message = ['All cells of parameters should '...
            'contain numbers.'];
        error_struct.identifier = ['quasi:noisyChannelParamSweep:'...
                'parameterNotNumeric'];
        error(error_struct)
    end
end

% Check for file path specified in argument.

writefile_flag = false;

if nargin > 4
    writefile_flag = true;
    result_file_path = [pwd '\' varargin{1}];
end

%%%% MAIN CODE %%%%


first_param = parameters{1};
first_param_name = param_names{1};

if num_params > 1
    multi_params_flag = true;
    num_values_2 = length(parameters{2});
else
    multi_params_flag = false;
    num_values_2 = 1;
end

num_values_1 = length(first_param);

rob_config_rows = size(robustness_config.optimisations,1);

%%% build header row of summary results table %%%

header_row = {'Channel type'};


% parameter columns
for mm = 1:num_params
    header_row = [header_row, param_names{mm}];
end
header_row_file = header_row;

% Result columns
for nn = 1:rob_config_rows
    quantity_name = robustness_config.optimisations{nn,3};
    val_header = [quantity_name ' value'];
    precision_header = [quantity_name ' precision'];
    status_header = [quantity_name ' status'];
    distrib_header = [quantity_name ' decomposition'];
    header_row = [header_row,val_header,precision_header,...
        status_header, distrib_header];
    header_row_file = [header_row_file,val_header,precision_header,...
        status_header];
end 
results = header_row

% open a file and write the header row
if writefile_flag
    fileid = fopen(result_file_path, 'a') ;
    fprintf(fileid, '"%s",', header_row_file{1,1:end-1}) ;
    fprintf(fileid, '"%s"\n', header_row_file{1,end}) ;
    fclose(fileid) ;
end

for pp = 1:num_values_1 % Loop over indices of first parameter sweep
    param_struct = struct(first_param_name,first_param(pp));
    init_param_string = [first_param_name '_' num2str(first_param(pp))];
    for qq = 1:num_values_2 % Loop over indices of subsequent sweeps
        result_row = {channel_type,first_param(pp)};
        param_string = init_param_string;
        % Pick additional parameters for this run        
        if multi_params_flag
            for kk = 2:num_params
                these_parameters = parameters{kk};
                this_parameter_val = these_parameters(qq);
                param_struct = setfield(param_struct,param_names{kk},...
                    this_parameter_val);
                result_row = [result_row, this_parameter_val];
                param_string = [param_string '_' param_names{kk} '_'...
                                num2str(this_parameter_val)];
            end
            result_row_file = result_row;
        end
        param_string
        param_struct
        construct_success = true;
        try
            kraus_ops = constructChannel(channel_type,param_struct);
        catch error_struct
            warning(['Error generating Kraus operators. Skipping this '...
                'set of parameters. Error thrown was: '...
                error_struct.message]);
            channel_error_ID = error_struct.identifier;
            construct_success = false;
            kraus_ops = 0;
        end
        
                
        for rr = 1:rob_config_rows % loop over quantities to be calculated.
                this_config =  robustness_config.optimisations(rr,:);
                [R_ID, A_mat_index, ~] = this_config{:};
                if length(this_config) > 3
                    solver_precision = this_config{1,4}
                else
                    solver_precision = 'default';
                end
                if length(this_config) > 4
                    solver_selection = this_config{1,5}
                else
                    solver_selection = 'default';
                end
                A_matrix = robustness_config.A_matrices{A_mat_index,1};
                if construct_success
                    % Try to calculate the robustness, catch the error and
                    % skip if something goes wrong.
                    try
                        [rob_val,rob_precision,solve_status, distrib] = ...
                            robustnessSwitch(R_ID,kraus_ops,A_matrix,...
                                solver_precision,solver_selection,...
                                param_string);
                    catch error_struct
                        warning(['Error calculating robustness-like '...
                            'quantity. Skipping. Error thrown was: '...
                            error_struct.message '.']);
                        rob_val = 0;
                        rob_precision = 0;
                        solve_status = ['error: ' error_struct.identifier];
                        distrib = [];
                    end
                    result_row = [result_row,...
                        rob_val,rob_precision,solve_status,distrib];
                    result_row_file = [result_row_file,...
                        rob_val,rob_precision,solve_status];
                else
                    error_status = ['channel error: ' channel_error_ID];
                    result_row = [result_row,0,0,error_status,[]];
                end
        end
        results = [results; result_row];
        
        
        % write result row to the output file
        if writefile_flag
            fileid = fopen(result_file_path, 'a') ;
            fprintf(fileid, '%s,', result_row_file{1,1}) ;
            fprintf(fileid, '%.8f,', result_row_file{1,2:(num_params + 1)});
            fprintf(fileid, '%.8f,%.8e,"%s",',...
                result_row_file{1,(num_params + 2):end});
            fprintf(fileid, '\n');
            fclose(fileid);
        end
    end
end

summary_result_array = results;
detail_result_array = 'placeholder';
    