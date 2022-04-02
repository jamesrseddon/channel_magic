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
% ROBUSTNESSSWITCH Calculate some robustness-like quantity based on R_ID
%  passed in.
function [ value,precision,status, dist ] = robustnessSwitch(R_ID,channel,...
                        A_matrix,solver_precision,solver_selection,...
                        varargin)


allowed_IDs = {'R_CHOI', 'R_CHOI_DIAG', 'CHANNEL_ROB',...
    'CHANNEL_ROB_DIAG', 'CAPACITY', 'CAPACITY_DIAG','R_CHOI_COMPLEX',...
    'GEN_ROB_CHOI'};

if ~any(strcmp(allowed_IDs,R_ID))
    error_struct.message = [R_ID ' is not a known identifier for a '...
        'robustness-type quantity. Allowed IDs are: '...
        sprintf('%s, ',allowed_IDs{:})];
    error_struct.identifier = ['quasi:robustnessSwitch:unknownQuantityID'];
    error(error_struct);
end
value = 0;
precision = 0;
status = 'unsolved';

calculated = false;


switch R_ID
    case 'R_CHOI'
        [optval,status,dist,dual_var,optbnd] =...
                            findChoiRobustness(channel,A_matrix,...
                                      solver_precision,solver_selection);
        calculated = true;
    case 'R_CHOI_COMPLEX'
        [optval,status,dist,dual_var,optbnd] =...
                            findChoiRobustness(channel,A_matrix,...
                             solver_precision,solver_selection,'complex');
        calculated = true;
    case 'R_CHOI_DIAG'
        [optval,status,dist,dual_var,optbnd] =...
                            findChoiRobustnessDiag(channel,A_matrix,...
                                      solver_precision,solver_selection);        
        calculated = true;
    case 'CHANNEL_ROB'
        [optval,status,dist,dual_var,optbnd] =...
                            findChoiCPTProbustness(channel,A_matrix,...
                                      solver_precision,solver_selection);
        calculated = true;
    case 'CHANNEL_ROB_DIAG'
        [optval,status,dist,dual_var,optbnd] =...
                            findChoiCPTProbustnessDiag(channel,A_matrix,...
                                      solver_precision,solver_selection);
        calculated = true;
    case 'CAPACITY'
        if nargin > 5
            param_string = varargin{1};
        else
            param_string = '';
        end
        logfilename = ['capacity_log_' param_string '.txt' ];
        [capacity, cap_results] = findCapacity(channel,A_matrix,A_matrix,...
                                      solver_precision,solver_selection,...
                                       logfilename);
        
        opt_val_vec = cell2mat(cap_results(2:end,1));
        opt_bnd_vec = cell2mat(cap_results(2:end,2));
        optval = capacity
        cap_index = find(opt_val_vec == capacity,1);
        optbnd = opt_bnd_vec(cap_index,1);
        calculated = true;
        status = [cap_results{(cap_index + 1),5} '(optimal channel)'];
        dist = cap_results{(cap_index + 1),6};
    case 'CAPACITY_DIAG'
        disp(['Capacity for diagonal states not yet coded.']);
    case 'GEN_ROB_CHOI'
        choi_state = makeChoiState(channel);
        [optval,dist,status] = RobMagG(choi_state,A_matrix);
        optbnd = optval;
end

value = mean([optval optbnd]);
precision = (optval - optbnd)/2;

end

