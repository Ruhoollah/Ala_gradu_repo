function [comp_match, method_data] = ctapeeg_recudetect_blink_ICs(EEG, varargin)
%CTAPEEG_RECUDETECT_BLINK_ICS measures the effect of removing each blink-like
% IC from the original data by projection.
% 
% Description: Similarity of each IC to the blink template
% is measured by eeglab_detect_icacomps_blinktemplate() and ICs are ordered 
% by this metric. Then the top Arg.test_pc of ICs are tested in turn: 
%  - blink-locked 1 sec ERPs are generated for the VEOG channel
%  - each IC is removed from the data by projection, and the value by which this 
%    operation moves the quantile band [2.5% 97.5%] toward zero is calculated.
%  - IF this 'toward-zero-movement' vector is significantly different to
%    zero by t-test, AND if its CIs are broader than 2uV (showing
%    meaningful amount of change), AND if the peak change occurs between 250ms 
%    and 750ms (change is blink-related and centre-locked)...
%    THEN the IC is labelled as bad
% 
% Syntax:
%   [comp_match, method_data] = ctapeeg_recudetect_blink_ICs(EEG, varargin)
% 
% Input:
%   EEG             struct, eeg data to test
%
% Varargin:
%   veog            string, channel name for VEOG, default = 'VEOG'
%   test_pc         scalar, percentage of ICs to test, default = 25
%   epoch_len_secs  scalar, half-length of epoch to centre on blink in seconds,
%                   default = 0.3, based on the estimate here - 
%           http://bionumbers.hms.harvard.edu/bionumber.aspx?s=y&id=100706&ver=0
% 
% Output:
%   comp_match      logical, vector of logicals indexing bad ICs
%   method_data     struct, fields are :
%                   'thArr' containing the blink template similarity of each IC, 
%                   'blinkERP' a dataframe with blink ERP of each channel
%
% See also: 
%   eeglab_detect_icacomps_blinktemplate()
%
%
% Version History:
% 10.01.2017 Created (Benjamin Cowley, FIOH)
%
% Copyright(c) 2015 FIOH:
% Benjamin Cowley (Benjamin.Cowley@ttl.fi), Jussi Korpela (jussi.korpela@ttl.fi)
%
% This code is released under the MIT License
% http://opensource.org/licenses/mit-license.php
% Please see the file LICENSE for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parse input arguments and set varargin defaults
p = inputParser;
p.addRequired('EEG', @isstruct);
p.addParameter('veog', 'VEOG', @iscell);
p.addParameter('test_pc', 25, @isscalar);
p.addParameter('epoch_len_secs', 0.3, @isscalar);

p.parse(EEG, varargin{:});
Arg = p.Results;
els = Arg.epoch_len_secs;


%% measure IC similarity to blink template
%Are blinks there? Fail if missing
if ( ~ismember('blink', unique({EEG.event.type})) )
    error('ctapeeg_recudetect_blink_ICs:noBlinkEvents',...
        ['No events of type ''blink'' found. Cannot proceed with'...
        ' blink template matching. Run CTAP_blink2event() to fix.'])
end

veogidx = ismember({EEG.chanlocs.labels}, Arg.veog);
if sum(veogidx) > 2
    error('ctapeeg_recudetect_blink_ICs:tooManyVeogChans',...
        '%d channels match the given VEOG info: can''t proceed.', sum(veogidx))
end    

% Detect components
[~, th_arr] = eeglab_detect_icacomps_blinktemplate(EEG, 'leqThreshold', 0);
comp_match = false(numel(th_arr), 1);


%% Compute blink ERP for comparison to detected ICs
EEGbl = pop_epoch(EEG, {'blink'}, [-els, els]);
EEGbl = pop_rmbase(EEGbl, [-els * 1000, 0]);
blERPdf = create_dataframe(mean(EEGbl.data, 3),...
    {'channel', 'time'},...
    {{EEGbl.chanlocs.labels}, EEGbl.times});

%get VEOG vector from one or more VEOG channels, create blink ERP + bands
if sum(veogidx) == 2
    blinkERP = diff(blERPdf.data(veogidx, :));
    veogdat = squeeze(diff(EEGbl.data(veogidx, :, :)));
else
    blinkERP = blERPdf.data(veogidx, :);
    veogdat = squeeze(EEGbl.data(veogidx, :, :));
end
blinkERP = sbf_make_ERPband(blinkERP, veogdat, 0.025);

%test if detected IC removal gives clean blink ERP
[~, test_idx] = sort(th_arr); %ascending sort: lower is better
for i = 1:ceil(numel(th_arr) * Arg.test_pc / 100)%consider only first X% ICs.
    %IC_match(test_idx(i)) = sbf_test_IC_removal(test_idx(i));
    comp_match(test_idx(i)) = sbf_test_IC_ALONE(test_idx(i));
end

method_data.thArr = th_arr;
method_data.blinkERP = blERPdf;


    %% Subfunctions
    % test blinkIC against given blink template IC
    function isblink = sbf_test_IC_ALONE(blinkIC)

        isblink = false;

        %Compute blinkIC-subtracted ERP for comparison
        EEGic = pop_subcomp(EEG, blinkIC);
        EEGic = pop_select(EEGic, 'channel', Arg.veog);
        EEGic = pop_epoch(EEGic, {'blink'}, [-els, els]);
        EEGic = pop_rmbase(EEGic, [-els * 1000, 0]);
        if sum(veogidx) == 2
            EEGic.data = diff(EEGic.data);
        end
        icERP = ...
            sbf_make_ERPband(mean(EEGic.data, 3), squeeze(EEGic.data), 0.025);

        % DEBUG
        % {
        %figure; subplot(1, 2, 1); plot(blinkERP'); subplot(1, 2, 2); plot(icERP')
        %pause(2)
        %close
        %}
        %numerically compare ERPs
        [corrected_erprop_min, corrected_erprop_mean, corrected_erprop_max] = sbf_get_erp_proportion(blinkERP, icERP);
         metric_min=  find(corrected_erprop_min >=0.5);
         metric_mean= find(corrected_erprop_mean >=0.5);
         metric_max=  find(corrected_erprop_max >=0.5);

        if numel(metric_min)+ numel (metric_mean)+ numel (metric_max) >= 75
        isblink = true;
        end
  
    end


    % In order to check if an IC is a blink-related IC, we need to compare the similarity between the IC ERP 
 % and Blink ERP. first we check if the coontribution of each component of IC ERP is more than 50%. 
 % to perform such task we attribute a proportional value to each componet OF IC ERP showing the extent of 
 % their contribution to the Blink ERP (Note that this value cannot be larger than one as there is no IC ERP 
 % which could contribute to Blink ERP by 100%). Min, mean and max extentions are representing the three chosen 
 % quantiles. After finding the proper values of the three ERPs, those with proportinal value bigger than 0.5 
 % will be considered as potential Bad IC cases. And  after finding the significant contributors, we need to check
 % HOW MANY of the ERP elements are having this characteristic, that is proprotional value >= 0.5. 
 % of course not all of the Blink ERP is interesting to us, but  we are more interested in significant 
 % changes/variations indicating the blink occurence. By inspecting all the blink ERPs (min&mean&max) we can see 
 % that the range of considerable variations  is between 100 and 250. in other words, it takes 150 elements 
 % for a blink to happen. Since ICA is an stochastic process, it would suffic to consider only half of this amount
 % (75 elements)to extract blink-related ICs. Thus, if ERP of an IC can fullfil the term that at least 75 elements
 % of that ERP should have a 50% contribution to the blink ERP, in this way that IC can be marked as blink-related 
 % or bad IC.
 
    % helper function to calculate the proportions of min & mean & max ERPs 
    function [corrected_erprop_min, corrected_erprop_mean, corrected_erprop_max] = sbf_get_erp_proportion(erq1, erq2) 
    erp_prop_min= ([abs(erq2(1,:)) ./ abs(erq1(1,:))]);
    wrong_min=erp_prop_min >=1;
    faulty_min= erp_prop_min (wrong_min);
    corrected_erprop_min = setdiff (erp_prop_min, faulty_min);


    erp_prop_mean= ([abs(erq2(2,:)) ./ abs(erq1(2,:))]);
    wrong_mean=erp_prop_mean >=1;
    faulty_mean= erp_prop_mean (wrong_mean);
    corrected_erprop_mean = setdiff (erp_prop_mean, faulty_mean);


    erp_prop_max= ([abs(erq2(3,:)) ./ abs(erq1(3,:))]);
    wrong_max=erp_prop_max >=1;
    faulty_max= erp_prop_max (wrong_max);
    corrected_erprop_max = setdiff (erp_prop_max, faulty_max);
   
    end

    % Helper function to compare quantiles by how much they move toward zero
    %function [erqdiff, h, p, ci] = sbf_compare_erp_quantiles(erq1, erq2)
        % working with [3 n] vectors of big quantiles, mean, little quantiles,
        % e.g. [2.5% M 97.5%]
        % measure how much erq2 moves the quantile band toward 0, i.e. denoises
        %erqdiff = abs(diff([abs(erq1(1, :)); abs(erq1(3, :))])) -...
            %abs(diff([abs(erq2(1, :)); abs(erq2(3, :))]));
        %[h, p, ci] = ttest(erqdiff,0,'alpha', 0.00005);
        %MAYBEDO - USE A TEST THAT ACCOMMODATES THE AUTOCORRELATION OF THIS
        %HIGHLY NON-INDEPENDENT SET OF OBSERVATIONS.
    %end

    function ERP = sbf_make_ERPband(ERPdata, qdata, qpc)
        ERPdata = ERPdata(:)';
        cols = size(ERPdata, 2);
        [r, c] = size(qdata);
        if ~any(cols == [r c])
            error('sbf_make_ERPband:dimension_mismatch', 'Sthg terribly wrong!')
        end
        if r == c
            dim = 2;
        else
            dim = find(cols ~= [r c]);
        end
        q1 = quantile(qdata, qpc, dim);
        q2 = quantile(qdata, 1 - qpc, dim);
        ERP = [q1(:)'; ERPdata; q2(:)'];
    end
    
end %ctapeeg_recudetect_blink_ICs()
  
