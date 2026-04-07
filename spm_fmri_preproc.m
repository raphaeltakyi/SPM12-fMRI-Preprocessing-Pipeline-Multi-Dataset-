%% SPM12 Preprocessing Script
% Datasets: Multiple datasets
% TR: Set specific TR for each dataset
% Slices: detected automatically per session from NIfTI header
% Steps: Realign+Unwarp > Slice Timing > Coregister > Segment > Normalise > Smooth

%% ========== PATHS ==========

dataset1_dir   = 'path/to/dataset1';
dataset2_dir    = 'path/to/dataset2';
output_dir = 'path/to/output';

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ========== SHARED PARAMETERS ==========

fwhm     = 6;            % smoothing kernel mm
vox_size = [2 2 2];      % normalised voxel size mm
bb       = [-78 -112 -70; 78 76 85];  % MNI bounding box

% Slice order *** CONFIRM WITH MRI PROTOCOL SHEET ***
% Options: 'ascending', 'descending', 'interleaved_ascending', 'interleaved_descending'
slice_order = 'ascending';

% TR per dataset (confirmed from NIfTI headers)
TR_dataset1 = XX;
TR_dataset2  = xx;

%% ========== INITIALISE SPM ==========

spm('defaults', 'FMRI');
spm_jobman('initcfg');

%% ========== COLLECT ALL SUBJECTS ==========

fprintf('=== Collecting subjects ===\n');

dataset1_subs = dir(fullfile(dataset1_dir, 'sub-*'));
dataset1_subs = dataset1_subs([dataset1_subs.isdir]);
dataset1_subs = dataset1_subs(arrayfun(@(x) ~startsWith(x.name, '.'), dataset1_subs));

dataset2_subs = dir(fullfile(dataset2_dir, 'dataset2_sub_*'));
dataset2_subs = dataset2_subs([dataset2_subs.isdir]);
dataset2_subs = dataset2_subs(arrayfun(@(x) ~startsWith(x.name, '.'), dataset2_subs));

all_subjects = {};

for s = 1:length(dataset1_subs)
    all_subjects{end+1} = struct(...
        'name',        dataset1_subs(s).name, ...
        'dataset',     'dataset1',            ...
        'base',        dataset1_dir,          ...
        'ses_pat',     'ses-*',           ...
        't1_pat',      '*T1w.nii',        ...
        'TR',          TR_dataset1,           ...
        'slice_order', slice_order);
end

for s = 1:length(dataset2_subs)
    all_subjects{end+1} = struct(...
        'name',        dataset2_subs(s).name,  ...
        'dataset',     'dataset2',             ...
        'base',        dataset2_dir,           ...
        'ses_pat',     'session_*',       ...
        't1_pat',      '*T1.nii',         ...
        'TR',          TR_dataset2,            ...
        'slice_order', slice_order);
end

fprintf('Total: %d subjects (%d dataset1 + %d dataset2)\n\n', ...
         length(all_subjects), length(dataset1_subs), length(dataset2_subs));

%% ========== PREPROCESSING LOG ==========

log_file = fullfile(output_dir, 'preprocessing_log.txt');
fid_log  = fopen(log_file, 'w');
fprintf(fid_log, 'SPM12 Preprocessing Log\n');
fprintf(fid_log, 'Started: %s\n\n', datestr(now));

%% ========== MAIN PREPROCESSING LOOP ==========

for s = 1:length(all_subjects)
    subj    = all_subjects{s};
    sub_dir = fullfile(subj.base, subj.name);

    fprintf('\n============================================================\n');
    fprintf('Subject %d/%d: %s (%s)\n', s, length(all_subjects), subj.name, subj.dataset);
    fprintf('============================================================\n');

    % ── Get sessions ─────────────────────────────────────────────────────
    ses_dirs = dir(fullfile(sub_dir, subj.ses_pat));
    ses_dirs = ses_dirs([ses_dirs.isdir]);
    ses_dirs = ses_dirs(arrayfun(@(x) ~startsWith(x.name, '.'), ses_dirs));
    ses_dirs = ses_dirs(1:min(2, end));

    if isempty(ses_dirs)
        msg = sprintf('WARNING: No sessions for %s — skipping\n', subj.name);
        fprintf(msg); fprintf(fid_log, msg);
        continue;
    end

    % ── Get structural ────────────────────────────────────────────────────
    anat_dir = fullfile(sub_dir, ses_dirs(1).name, 'anat');
    t1_files = dir(fullfile(anat_dir, subj.t1_pat));
    t1_files = t1_files(arrayfun(@(x) ~startsWith(x.name, '._'), t1_files));

    if isempty(t1_files)
        msg = sprintf('WARNING: No T1 for %s — skipping\n', subj.name);
        fprintf(msg); fprintf(fid_log, msg);
        continue;
    end

    t1_path = fullfile(anat_dir, t1_files(1).name);
    fprintf('  Structural: %s\n', t1_files(1).name);

    % ── Collect 4D functional files + detect slices per session ──────────
    all_session_vols = {};
    session_nslices  = [];   % n_slices detected per session from header
    session_TR       = [];   % TR per session

    for ss = 1:length(ses_dirs)
        ses_name = ses_dirs(ss).name;
        func_dir = fullfile(sub_dir, ses_name, 'func');

        bold = dir(fullfile(func_dir, '*bold*.nii'));
        bold = bold(arrayfun(@(x) ~startsWith(x.name, '._'), bold));
        bold = bold(~[bold.isdir]);
        bold = bold([bold.bytes] > 100e6);   % only 4D files

        if isempty(bold)
            msg = sprintf('  WARNING: No 4D bold for %s/%s\n', subj.name, ses_name);
            fprintf(msg); fprintf(fid_log, msg);
            all_session_vols{ss} = '';
            session_nslices(ss)  = NaN;
            session_TR(ss)       = subj.TR;
            continue;
        end

        bold_path = fullfile(func_dir, bold(1).name);

        % Read header to get actual slice count for THIS session
        V_hdr               = spm_vol(bold_path);
        n_sl                = V_hdr(1).dim(3);   % z dimension = n_slices
        session_nslices(ss) = n_sl;
        session_TR(ss)      = subj.TR;

        all_session_vols{ss} = bold_path;
        fprintf('  Session %d: %s | %d slices | TR=%.4fs\n', ...
                 ss, bold(1).name, n_sl, subj.TR);
    end

    % ================================================================
    % STEP 1 — REALIGN AND UNWARP (one session at a time)
    % Sessions are processed separately to handle different slice counts
    % ================================================================
    fprintf('\n  >> Step 1: Realign and Unwarp\n');

    unwarped_vols = {};
    mean_func     = '';

    for ss = 1:length(all_session_vols)
        if isempty(all_session_vols{ss}) || strcmp(all_session_vols{ss}, '')
            unwarped_vols{ss} = '';
            continue;
        end

        fprintf('     Session %d (%d slices)...\n', ss, session_nslices(ss));

        vols = expand_4d(all_session_vols{ss});

        matlabbatch = {};
        matlabbatch{1}.spm.spatial.realignunwarp.data(1).scans  = vols;
        matlabbatch{1}.spm.spatial.realignunwarp.data(1).pmscan = {''};

        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality    = 0.9;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep        = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm       = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm        = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp    = 2;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap      = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight     = {''};

        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn   = [8 8];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda   = 100000;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm       = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot      = [4 5];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot      = [];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm   = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem      = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi      = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';

        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich  = [2 1];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp  = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap     = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask     = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix   = 'u';

        try
            spm_jobman('run', matlabbatch);

            [fp, fn, fe]      = fileparts(all_session_vols{ss});
            unwarped_vols{ss} = fullfile(fp, ['u' fn fe]);

            % Save mean functional from first valid session for coregistration
            if isempty(mean_func)
                mean_func = fullfile(fp, ['meanu' fn '.nii']);
            end

            fprintf('     Session %d done\n', ss);

        catch ME
            msg = sprintf('  WARNING: Realign failed %s ses-%d: %s\n', ...
                           subj.name, ss, ME.message);
            fprintf(msg); fprintf(fid_log, msg);
            unwarped_vols{ss} = '';
        end
    end

    % Check at least one session succeeded
    valid_unwarped = ~cellfun(@(x) isempty(x) || strcmp(x,''), unwarped_vols);
    if ~any(valid_unwarped)
        msg = sprintf('WARNING: All sessions failed realign for %s — skipping\n', subj.name);
        fprintf(msg); fprintf(fid_log, msg);
        continue;
    end

    fprintf('  >> Step 1 done\n');

    % ================================================================
    % STEP 2 — SLICE TIMING CORRECTION (one session at a time)
    % Uses n_slices detected from header for each session
    % ================================================================
    fprintf('\n  >> Step 2: Slice Timing Correction\n');

    st_vols = {};

    for ss = 1:length(unwarped_vols)
        if isempty(unwarped_vols{ss}) || strcmp(unwarped_vols{ss}, '')
            st_vols{ss} = '';
            continue;
        end

        n_sl  = session_nslices(ss);
        TR_ss = session_TR(ss);
        r_sl  = round(n_sl / 2);

        fprintf('     Session %d: n_slices=%d | ref=%d | TR=%.4f\n', ...
                 ss, n_sl, r_sl, TR_ss);

        s_order = build_slice_order(subj.slice_order, n_sl);
        vols    = expand_4d(unwarped_vols{ss});

        matlabbatch = {};
        matlabbatch{1}.spm.temporal.st.scans    = {vols};
        matlabbatch{1}.spm.temporal.st.nslices  = n_sl;
        matlabbatch{1}.spm.temporal.st.tr       = TR_ss;
        matlabbatch{1}.spm.temporal.st.ta       = TR_ss - (TR_ss / n_sl);
        matlabbatch{1}.spm.temporal.st.so       = s_order;
        matlabbatch{1}.spm.temporal.st.refslice = r_sl;
        matlabbatch{1}.spm.temporal.st.prefix   = 'a';

        try
            spm_jobman('run', matlabbatch);

            [fp, fn, fe]  = fileparts(unwarped_vols{ss});
            st_vols{ss}   = fullfile(fp, ['a' fn fe]);
            fprintf('     Session %d done\n', ss);

        catch ME
            msg = sprintf('  WARNING: Slice timing failed %s ses-%d: %s\n', ...
                           subj.name, ss, ME.message);
            fprintf(msg); fprintf(fid_log, msg);
            st_vols{ss} = '';
        end
    end

    fprintf('  >> Step 2 done\n');

    % ================================================================
    % STEP 3 — COREGISTRATION
    % Aligns T1 to mean functional from session 1
    % ================================================================
    fprintf('\n  >> Step 3: Coregistration\n');

    if ~exist(mean_func, 'file')
        % Fallback to first volume of first valid session
        first_st = st_vols{find(~cellfun(@(x) isempty(x)||strcmp(x,''), st_vols), 1)};
        vols      = expand_4d(first_st);
        mean_func = vols{1};
        fprintf('     Using first volume as reference (mean not found)\n');
    end

    matlabbatch = {};
    matlabbatch{1}.spm.spatial.coreg.estimate.ref               = {mean_func};
    matlabbatch{1}.spm.spatial.coreg.estimate.source            = {t1_path};
    matlabbatch{1}.spm.spatial.coreg.estimate.other             = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = ...
        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];

    try
        spm_jobman('run', matlabbatch);
        fprintf('  >> Step 3 done\n');
    catch ME
        msg = sprintf('  WARNING: Coregistration failed %s: %s\n', subj.name, ME.message);
        fprintf(msg); fprintf(fid_log, msg);
        continue;
    end

    % ================================================================
    % STEP 4 — SEGMENTATION
    % Segments T1 into GM, WM, CSF + estimates deformation field
    % ================================================================
    fprintf('\n  >> Step 4: Segmentation\n');

    tpm_path = fullfile(spm('dir'), 'tpm', 'TPM.nii');

    matlabbatch = {};
    matlabbatch{1}.spm.spatial.preproc.channel.vols     = {t1_path};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 1];

    ngaus_vals = [1, 1, 2, 3, 4, 2];
    for c = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(c).tpm   = ...
            {[tpm_path, ',', num2str(c)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(c).ngaus = ngaus_vals(c);
        if c <= 3
            matlabbatch{1}.spm.spatial.preproc.tissue(c).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(c).warped = [0 0];
        else
            matlabbatch{1}.spm.spatial.preproc.tissue(c).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(c).warped = [0 0];
        end
    end

    matlabbatch{1}.spm.spatial.preproc.warp.mrf        = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup    = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg        = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg     = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm       = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp       = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write      = [0 1];

    try
        spm_jobman('run', matlabbatch);
        fprintf('  >> Step 4 done\n');
    catch ME
        msg = sprintf('  WARNING: Segmentation failed %s: %s\n', subj.name, ME.message);
        fprintf(msg); fprintf(fid_log, msg);
        continue;
    end

    % Get forward deformation field
    [anat_path_only, t1_name, ~] = fileparts(t1_path);
    def_field = fullfile(anat_path_only, ['y_' t1_name '.nii']);

    if ~exist(def_field, 'file')
        msg = sprintf('WARNING: Deformation field not found for %s\n', subj.name);
        fprintf(msg); fprintf(fid_log, msg);
        continue;
    end

    % ================================================================
    % STEP 5 — NORMALISATION
    % Warps all functionals + structural to MNI space
    % ================================================================
    fprintf('\n  >> Step 5: Normalisation\n');

    all_vols_to_norm = {};
    for ss = 1:length(st_vols)
        if ~isempty(st_vols{ss}) && ~strcmp(st_vols{ss}, '')
            vols             = expand_4d(st_vols{ss});
            all_vols_to_norm = [all_vols_to_norm; vols];
        end
    end

    % Include bias-corrected structural
    bias_corr_t1 = fullfile(anat_path_only, ['m' t1_name '.nii']);
    if exist(bias_corr_t1, 'file')
        all_vols_to_norm{end+1} = bias_corr_t1;
    end

    matlabbatch = {};
    matlabbatch{1}.spm.spatial.normalise.write.subj.def        = {def_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample   = all_vols_to_norm;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb     = bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = vox_size;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

    try
        spm_jobman('run', matlabbatch);
        fprintf('  >> Step 5 done\n');
    catch ME
        msg = sprintf('  WARNING: Normalisation failed %s: %s\n', subj.name, ME.message);
        fprintf(msg); fprintf(fid_log, msg);
        continue;
    end

    % Get normalised functional paths
    norm_vols = {};
    for ss = 1:length(st_vols)
        if ~isempty(st_vols{ss}) && ~strcmp(st_vols{ss}, '')
            [fp, fn, fe]  = fileparts(st_vols{ss});
            norm_vols{ss} = fullfile(fp, ['w' fn fe]);
        else
            norm_vols{ss} = '';
        end
    end

    % ================================================================
    % STEP 6 — SMOOTHING
    % Final step — output swau*.nii files ready for CONN
    % ================================================================
    fprintf('\n  >> Step 6: Smoothing (FWHM=%dmm)\n', fwhm);

    all_norm_vols = {};
    for ss = 1:length(norm_vols)
        if ~isempty(norm_vols{ss}) && ~strcmp(norm_vols{ss}, '')
            vols          = expand_4d(norm_vols{ss});
            all_norm_vols = [all_norm_vols; vols];
        end
    end

    matlabbatch = {};
    matlabbatch{1}.spm.spatial.smooth.data   = all_norm_vols;
    matlabbatch{1}.spm.spatial.smooth.fwhm   = [fwhm fwhm fwhm];
    matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{1}.spm.spatial.smooth.im     = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

    try
        spm_jobman('run', matlabbatch);
        fprintf('  >> Step 6 done\n');
    catch ME
        msg = sprintf('  WARNING: Smoothing failed %s: %s\n', subj.name, ME.message);
        fprintf(msg); fprintf(fid_log, msg);
        continue;
    end

    msg = sprintf('COMPLETE: %s (%s) — %s\n', subj.name, subj.dataset, datestr(now));
    fprintf('\n  %s', msg);
    fprintf(fid_log, msg);

end

%% ========== DONE ==========

fclose(fid_log);

fprintf('\n\n============================================================\n');
fprintf('  ALL SUBJECTS COMPLETE\n');
fprintf('  Final files prefix: swau\n');
fprintf('    s = smoothed\n');
fprintf('    w = normalised to MNI\n');
fprintf('    a = slice timing corrected\n');
fprintf('    u = realigned + unwarped\n');
fprintf('  Log: %s\n', log_file);
fprintf('============================================================\n');

%% ========== HELPER FUNCTIONS ==========

function vols = expand_4d(fname)
    V    = spm_vol(fname);
    vols = cell(length(V), 1);
    for i = 1:length(V)
        vols{i} = sprintf('%s,%d', fname, i);
    end
end

function s_order = build_slice_order(order_name, n_slices)
    switch lower(order_name)
        case 'ascending'
            s_order = 1:n_slices;
        case 'descending'
            s_order = n_slices:-1:1;
        case 'interleaved_ascending'
            if mod(n_slices, 2) == 0
                s_order = [2:2:n_slices, 1:2:n_slices];
            else
                s_order = [1:2:n_slices, 2:2:n_slices];
            end
        case 'interleaved_descending'
            if mod(n_slices, 2) == 0
                s_order = [n_slices:-2:2, n_slices-1:-2:1];
            else
                s_order = [n_slices:-2:1, n_slices-1:-2:2];
            end
        otherwise
            warning('Unknown slice order: %s — defaulting to ascending', order_name);
            s_order = 1:n_slices;
    end
end