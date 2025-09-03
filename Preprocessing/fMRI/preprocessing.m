function [ok,matlabbatch,outputfiles,job_id]=conn_setup_preproc(STEPS,varargin)
% CONN_SETUP_PREPROC
% Runs individual preprocessing steps
%
% conn_setup_preproc(steps)
% runs preprocessing pipeline (default_*) or one/multiple individual preprocessing steps (structural_* and functional_*). Valid step names are (enter as cell array to run multiple sequential steps):
%
% conn_setup_preproc(steps,'param1_name',param1_value,'param2_name',param2_value,...)
% defines additional non-default values for parameters specific to individual steps to be used in this preprocessing run
%
% conn_setup_preproc('settings', 'param1_name', param1_value, 'param2_name', param2_value, ...)
% defines additional non-default values for parameters specific to individual steps to be used in ALL future calls to conn_setup_preproc (during the current Matlab session, or until a new conn_setup_preproc('settings',...) command) (e.g. use "conn_setup_preproc('settings')" to clear previous values)
%
% conn_setup_preproc('steps')
% returns the full list of valid preprocessing-step names
%
% see "help conn_batch" for details about available preprocessing steps as well as a full list of additional preprocessing parameters
%


global CONN_x CONN_gui;
persistent fixed_options;

PREFERSPM8OVERSPM12=false; % set to true if you prefer to use SPM8 procedures over SPM12 ones (when SPM12 is installed)
ALLSETSPERMISSIONS=false;  % set to true if you want to allow dataset-1 or above preprocessing steps to import new ROIs, first-level covariates, or second-level covariates into your project
if isdeployed, spmver12=true;
else spmver12=str2double(regexp(spm('ver'),'SPM(\d+)','tokens','once'))>=12;
end
if isfield(CONN_gui,'font_offset'),font_offset=CONN_gui.font_offset; else font_offset=0; end
%if isfield(CONN_x,'pobj')&&isfield(CONN_x.pobj,'readonly')&&CONN_x.pobj.readonly, error('This procedure cannot be run while in view-only mode. Please re-load your project to enable edits'); end
if ~nargin, STEPS=''; varargin={'multiplesteps',1}; end
if isequal(STEPS,'settings'), fixed_options=varargin; return; end
if isempty(fixed_options), options=varargin;
else options=[fixed_options, varargin];
end
steps={'default_mni','default_mnifield','default_mnidirectfield','default_ss','default_ssfield','default_ssnl',...
    'functional_surface_coreg&resample',...
    'structural_manualorient','structural_center','structural_segment',...
    'structural_normalize','structural_segment&normalize',...
    'structural_segment&normalize_withlesion',...
    'functional_segment&normalize_indirect_withlesion',...
    'structural_normalize_preservemasks',...
    'structural_manualspatialdef', ...
    'structural_mask',...
    'functional_removescans','functional_manualorient','functional_center','functional_centertostruct'...
    'functional_slicetime','functional_bandpass','functional_regression','functional_roiextract','functional_mask','functional_realign','functional_realign&unwarp',...
    'functional_realign&unwarp&fieldmap','functional_art','functional_coregister_affine_reslice',...
    'functional_segment',...
    'functional_manualspatialdef',...
    'functional_smooth','functional_motionmask',...
    'functional_segment&normalize_indirect','functional_normalize_indirect','functional_normalize_indirect_preservemasks', ...
    'functional_segment&normalize_direct','functional_normalize_direct', ...
    'functional_realign_noreslice', ...
    'functional_coregister_nonlinear', ...
    'functional_coregister_affine_noreslice', ...
    'functional_label', ...
    'functional_label_as_original', ...
    'functional_label_as_subjectspace', ...
    'functional_label_as_mnispace', ...
    'functional_label_as_surfacespace', ...
    'functional_label_as_smoothed', ...
    'functional_load', ...
    'functional_load_from_original', ...
    'functional_load_from_subjectspace', ...
    'functional_load_from_mnispace', ...
    'functional_load_from_surfacespace', ...
    'functional_load_from_smoothed', ...
    'functional_smooth_masked',...
    'functional_surface_resample', ...
    'functional_surface_smooth', ...
    'functional_vdm_apply', ...
    'functional_vdm_create' ...
    };
%'functional_normalize','functional_segment&normalize',...
steps_names={'<HTML><b>default preprocessing pipeline</b> for volume-based analyses (direct normalization to MNI-space)</HTML>','<HTML><b>preprocessing pipeline</b> for volume-based analyses (indirect normalization to MNI-space) when FieldMaps are available</HTML>','<HTML><b>preprocessing pipeline</b> for volume-based analyses (direct normalization to MNI-space) when FieldMaps are available</HTML>','<HTML><b>preprocessing pipeline</b> for surface-based analyses (in subject-space)</HTML>','<HTML><b>preprocessing pipeline</b> for surface-based analyses (in subject-space) when FieldMaps are available</HTML>','<HTML><b>preprocessing pipeline</b> for surface-based analyses (in subject-space) using nonlinear coregistration</HTML>',...
    'functional Direct coregistration to structural followed by resampling of functional data within the cortical surface (converts volume- to surface- level data using FreeSurfer subject-specific surfaces)', ...
    'structural Manual transformation (rotation/flip/translation/affine of structural volumes)','structural Center to (0,0,0) coordinates (translation)','structural Segmentation (Grey/White/CSF tissue estimation)',...
    'structural Normalization (MNI space normalization)',...
    'structural Segmentation and normalization (simultaneous Grey/White/CSF segmentation and MNI normalization)',...
    'structural Segmentation and normalization with structural lesion mask (simultaneous Grey/White/CSF segmentation and MNI normalization, creating a new subject-specific TPM with a lesion tissue class)',...
    'functional Indirect segmentation and normalization with structural lesion mask (coregister functional/structural; structural segmentation & normalization; apply same deformation field to functional; create a new subject-specific TPM with a lesion tissue class)',...
    'structural Normalization preserving Grey/White/CSF masks (MNI space normalization of structural, applying same transformation to existing Grey/White/CSF masks)',...
    'structural Manual deformation (non-linear transformation of structural volumes)', ...
    'structural Masking (apply mask to structural data)',...
    'functional Removal of initial scans (disregard initial functional scans)','functional Manual transformation (rotation/flip/translation/affine of functional volumes)','functional Center to (0,0,0) coordinates (translation)','functional Center to structural coordinates (translation)'...
    'functional Slice timing correction (STC; correction for inter-slice differences in acquisition time)','functional Band-pass filtering (temporal filtering of BOLD data)','functional Regression of temporal components (keep residuals of linear model to BOLD timeseries)','functional ROI extraction (compute BOLD timeseres within ROI)','functional Masking (apply mask to functional data)','functional Realignment (subject motion estimation and correction)','functional Realignment with correction of susceptibility distortion interactions (subject motion estimation and correction)',...
    'functional Realignment with susceptibility distortion correction using fieldmaps (subject motion estimation and correction)','functional Outlier detection (ART-based identification of outlier scans for scrubbing)','functional Direct coregistration to structural (rigid body transformation)',...
    'functional Segmentation (Grey/White/CSF segmentation)',...
    'functional Manual deformation (non-linear transformation of functional volumes)',...
    'functional Smoothing (spatial convolution with Gaussian kernel)','functional Motion-mask estimation (BOLD signal derivative wrt movement parameters)',...
    'functional Indirect segmentation and MNI-space normalization (coregister functional/structural; structural segmentation & normalization; apply same deformation field to functional)', ...
    'functional Indirect MNI-space normalization (coregister functional/structural; structural normalization; apply same deformation field to functional)',...
    'functional Indirect MNI-space normalization preserving Grey/White/CSF masks (coregister functional/structural; structural normalization; apply same deformation field to functional and to existing Grey/White/CSF masks)',...
    'functional Direct segmentation and MNI-space normalization (simultaneous Grey/White/CSF segmentation and MNI normalization)',...
    'functional Direct MNI-space normalization (intersubject coregistration)', ...
    'functional Realignment (without reslicing; subject motion estimation and correction)', ...
    'functional Indirect coregistration to structural (non-linear transformation)', ...
    'functional Direct coregistration to structural (without reslicing; rigid body transformation)', ...
    'functional Label current functional files as new secondary dataset (custom label)', ...
    'functional Label current functional files as "original data"', ...
    'functional Label current functional files as "subject-space data"', ...
    'functional Label current functional files as "mni-space data"', ...
    'functional Label current functional files as "surface-space data"', ...
    'functional Label current functional files as "smoothed data"', ...
    'functional Load functional data from previously labeled dataset (custom label)', ...
    'functional Load functional data from "original data" dataset', ...
    'functional Load functional data from "subject-space data" dataset', ...
    'functional Load functional data from "mni-space data" dataset', ...
    'functional Load functional data from "surface-space data" dataset', ...
    'functional Load functional data from "smoothed data" dataset', ...
    'functional Masked smoothing (spatial convolution with Gaussian kernel restricted to voxels within a custom mask)', ...
    'functional Resampling of functional data within the cortical surface (converts volume- to surface- level data using FreeSurfer subject-specific surfaces)', ...
    'functional Smoothing of surface-level functional data (spatial diffusion on surface tessellation)', ...
    'functional Susceptibility distortion correction using voxel-displacement maps (VDM)', ...
    'functional Creation of voxel-displacement maps (VDM) for Susceptibility Distortion Correction' ...
    };
%'functional Normalization (MNI space normalization)','functional Segmentation & Normalization (simultaneous Grey/White/CSF segmentation and MNI normalization)',...
steps_descr={{'INPUT: structural&functional volumes','OUTPUT (all in MNI-space): skull-stripped normalized structural volume, Grey/White/CSF normalized masks, realigned slice-time corrected normalized smoothed functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional&FieldMap volumes',  'OUTPUT (all in MNI-space): skull-stripped normalized structural volume, Grey/White/CSF normalized masks, realigned&unwarp slice-time corrected normalized smoothed functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional&FieldMap volumes',  'OUTPUT (all in MNI-space): skull-stripped normalized structural volume, Grey/White/CSF normalized masks, realigned&unwarp slice-time corrected normalized smoothed functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional volumes','OUTPUT (all in subject-space): skull-stripped structural volume, Grey/White/CSF masks, realigned slice-time corrected coregistered functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional&VDM volumes','OUTPUT (all in subject-space): skull-stripped structural volume, Grey/White/CSF masks, realigned&unwarp slice-time corrected coregistered functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional volumes','OUTPUT (all in subject-space): skull-stripped structural volume, Grey/White/CSF masks, realigned slice-time corrected coregistered functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},...
    {'INPUT: functional data (volume files); structural volume; FreeSurfer-processed structural volume','OUTPUT: functional data (surface files)'}, ...
    {'INPUT: structural volume','OUTPUT: structural volume (same files re-oriented, not resliced)'}, {'INPUT: structural volume','OUTPUT: structural volume (same files translated, not resliced)'}, {'INPUT: structural volume','OUTPUT: skull-stripped structural volume, Grey/White/CSF masks (in same space as structural)'},...
    {'INPUT: structural volume','OUTPUT: skull-stripped normalized structural volume (in MNI space)'},...
    {'INPUT: structural volume','OUTPUT: skull-stripped normalized structural volume, normalized Grey/White/CSF masks (all in MNI space)'},...
    {'INPUT: structural volume & coregistered structural lesion ROI','OUTPUT: skull-stripped normalized structural volume, normalized Grey/White/CSF masks, normalized lesion ROI, [TPM+lesion] maps in "tpm" dataset (all in MNI space)'},...
    {'INPUT: functional volumes; structural volume & coregistered structural lesion ROI','OUTPUT: skull-stripped normalized structural volume, normalized Grey/White/CSF masks, normalized lesion ROI, normalized functional volumes, [TPM+lesion] maps in "tpm" dataset (all in MNI space)'},...
    {'INPUT: structural volume; Grey/White/CSF masks (in same space as structural)','OUTPUT: skull-stripped normalized structural volume, normalized Grey/White/CSF masks (all in MNI space)'},...
    {'INPUT: structural volume; user-defined spatial deformation file (e.g. y_#.nii file)','OUTPUT: resampled structural volumes'}, ...
    {'INPUT: structural volume; ROIs (in same space as structural volume)','OUTPUT: structural volumes masked with ROI (or union of multiple ROIs)'}, ...
    {'INPUT: functional volumes','OUTPUT: temporal subset of functional volumes; temporal subset of first-level covariates (if already defined)'},{'INPUT: functional volumes','OUTPUT: functional volumes (same files re-oriented, not resliced)'},{'INPUT: functional volumes','OUTPUT: functional volumes (same files translated, not resliced)'},{'INPUT: structural and functional volumes','OUTPUT: functional volumes (same files translated, not resliced)'}, ...
    {'INPUT: functional volumes','OUTPUT: slice-timing corrected functional volumes'},{'INPUT: functional volumes','OUTPUT: band-pass filtered functional volumes'},{'INPUT: functional volumes; first-level covariates','OUTPUT: functional volumes with selected covariates regressed-out'},{'INPUT: functional volumes; ROIs (in same space as functional volumes)','OUTPUT: QC_rois first-level covariate with BOLD timeseries within ROIs'},{'INPUT: functional volumes; ROIs (in same space as functional volumes)','OUTPUT: functional volumes masked with ROI (or union of multiple ROIs)'},{'INPUT: functional volumes','OUTPUT: realigned functional volumes, mean functional image, subject movement ''realignment'' 1st-level covariate'},{'INPUT: functional volumes','OUTPUT: realigned&unwarp functional volumes, mean functional image, subject movement ''realignment'' 1st-level covariate'},...
    {'INPUT: functional volumes & VDM maps','OUTPUT: realigned&unwarp functional volumes, mean functional image, subject movement ''realignment'' 1st-level covariate'},{'INPUT: functional volumes, realignment parameters','OUTPUT: outlier scans 1st-level covariate, mean functional image, QA 2nd-level covariates'},{'INPUT: structural and mean functional volume (or first functional)','OUTPUT: coregistered functional volumes'},...
    {'INPUT: mean functional volume (or first functional)','OUTPUT: Grey/White/CSF masks (in same space as functional volume)'},...
    {'INPUT: functional volumes; user-defined spatial deformation file (e.g. y_#.nii file)','OUTPUT: resampled functional volumes'},...
    {'INPUT: functional volumes','OUTPUT: smoothed functional volumes'},{'INPUT: functional volumes','OUTPUT: motion masks'},...
    {'INPUT: structural volume; functional volumes','OUTPUT: skull-stripped normalized structural volume, normalized Grey/White/CSF masks; normalized functional volumes (all in MNI space)'},...
    {'INPUT: structural volume; functional volumes','OUTPUT: skull-stripped normalized structural volume, normalized functional volumes (all in MNI space)'},...
    {'INPUT: structural volume; functional volumes; Grey/White/CSF ROIs (in same space as structural volumes)','OUTPUT: skull-stripped normalized structural volume, normalized functional volumes, normalized Grey/White/CSF masks (all in MNI space)'},...
    {'INPUT: mean functional volume (or first functional)','OUTPUT: normalized functional volumes, normalized Grey/White/CSF masks'},...
    {'INPUT: mean functional volume (or first functional)','OUTPUT: normalized functional volumes'}, ...
    {'INPUT: functional volumes','OUTPUT: realigned functional volumes (same files re-oriented, not resliced), mean functional image, subject movement ''realignment'' 1st-level covariate'}, ...
    {'INPUT: structural and mean functional volume (or first functional)','OUTPUT: functional volumes coregistered to structural (direct normalization to MNI space + inverse deformation field transformation); Grey/White/CSF masks (in same space as functional volume)'}, ...
    {'INPUT: structural and mean functional volume (or first functional)','OUTPUT: functional volumes (all functional volumes are coregistered but not resliced)'}, ...
    {'INPUT: functional volumes','OUTPUT: none (one of the secondary datasets will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("original data" secondary dataset will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("subject-space data" secondary dataset will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("mni-space data" secondary dataset will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("surface-space data" secondary dataset will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("smoothed data" datasets will point to current version of functional volumes)'}, ...
    {'INPUT: other imaging volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes will point to same files as one of the secondary datasets)'}, ...
    {'INPUT: other imaging volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes will point to same files as "original data" secondary dataset)'}, ...
    {'INPUT: other imaging volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes will point to same files as "subject-space data" secondary dataset)'}, ...
    {'INPUT: other imaging volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes will point to same files as "mni-space data" secondary dataset)'}, ...
    {'INPUT: other imaging volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes will point to same files as "surface-space data" secondary dataset)'}, ...
    {'INPUT: other imaging volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes will point to same files as "smoothed data" secondary dataset)'}, ...
    {'INPUT: functional volumes; ROIs (in same space as functional volumes)','OUTPUT: smoothed functional volumes (smoothing restricted within ROI mask -or union of multiple ROIs-)'}, ...
    {'INPUT: functional data (volume files coregistered to structural); FreeSurfer-processed structural volume','OUTPUT: functional data (surface files)'}, ...
    {'INPUT: functional data (surface files)','OUTPUT: smoothed functional data (surface files)'}, ...
    {'INPUT: functional volumes & VDM maps','OUTPUT: Distortion corrected functional volumes'}, ...
    {'INPUT: functional volumes; double-echo FieldMap acquisition files (e.g. Magnitude+PhaseDiff volumes) in "fmap" dataset','OUTPUT: SPM VDM maps in "vdm" dataset'} ...
    };
%{'INPUT: mean functional volume (or first functional)','OUTPUT: normalized functional volumes'},{'INPUT: mean functional volume (or first functional)','OUTPUT: normalized functional volumes, normalized Grey/White/CSF masks '},...
steps_index=num2cell(1:numel(steps));
steps_combinedpipelines={...
    {'functional_label_as_original','functional_realign&unwarp','functional_center','functional_slicetime','functional_art','functional_segment&normalize_direct','functional_label_as_mnispace','structural_center','structural_segment&normalize','functional_smooth','functional_label_as_smoothed'},...
    {'functional_label_as_original','functional_vdm_create','functional_realign&unwarp&fieldmap','functional_center','functional_slicetime','functional_art','structural_center','functional_segment&normalize_indirect','functional_label_as_mnispace','functional_smooth','functional_label_as_smoothed'},...
    {'functional_label_as_original','functional_vdm_create','functional_realign&unwarp&fieldmap','functional_center','functional_slicetime','functional_art','functional_segment&normalize_direct','functional_label_as_mnispace','structural_center','structural_segment&normalize','functional_smooth','functional_label_as_smoothed'},...
    {'functional_label_as_original','functional_realign&unwarp','functional_slicetime','functional_art','functional_coregister_affine_noreslice','functional_label_as_subjectspace','functional_surface_resample','functional_label_as_surfacespace','functional_surface_smooth','functional_label_as_smoothed','structural_segment'},...
    {'functional_label_as_original','functional_vdm_create','functional_realign&unwarp&fieldmap','functional_slicetime','functional_art','functional_coregister_affine_noreslice','functional_label_as_subjectspace','functional_surface_resample','functional_label_as_surfacespace','functional_surface_smooth','functional_label_as_smoothed','structural_segment'},...
    {'functional_label_as_original','functional_realign&unwarp','functional_slicetime','functional_art','functional_coregister_nonlinear','functional_label_as_subjectspace','functional_surface_resample','functional_label_as_surfacespace','functional_surface_smooth','functional_label_as_smoothed'},...
    {'functional_coregister_affine_noreslice','functional_surface_resample'} ...
    };
for n=1:numel(steps_combinedpipelines),
    [ok,idx]=ismember(steps_combinedpipelines{n},steps);
    if ~all(ok), error('preprocessing step names have changed'); end
    steps_index{n}=idx(:)';
end
if nargin>0&&ischar(STEPS)&&strcmp(STEPS,'steps'), [ok,idx]=sort(steps); matlabbatch=steps_names(idx); return; end
steps_pipelines=cellfun('length',steps_index)>1;
steps_default=cellfun('length',regexp(steps,'^default_'))>0;
dogui=false;
sets=0;
subjects=1:CONN_x.Setup.nsubjects;
sessions=1:max(CONN_x.Setup.nsessions);
doimport=true;
typeselect='';
multiplesteps=iscell(STEPS)&numel(STEPS)>1;
showallsteps=false;
voxelsize_anat=1;
voxelsize_func=2;
boundingbox=[-90,-126,-72;90,90,108]; % default bounding-box
interp=[];
fwhm=[];
diffusionsteps=[];
sliceorder=[];
sliceorder_select=[];
label={};
load_label={};
ta=[];
unwarp=[];
removescans=[];
bp_filter=[];
bp_keep0=1;
reg_names={};
reg_dimensions=[];
reg_deriv=[];
reg_filter=[];
reg_detrend=1;
reg_lag=[];
reg_lagmax=8;
reg_skip=0;
roi_names={};
roi_dimensions=[];
roi_deriv=[];
roi_filter=[];
roi_detrend=0;
roi_scale=1;
mask_names_func={};
mask_inclusive_func=1;
mask_names_anat={};
mask_inclusive_anat=1;
reorient=[];
respatialdef=[];
coregtomean=1;
rtm=0;
rmask=1;
coregsource={};
applytofunctional=false;
tpm_template=[];
tpm_structlesion=[];
tpm_overwrite=false;
tpm_lesionscale=1;
tpm_lesionprobabilistic=false;
tpm_ngaus=[]; % e.g. [1 1 2 3 4 2]
affreg='mni';
warpreg=[]; % e.g. [0 0.001 0.5 0.05 0.2];
vdm_et1=[]; % eg. 2.84, 4.37;
vdm_et2=[]; % eg. 5.30, 6.83
vdm_ert=[]; % eg. 37.6
vdm_blip=[];% eg. -1
vdm_type=[];
vdm_fmap=[];
art_thresholds=[];
art_useconservative=1;
art_global_thresholds=[9 5 3];
art_motion_thresholds=[2 .9 .5];
art_global_threshold=art_global_thresholds(1+art_useconservative); % default art scan-to-scan global signal z-value thresholds
art_motion_threshold=art_motion_thresholds(1+art_useconservative); % default art scan-to-scan composite motion mm thresholds
art_use_diff_motion=1;
art_use_diff_global=1;
art_use_norms=1;
art_force_interactive=0;
art_drop_flag=0;
art_gui_display=true;
parallel_profile=[];
parallel_N=0;
functional_template=fullfile(fileparts(which('spm')),'templates','EPI.nii');
if ~conn_existfile(functional_template), functional_template=fullfile(fileparts(which('spm')),'toolbox','OldNorm','EPI.nii'); end
structural_template=fullfile(fileparts(which('spm')),'templates','T1.nii');
if ~conn_existfile(structural_template), structural_template=fullfile(fileparts(which('spm')),'toolbox','OldNorm','T1.nii'); end
selectedstep=1; if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), subjects=CONN_x.pobj.subjects; end % this field overwrites user-defined options
if ~isempty(sets)&&ischar(sets), sets=conn_datasetlabel(sets,'error'); end

if ~nargin||isempty(STEPS)||dogui,
    dogui=true;
    if showallsteps, idx=1:numel(steps);
    else idx=find(cellfun('length',regexp(steps,'^structural|^functional|^default')));
    end
    if ~isempty(typeselect)
        switch(typeselect)
            case 'structural', idx=find(cellfun('length',regexp(steps,'^structural'))&~steps_pipelines);
            case 'functional', idx=find(cellfun('length',regexp(steps,'^functional'))&~steps_pipelines);
        end
    end
    steps0=steps;
    steps=steps(idx);
    steps_names=steps_names(idx);
    steps_descr=steps_descr(idx);
    steps_pipelines=steps_pipelines(idx);
    steps_default=steps_default(idx);
    steps_index=steps_index(idx);
    [nill,steps_order]=sort(steps_names);
    steps_order=[sort(steps_order(steps_default(steps_order))) steps_order(~steps_default(steps_order))];
    scalefig=1+multiplesteps;
    dlg.steps=steps;
    dlg.steps_names=steps_names;
    dlg.steps_descr=steps_descr;
    dlg.steps_index=steps_index;
    dlg.steps_order=steps_order;
    dlg.fig=figure('units','norm','position',[.2,.3,.5+.2*(1|multiplesteps),.6],'menubar','none','numbertitle','off','name','SPM data preprocessing step','color',1*[1 1 1]);
    if multiplesteps,
        uicontrol('style','frame','units','norm','position',[0,.57,1,.43],'backgroundcolor',.9*[1 1 1],'foregroundcolor',.9*[1 1 1]);
        uicontrol('style','frame','units','norm','position',[.05,.2,.9,.33],'backgroundcolor',1*[1 1 1],'foregroundcolor',.8*[1 1 1]);
        %uicontrol('style','frame','units','norm','position',[.025,.025,.95,.55],'backgroundcolor',1*[1 1 1],'foregroundcolor',.75*[1 1 1],'fontsize',9+font_offset);
    end
    htm0=uicontrol('style','text','units','norm','position',[.05,.9,.85,.05],'backgroundcolor',1*[1 1 1],'foregroundcolor','k','horizontalalignment','left','string','Select individual preprocessing step:','fontweight','bold','fontsize',9+font_offset);
    dlg.m0=uicontrol('style','popupmenu','units','norm','position',[.05,.85,.85,.05],'string',steps_names(steps_order),'value',find(ismember(steps_order,selectedstep)),'backgroundcolor',1*[1 1 1],'foregroundcolor','k','tooltipstring','Select a data preprocessing step','callback',@(varargin)conn_setup_preproc_update,'fontsize',9+font_offset);
    if multiplesteps,
        set(htm0,'string','List of all available preprocessing steps:');
        set(dlg.m0,'tooltipstring','Select a data preprocessing step or pipeline and click ''Add'' to add it to your data preprocessing pipeline');
    end
    dlg.m6=uicontrol('style','text','units','norm','position',[.05,.725,.85,.1],'max',2,'string','','backgroundcolor',1*[1 1 1],'enable','inactive','horizontalalignment','left','fontsize',9+font_offset);
    dlg.m4=uicontrol('style','checkbox','units','norm','position',[.05,.65,.85,.05],'value',~coregtomean,'string','First functional volume as reference','backgroundcolor',1*[1 1 1],'tooltipstring','<HTML>Uses firts functional volume as reference in coregistration/normalization step <br/> - if unchecked coregistration/normalization uses mean-volume as reference instead<br/> - note: mean volume is created during realignment</HTML>','visible','off','fontsize',9+font_offset);
    %dlg.m3=uicontrol('style','checkbox','units','norm','position',[.1,.5,.8/scalefig,.05],'value',applytofunctional,'string','Apply structural deformation field to functional data as well','backgroundcolor',1*[1 1 1],'tooltipstring','Apply structural deformation field computed during structural normalization/segmentation step to coregistered functional data as well','visible','off','fontsize',9+font_offset);
    dlg.m2=uicontrol('style','popupmenu','units','norm','position',[.05,.55,.85,.05],'value',1,'string',{'Run process and import results to CONN project','Run process only (do not import results)','Interactive SPM batch editor only (do not run process)'}','backgroundcolor',1*[1 1 1],'fontsize',9+font_offset);
    dlg.m1=uicontrol('style','checkbox','units','norm','position',[.05,.08,.3,.05],'value',1,'string','Process all subjects','backgroundcolor',1*[1 1 1],'tooltipstring','Apply this preprocessing to all subjects in your curent CONN project','callback',@(varargin)conn_setup_preproc_update,'fontsize',9+font_offset);
    dlg.m5=uicontrol('style','listbox','units','norm','position',[.35,.11,.15,.08],'max',2,'string',arrayfun(@(n)sprintf('Subject%d',n),1:CONN_x.Setup.nsubjects,'uni',0),'backgroundcolor',1*[1 1 1],'tooltipstring','Select subjects','visible','off','fontsize',9+font_offset);
    dlg.m1b=uicontrol('style','checkbox','units','norm','position',[.05,.03,.3,.05],'value',1,'string','Process all sessions','backgroundcolor',1*[1 1 1],'tooltipstring','Apply this preprocessing to all sessions in your curent CONN project','callback',@(varargin)conn_setup_preproc_update,'fontsize',9+font_offset);
    dlg.m5b=uicontrol('style','listbox','units','norm','position',[.35,.03,.15,.08],'max',2,'string',arrayfun(@(n)sprintf('Session%d',n),1:max(CONN_x.Setup.nsessions),'uni',0),'backgroundcolor',1*[1 1 1],'tooltipstring','Select sessions','visible','off','fontsize',9+font_offset);
    if any(cellfun('length',regexp(dlg.steps_names,'^functional'))), dlg.m10=uicontrol('style','popupmenu','units','norm','position',[.05,.13,.3,.05],'value',1+sets,'string',[{'Process functional/structural data'},arrayfun(@(n)sprintf('Process secondary dataset #%d %s',n,regexprep(CONN_x.Setup.secondarydataset(n).label,'(.+)','($1)')),1:numel(CONN_x.Setup.secondarydataset),'uni',0)],'tooltipstring','<HTML>Apply this preprocessing to selected functional dataset in your curent CONN project (the same dataset will hold the OUTPUT files of each preprocessing step)<br/> - note: only when preprocessing the functional data (dataset-0) non-essential OUTPUTS of each preprocessing step (ROIs, first- and <br/>second- level covariates) will be automatically imported into your CONN project</HTML>','visible','on','fontsize',9+font_offset);
    else dlg.m10=[];
    end
    [tstr,tidx]=conn_jobmanager('profiles');
    tnull=find(strcmp('Null profile',conn_jobmanager('profiles')));
    tlocal=find(strcmp('Background process (Unix,Mac)',tstr),1);
    tvalid=setdiff(1:numel(tstr),tnull);
    tstr=cellfun(@(x)sprintf('distributed processing (run on %s)',x),tstr,'uni',0);
    if 1, tvalid=tidx; if isunix&&~isempty(tlocal)&&~ismember(tlocal,tvalid), tvalid=[tvalid(:)' tlocal]; end
    elseif 1, tvalid=tidx; % show only default scheduler
    else tstr{tidx}=sprintf('<HTML><b>%s</b></HTML>',tstr{tidx});
    end
    toptions=[{'local processing (run on this computer)' 'queue/script it (save as scripts to be run later)'} tstr(tvalid)];
    if CONN_gui.isremote
        info=conn_remotely('info');
        if isfield(info,'host')&&~isempty(info.host), tnameserver=info.host;
        elseif isfield(info,'remote_ip')&&~isempty(info.remote_ip), tnameserver=info.remote_ip;
        else tnameserver='CONN server';
        end
        toptions=regexprep(toptions,'\<run on (this computer)?',['run on ',tnameserver,' ']);
    end    
    dlg.m9=uicontrol('style','popupmenu','units','norm','position',[.55,.12,.40,.05],'string',toptions,'value',1,'backgroundcolor',1*[1 1 1],'fontsize',8+CONN_gui.font_offset);
    if multiplesteps, dlg.m11=uicontrol('style','pushbutton','units','norm','position',[.55,.04,.2,.07],'string','Start','tooltipstring','Accept changes and run data preprocessing pipeline','callback','set(gcbf,''userdata'',0); uiresume(gcbf)','fontsize',9+font_offset,'fontweight','bold');
    else              dlg.m11=uicontrol('style','pushbutton','units','norm','position',[.55,.04,.2,.07],'string','Start','tooltipstring','Accept changes and run data preprocessing step','callback','set(gcbf,''userdata'',0); uiresume(gcbf)','fontsize',9+font_offset);
    end
    dlg.m12=uicontrol('style','pushbutton','units','norm','position',[.75,.04,.2,.07],'string','Cancel','callback','delete(gcbf)','fontsize',9+font_offset);
    if multiplesteps
        set([htm0 dlg.m0],'visible','off');
        dlg.m0b=uicontrol('style','text','units','norm','position',[.07,.45,.86,.06],'backgroundcolor',1*[1 1 1],'foregroundcolor',0*[1 1 1],'horizontalalignment','left','string','','fontweight','bold','fontsize',9+font_offset);
        set(dlg.m6,'position',[.07,.26,.86,.18],'backgroundcolor',1*[1 1 1]);
        set(dlg.m4,'position',[.07,.215,.86,.04],'backgroundcolor',1*[1 1 1]);
        set(dlg.m2,'visible','off');%'string',{'Run process and import results to CONN project'});
        %set(dlg.m3,'position',get(dlg.m3,'position')-[0 .075 0 0]);
        %set(dlg.m4,'position',get(dlg.m4,'position')-[0 .075 0 0]);
        set(dlg.fig,'name','SPM data preprocessing pipeline');
        uicontrol('style','text','units','norm','position',[.05,.915,.85,.05],'backgroundcolor',.9*[1 1 1],'foregroundcolor','k','horizontalalignment','left','string','Data preprocessing pipeline:','fontweight','bold','fontsize',11+font_offset);
        dlg.m7=uicontrol('style','listbox','units','norm','position',[.05,.59,.78,.325],'max',2,'string',{},'backgroundcolor',.9*[1 1 1],'tooltipstring','Define sequence of preprocessing steps','fontsize',9+font_offset,'callback','dlg=get(gcbo,''userdata''); str=get(gcbo,''string''); val=get(gcbo,''value''); if numel(val)==1, idx=find(strcmp(dlg.steps_names(dlg.steps_order),str{val})); if numel(idx)==1, set(dlg.m0,''value'',idx); feval(get(dlg.m0,''callback'')); end; end');
        dlg.m8a=uicontrol('style','pushbutton','units','norm','position',[.84,.87,.11,.045],'string','Add','fontweight','bold','tooltipstring','Adds new data preprocessing step to this list','callback',@conn_setup_preproc_update_add,'fontsize',9+font_offset);
        dlg.m8b=uicontrol('style','pushbutton','units','norm','position',[.84,.825,.11,.045],'string','Remove','tooltipstring','Removes selected preprocessing step from this list','callback','dlg=get(gcbo,''userdata''); str=get(dlg.m7,''string''); str=str(setdiff(1:numel(str),get(dlg.m7,''value''))); set(dlg.m7,''string'',str,''value'',[]); feval(get(dlg.m0,''callback'')); ','fontsize',9+font_offset);
        dlg.m8g=uicontrol('style','pushbutton','units','norm','position',[.84,.78,.11,.045],'string','Clear','tooltipstring','Removes all preprocessing steps from this list','callback','dlg=get(gcbo,''userdata''); set(dlg.m7,''string'',{},''value'',[]); feval(get(dlg.m0,''callback'')); ','fontsize',9+font_offset);
        dlg.m8c=uicontrol('style','pushbutton','units','norm','position',[.84,.735,.11,.045],'string','Move up','tooltipstring','Moves selected preprocessing step up in this list','callback','dlg=get(gcbo,''userdata''); str=get(dlg.m7,''string''); val=get(dlg.m7,''value''); idx=1:numel(str); idx(val)=min(idx(val))-1.5; [nill,idx]=sort(idx); str=str(idx); set(dlg.m7,''string'',str,''value'',find(rem(nill,1)~=0));','fontsize',9+font_offset);
        dlg.m8d=uicontrol('style','pushbutton','units','norm','position',[.84,.69,.11,.045],'string','Move down','tooltipstring','Moves selected preprocessing step down this list','callback','dlg=get(gcbo,''userdata''); str=get(dlg.m7,''string''); val=get(dlg.m7,''value''); idx=1:numel(str); idx(val)=max(idx(val))+1.5; [nill,idx]=sort(idx); str=str(idx); set(dlg.m7,''string'',str,''value'',find(rem(nill,1)~=0));','fontsize',9+font_offset);
        dlg.m8e=uicontrol('style','pushbutton','units','norm','position',[.84,.635,.11,.045],'string','Save','tooltipstring','Saves this data preprocessing pipeline list for future use','callback',@conn_setup_preproc_save,'fontsize',9+font_offset);
        dlg.m8f=uicontrol('style','pushbutton','units','norm','position',[.84,.59,.11,.045],'string','Load','tooltipstring','Loads data preprocessing pipeline list from file','callback',@conn_setup_preproc_load,'fontsize',9+font_offset);
        set([dlg.m7 dlg.m8a dlg.m8b dlg.m8c dlg.m8d dlg.m8e dlg.m8f dlg.m8g],'userdata',dlg);
    else dlg.m7=[]; dlg.m0b=[];
    end
    set([dlg.m0 dlg.m1 dlg.m1b],'userdata',dlg);
    %if isempty(STEPS)&&multiplesteps, STEPS=steps(steps_index{1}); end
    if ~isempty(STEPS)
        [tok,idx]=ismember(STEPS,steps);
        if multiplesteps, set(dlg.m7,'string',steps_names(idx(tok>0))');
        else set(dlg.m0, 'value',find(ismember(steps_order,idx(tok>0)),1));
        end
    end
    conn_setup_preproc_update(dlg.m0);
    if multiplesteps,
        conn_setup_preproc_load(dlg.m8f);
        if isempty(get(dlg.m7,'string')), conn_setup_preproc_update_add(dlg.m8a); end
    end
    uiwait(dlg.fig);
   
    if ~ishandle(dlg.fig), return; end
    pressedok=get(dlg.fig,'userdata');
    if isempty(pressedok), return; end
    if multiplesteps
        STEPS=get(dlg.m7,'string');
        [tok,idx]=ismember(STEPS,steps_names);
        STEPS=steps(idx(tok>0));
    else
        STEPS=steps(dlg.steps_order(get(dlg.m0,'value')));
        %idx0=find(steps_pipelines);
        %[ok,idx]=ismember(lower(STEPS),steps(idx0));
        %if ok, STEPS=steps0(steps_index{idx0(idx)}); end
    end
    %STEP_name=steps_names{get(dlg.m0,'value')};
    %if any(ismember(STEPS,{'structural_segment&normalize','structural_normalize'})), applytofunctional=get(dlg.m3,'value'); end
    if any(cellfun('length',regexp(STEPS,'^functional_coregister|^functional_normalize|functional_segment|functional_center'))), coregtomean=~get(dlg.m4,'value'); end
    if ~get(dlg.m1,'value'), subjects=get(dlg.m5,'value'); end
    if ~get(dlg.m1b,'value'), sessions=get(dlg.m5b,'value'); end
    if ~isempty(dlg.m10), sets=get(dlg.m10,'value')-1; end
    dorun=get(dlg.m2,'value');
    doparallel=get(dlg.m9,'value');
    if multiplesteps, conn_setup_preproc_save(dlg.m8f); end
    delete(dlg.fig);
    switch(dorun)
        case 1, STEPS=cellfun(@(x)['run_',x],STEPS,'uni',0); doimport=true;
        case 2, STEPS=cellfun(@(x)['run_',x],STEPS,'uni',0); doimport=false;
        case 3, STEPS=cellfun(@(x)['interactive_',x],STEPS,'uni',0); doimport=false;
        case 4, STEPS=cellfun(@(x)['update_',x],STEPS,'uni',0); doimport=true;
    end
   
    if doparallel>1
        if doparallel==2, parallel_profile=find(strcmp('Null profile',conn_jobmanager('profiles')));
        else parallel_profile=tvalid(doparallel-2);
            if conn_jobmanager('ispending')
                answ=conn_questdlg({'There are previous pending jobs associated with this project','This job cannot be submitted until all pending jobs finish',' ','Would you like to queue this job for later?','(pending jobs can be seen at Tools.Cluster/HPC.View pending jobs'},'Warning','Queue','Cancel','Queue');
                if isempty(answ)||strcmp(answ,'Cancel'), ok=false; end
                parallel_profile=find(strcmp('Null profile',conn_jobmanager('profiles')));
            end
        end
        if numel(subjects)>1,
            answer=conn_menu_inputdlg(sprintf('Number of parallel jobs? (1-%d)',numel(subjects)),'CONN HPC',1,{num2str(1)});
            if isempty(answer)||isempty(str2num(answer{1})), return; end
            parallel_N=str2num(answer{1});
        else parallel_N=1;
        end
    end
end

lSTEPS=regexprep(lower(STEPS),'^run_|^update_|^interactive_','');
sliceorder_select_options={'ascending','descending','interleaved (middle-top)','interleaved (bottom-up)','interleaved (top-down)','interleaved (Siemens)','interleaved (Philips)','BIDS'};
sliceorder_select_options_extended={'ascending (e.g. 1,2,3...9,10)','descending (e.g. 10,9,8...2,1)','interleaved middle-top (e.g. 10,5,9,4...6 1)','interleaved bottom-up (e.g. 1,3,5...8,10)','interleaved top-down (e.g. 10,8,6...3,1)','interleaved (Siemens) (e.g. 2,4,6...7,9)','interleaved (Philips) (e.g. 1,4,7...6,9)','BIDS (from functional .json metadata)'};
if any(ismember('functional_slicetime',lSTEPS))
    if ischar(sliceorder),
        [slok,sliceorder_select]=ismember(sliceorder,sliceorder_select_options);
        if ~slok, conn_disp(sprintf('Warning: incorrect sliceorder name %s',sliceorder)); sliceorder_select=[]; end
        sliceorder=[];
    end
    if isempty(sliceorder)&&(isempty(sliceorder_select)||dogui)
        [sliceorder_select,tok] = listdlg('PromptString','Select slice order:','ListSize',[400 200],'SelectionMode','single','InitialValue',sliceorder_select,'ListString',[sliceorder_select_options_extended,{'manually define','do not know (skip slice timing correction)'}]);
        if isempty(sliceorder_select), return; end
        if sliceorder_select==numel(sliceorder_select_options)+1
            sliceorder=inputdlg(['Slice order? (enter slice indexes from z=1 -first slice in image- to z=? -last slice- in the order they were acquired). Alternatively enter acquisition time of each slice in milliseconds (e.g. for multiband sequences). Press Cancel to enter this information at a later point (e.g. separately for each subject)'],'conn_setup_preproc',1,{' '});
            if ~isempty(sliceorder), sliceorder=str2num(sliceorder{1});
            else sliceorder=[];
            end
        elseif sliceorder_select==numel(sliceorder_select_options)+2
            sliceorder_select=[]; STEPS=STEPS(~ismember(lSTEPS,'functional_slicetime'));
        end
    end
end

if any(ismember('functional_removescans',lSTEPS))
    if isempty(removescans)||dogui
        if isempty(removescans), removescans=0; end
        removescans=conn_menu_inputdlg('Enter number of initial scans to remove','conn_setup_preproc',1,{num2str(removescans)});
        if isempty(removescans), return; end
        removescans=str2num(removescans{1});
    end
end

if any(ismember('functional_bandpass',lSTEPS))
    if isempty(bp_filter)||dogui
        if isempty(bp_filter), bp_filter=[0.01 0.10]; end
        bp_filter=conn_menu_inputdlg('Enter band-pass filter thresholds in Hz','conn_setup_preproc',1,{num2str(bp_filter)});
        if isempty(bp_filter), return; end
        bp_filter=str2num(bp_filter{1});
    end
end

if any(ismember('functional_regression',lSTEPS))
    if isempty(reg_names)||dogui
        temp_reg_names=CONN_x.Setup.l1covariates.names(1:end-1);
        if isempty(reg_names), answ=~cellfun('length',regexp(temp_reg_names,'^QC_'));
        else answ=reshape(ismember(temp_reg_names,reg_names),1,[]);
        end
        if isempty(reg_deriv), tansw=~cellfun('length',regexpi(temp_reg_names,'^effect of|realign|movement|motion'));
        else tansw=reshape(reg_deriv==0,1,[]);
        end
        answ=reshape([answ&tansw;answ&~tansw],1,[]);
        temp_reg_names=reshape([temp_reg_names;cellfun(@(x)sprintf('%s + 1st order temporal derivative',x),temp_reg_names,'uni',0)],1,[]);
        temp_reg_names=[{'time (detrending)'},temp_reg_names];
        if isempty(reg_detrend), answ=[true, answ];
        else answ=[reg_detrend, answ];
        end
        answ=listdlg('Promptstring','Select model regressors','selectionmode','multiple','liststring',temp_reg_names,'initialvalue',find(answ),'ListSize',[320 300]);
        if numel(answ)>=1,
            reg_detrend=any(answ==1);
            answ=answ(answ>1);
            uansw=unique(2*floor(answ/2));
            reg_names=temp_reg_names(uansw);
            reg_deriv=double(ismember(uansw+1,answ));
        else return;
        end
    end
end

if any(ismember('functional_roiextract',lSTEPS))
    if isempty(roi_names)||dogui
        temp_roi_names=reshape(CONN_x.Setup.rois.names(1:end-1),1,[]);
        if isempty(roi_names),
            answ=false(size(temp_roi_names));
        else
            temp_roi_names=[reshape(roi_names,1,[]), reshape(temp_roi_names(~ismember(temp_roi_names,roi_names)),1,[])];
            answ=ismember(temp_roi_names,roi_names);
        end
        answ=listdlg('Promptstring','Select ROIs','selectionmode','multiple','liststring',temp_roi_names,'initialvalue',find(answ),'ListSize',[320 300]);
        if numel(answ)>=1,
            roi_names=temp_roi_names(answ);
        else return;
        end
    end
end

if any(ismember({'functional_mask','functional_smooth_masked'},lSTEPS))
    if isempty(mask_names_func)||dogui
        temp_mask_names_func=reshape(CONN_x.Setup.rois.names(1:end-1),1,[]);
        temp_mask_names_func0=reshape([temp_mask_names_func; temp_mask_names_func],1,[]);
        temp_mask_names_func1=reshape([cellfun(@(x)[x ' (inclusive mask; keep only voxels within this ROI)'],temp_mask_names_func,'uni',0); cellfun(@(x)[x ' (exclusive mask; exclude only voxels within this ROI)'], temp_mask_names_func,'uni',0)],1,[]);
        if isempty(mask_names_func),
            answ=false(size(temp_mask_names_func0));
        else
            answ=ismember(temp_mask_names_func0,mask_names_func);
        end
        answ=listdlg('Promptstring','Select ROI(s) for functional masking','selectionmode','single','liststring',temp_mask_names_func1,'initialvalue',find(answ),'ListSize',[320 300]);
        if numel(answ)>=1,
            mask_names_func=temp_mask_names_func0(answ);
            mask_inclusive_func=rem(answ,2);
        else return;
        end
    end
end

if any(ismember('structural_mask',lSTEPS))
    if isempty(mask_names_anat)||dogui
        temp_mask_names_anat=reshape(CONN_x.Setup.rois.names(1:end-1),1,[]);
        temp_mask_names_anat0=reshape([temp_mask_names_anat; temp_mask_names_anat],1,[]);
        temp_mask_names_anat1=reshape([cellfun(@(x)[x ' (inclusive mask; keep only voxels within this ROI)'],temp_mask_names_anat,'uni',0); cellfun(@(x)[x ' (exclusive mask; exclude only voxels within this ROI)'], temp_mask_names_anat,'uni',0)],1,[]);
        if isempty(mask_names_anat),
            answ=false(size(temp_mask_names_anat0));
        else
            answ=ismember(temp_mask_names_anat0,mask_names_anat);
        end
        answ=listdlg('Promptstring','Select ROI(s) for structural masking','selectionmode','single','liststring',temp_mask_names_anat1,'initialvalue',find(answ),'ListSize',[320 300]);
        if numel(answ)>=1,
            mask_names_anat=temp_mask_names_anat0(answ);
            mask_inclusive_anat=rem(answ,2);
        else return;
        end
    end
end

if dogui&&any(ismember(lSTEPS,{'functional_vdm_create'}))
    thfig=figure('units','norm','position',[.4,.4,.35,.3],'color',1*[1 1 1],'name','VDM create settings','numbertitle','off','menubar','none');
    ht1=uicontrol('style','popupmenu','units','norm','position',[.1,.85,.8,.1],'string',arrayfun(@(n)sprintf('Fieldmap location: secondary dataset #%d %s',n,regexprep(CONN_x.Setup.secondarydataset(n).label,'(.+)','($1)')),1:numel(CONN_x.Setup.secondarydataset),'uni',0),'value',1,'backgroundcolor',1*[1 1 1],'tooltipstring','defines location of available fieldmap-sequence files');
    ht2=uicontrol('style','popupmenu','units','norm','position',[.1,.75,.8,.1],'string',{'Fieldmap type: automatically determine','Fieldmap type: Magnitude,Phasediff files','Fieldmap type: Real1,Imag1,Real2,Imag2 files','Fieldmap type: Pre-computed fieldmap file (Hz)'},'value',1,'backgroundcolor',1*[1 1 1],'tooltipstring','defines type of available fieldmap-sequence files');
    ht3=uicontrol('style','checkbox','units','norm','position',[.1,.64,.8,.1],'string','Read double-echo timing from BIDS / .json files','value',1,'backgroundcolor',1*[1 1 1],'tooltipstring','use information in .json sidecar files to estimate EchoTime and EPI Total Readout Time values');
    ht4=[];ht5=[];ht6=[];
    ht4a=uicontrol('style','text','units','norm','position',[.1,.5,.6,.1],'string','Fieldmap''s Short Echo Time (in ms)','horizontalalignment','left','backgroundcolor',1*[1 1 1],'enable','off');
    ht4=uicontrol('style','edit','units','norm','position',[.7,.5,.2,.1],'string',num2str(vdm_et1),'tooltipstring','defines Echo Time (in ms units) of first dual-echo acquisition (leave empty to import from .json / BIDS file)','enable','off');
    ht5a=uicontrol('style','text','units','norm','position',[.1,.4,.6,.1],'string','Fieldmap''s Long Echo Time (in ms)','horizontalalignment','left','backgroundcolor',1*[1 1 1],'enable','off');
    ht5=uicontrol('style','edit','units','norm','position',[.7,.4,.2,.1],'string',num2str(vdm_et2),'tooltipstring','defines Echo Time (in ms units) of second dual-echo acquisition (leave empty to import from .json / BIDS file)','enable','off');
    ht6a=uicontrol('style','text','units','norm','position',[.1,.3,.6,.1],'string','Functional''s Total Readout Time (in ms)','horizontalalignment','left','backgroundcolor',1*[1 1 1],'enable','off');
    ht6=uicontrol('style','edit','units','norm','position',[.7,.3,.2,.1],'string',num2str(vdm_ert),'tooltipstring','defines EPI Total Readout Time (in ms units) of functional data (note: equal to 1000/BandwidthPerPixelPhaseEncode;  leave empty to import from .json / BIDS file)','enable','off');
    ht7a=uicontrol('style','text','units','norm','position',[.1,.2,.6,.1],'string','Functional''s Blip direction (+1,-1,S,R)','horizontalalignment','left','backgroundcolor',1*[1 1 1],'enable','off');
    if isempty(vdm_blip), tvdm_blip=-1; else tvdm_blip=vdm_blip; end
    if ~ischar(tvdm_blip), tvdm_blip=num2str(tvdm_blip); end
    ht7=uicontrol('style','edit','units','norm','position',[.7,.2,.2,.1],'string',tvdm_blip,'tooltipstring','defines k-space traversal blip direction: +1 for positive direction, -1 for negative direction, leave empty or set to ''S'' to derive this information from the PhaseEncodingDirection field in .json/BIDS file, set to ''R'' to try the reverse direction of ''S''','enable','off');
    uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.15],'callback','uiresume');
    uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.15],'callback','delete(gcbf)');
    onoff={'on','off'};
    if isempty(vdm_fmap), vdm_fmap=conn_datasetlabel('fmap'); end
    if ischar(vdm_fmap), vdm_fmap=conn_datasetlabel(vdm_fmap); end
    if isempty(vdm_fmap), vdm_fmap=1; end
    set(ht1,'value',vdm_fmap);
    if isempty(vdm_type), set(ht2,'value',1);
    else set(ht2,'value',1+vdm_type); set([ht4 ht4a ht5 ht5a],'visible',onoff{1+(get(ht2,'value')==4)});
    end
    if isempty(vdm_et1)&&isempty(vdm_et2)&&isempty(vdm_ert)&&isempty(vdm_blip), set(ht3,'value',1);
    else set(ht3,'value',0); set([ht4 ht4a ht5 ht5a ht6 ht6a ht6 ht6a ht7 ht7a],'enable','on');
    end
    set(ht2,'userdata',[],'callback',@(varargin)set([ht4 ht4a ht5 ht5a],'visible',onoff{1+(get(ht2,'value')==4)}));
    set(ht3,'userdata',[],'callback',@(varargin)set([ht4 ht4a ht5 ht5a ht6 ht6a ht7 ht7a],'enable',onoff{1+get(ht3,'value')}));
    uiwait(thfig);
    if ~ishandle(thfig), return; end
    vdm_fmap=get(ht1,'value');
    vdm_type=get(ht2,'value')-1; if ~vdm_type, vdm_type=[]; end
    if get(ht3,'value'), vdm_et1=[]; vdm_et2=[]; vdm_ert=[]; vdm_blip=[];
    else
        temp=get(ht4,'string'); if isempty(temp), vdm_et1=temp; elseif ~isempty(str2num(temp)), vdm_et1=str2num(temp); else error('unable to interpret vdm_et1 input %s',temp); end
        temp=get(ht5,'string'); if isempty(temp), vdm_et2=temp; elseif ~isempty(str2num(temp)), vdm_et2=str2num(temp); else error('unable to interpret vdm_et2 input %s',temp); end
        temp=get(ht6,'string'); if isempty(temp), vdm_ert=temp; elseif ~isempty(str2num(temp)), vdm_ert=str2num(temp); else error('unable to interpret vdm_ert input %s',temp); end
        temp=get(ht7,'string'); if isempty(temp)||isequal(temp,'s')||isequal(temp,'S'), vdm_blip=[]; elseif ~isempty(str2num(temp)), vdm_blip=str2num(temp); elseif isequal(temp,'r')||isequal(temp,'R'), vdm_blip=0; else error('unable to interpret vdm_blip input %s',temp); end
    end
    delete(thfig);
    drawnow;
end

if any(ismember({'structural_manualorient','functional_manualorient'},lSTEPS))
    if isempty(reorient)||dogui
        ntimes=sum(ismember(lSTEPS,{'structural_manualorient','functional_manualorient'}));
        reorient={};
        opts={'translation to 0/0/0 coordinates',nan;
            '90-degree rotation around x-axis (x/y/z to x/-z/y)',[1 0 0;0 0 1;0 -1 0];
            '90-degree rotation around x-axis (x/y/z to x/z/-y)',[1 0 0;0 0 -1;0 1 0];
            '90-degree rotation around y-axis (x/y/z to -z/y/x)',[0 0 1;0 1 0;-1 0 0];
            '90-degree rotation around y-axis (x/y/z to z/y/-x)',[0 0 -1;0 1 0;1 0 0];
            '90-degree rotation around z-axis (x/y/z to y/-x/z)',[0 -1 0;1 0 0;0 0 1];
            '90-degree rotation around z-axis (x/y/z to -y/x/z)',[0 1 0;-1 0 0;0 0 1];
            '180-degree rotation around x-axis (x/y/z to x/-y/-z)',[1 0 0;0 -1 0;0 0 -1];
            '180-degree rotation around y-axis (x/y/z to -x/y/-z)',[-1 0 0;0 1 0;0 0 -1];
            '180-degree rotation around z-axis (x/y/z to -x/-y/z)',[-1 0 0;0 -1 0;0 0 1];
            'clockwise rotation around x-axis (arbitrary angle)',@(a)[1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)];
            'clockwise rotation around y-axis (arbitrary angle)',@(a)[cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
            'clockwise rotation around z-axis (arbitrary angle)',@(a)[cos(a) sin(a) 0;-sin(a) cos(a) 0;0 0 1];
            'non-rigid reflection along x-axis (x/y/z/ to -x/y/z)', [-1 0 0;0 1 0;0 0 1];
            'non-rigid reflection along y-axis (x/y/z/ to x/-y/z)', [1 0 0;0 -1 0;0 0 1];
            'non-rigid reflection along z-axis (x/y/z/ to x/y/-z)', [1 0 0;0 1 0;0 0 -1];
            'arbitrary affine transformation matrix (manually define 4x4 matrix)', 1;
            'arbitrary affine transformation matrix (load 4x4 matrix from file)', 2};
        for ntime=1:ntimes
            if ntimes>1 [treorient,tok] = listdlg('PromptString',sprintf('Select re-orientation transformation for STEP %d/%d:',ntime,ntimes),'ListSize',[300 200],'SelectionMode','single','ListString',opts(:,1));
            else [treorient,tok] = listdlg('PromptString','Select re-orientation transformation:','ListSize',[300 200],'SelectionMode','single','ListString',opts(:,1));
            end
            if isempty(treorient), return; end
            reorient{ntime}=opts{treorient,2};
            if isequal(reorient{ntime},1)
                answ=conn_menu_inputdlg('Enter affine transformation matrix (4x4 values)','conn_setup_preproc',1,{mat2str(eye(4))});
                if isempty(answ), return; end
                answ=str2num(answ{1});
                reorient{ntime}=answ;
            elseif isequal(reorient{ntime},2)
                [tfilename1,tfilename2]=conn_fileutils('uigetfile','*.mat','Select file');
                if ~ischar(tfilename1), return; end
                filename=fullfile(tfilename2,tfilename1);
                reorient{ntime}=filename;
            elseif isequal(reorient{ntime},3)
                [tfilename1,tfilename2]=conn_fileutils('uigetfile','*.nii','Select file');
                if ~ischar(tfilename1), return; end
                filename=fullfile(tfilename2,tfilename1);
                reorient{ntime}=filename;
            elseif isa(reorient{ntime},'function_handle'),
                answ=conn_menu_inputdlg('Angular rotation (in degrees)','conn_setup_preproc',1,{num2str(90)});
                if isempty(answ), return; end
                answ=str2num(answ{1});
                reorient{ntime}=reorient{ntime}(answ/180*pi);
            end
        end
    end
end

if any(ismember('functional_art',lSTEPS))
    if ~isempty(art_thresholds)
        art_global_threshold=art_thresholds(1);
        art_motion_threshold=art_thresholds(2);
        if numel(art_thresholds)>=3, art_use_diff_global=art_thresholds(3); end
        if numel(art_thresholds)>=4, art_use_diff_motion=art_thresholds(4); end
        if numel(art_thresholds)>=5, art_use_norms=art_thresholds(5); end
        if numel(art_thresholds)>=6, art_force_interactive=art_thresholds(6); end
        if numel(art_thresholds)>=7&&~isnan(art_thresholds(7)), art_motion_threshold(2)=art_thresholds(7); end
        if numel(art_thresholds)>=8, art_drop_flag=art_thresholds(8); end
    end
    if isempty(art_thresholds)||dogui
        thfig=figure('units','norm','position',[.4,.4,.3,.4],'color',1*[1 1 1],'name','Functional outlier detection settings','numbertitle','off','menubar','none');
        ht0=uicontrol('style','popupmenu','units','norm','position',[.05,.8,.9,.1],'string',{'Use liberal settings (99th percentiles in normative sample)','Use intermediate settings (97th percentiles in normative sample)','Use conservative settings (95th percentiles in normative sample)','Edit settings','Edit settings interactively (ART gui)'},'value',1+art_useconservative,'backgroundcolor',1*[1 1 1]);
        ht1a=uicontrol('style','text','units','norm','position',[.05,.7,.9,.05],'string','Global-signal z-value threshold','backgroundcolor',1*[1 1 1]);
        ht1=uicontrol('style','edit','units','norm','position',[.05,.6,.9,.1],'string',num2str(art_global_threshold));
        ht2a=uicontrol('style','text','units','norm','position',[.05,.5,.9,.05],'string','Subject-motion mm threshold','backgroundcolor',1*[1 1 1]);
        ht2=uicontrol('style','edit','units','norm','position',[.05,.4,.9,.1],'string',num2str(art_motion_threshold));
        ht3a=uicontrol('style','checkbox','units','norm','position',[.05,.3,.4,.05],'string','Use diff global','value',art_use_diff_global,'backgroundcolor',1*[1 1 1],'tooltipstring','Global-signal threshold based on scan-to-scan changes in global BOLD signal');
        ht3b=uicontrol('style','checkbox','units','norm','position',[.05,.25,.4,.05],'string','Use abs global','value',~art_use_diff_global,'backgroundcolor',1*[1 1 1],'tooltipstring','Global-signal threshold based on absolute global BOLD signal values');
        ht3c=uicontrol('style','checkbox','units','norm','position',[.05,.20,.4,.05],'string','Drop first scan(s)','value',art_drop_flag>0,'backgroundcolor',1*[1 1 1],'userdata',art_drop_flag,'tooltipstring','Flags first scan(s) in each session for removal');
        ht4a=uicontrol('style','checkbox','units','norm','position',[.55,.3,.4,.05],'string','Use diff motion','value',art_use_diff_motion,'backgroundcolor',1*[1 1 1],'tooltipstring','Subject-motion threshold based on scan-to-scan changes in motion parameters');
        ht4b=uicontrol('style','checkbox','units','norm','position',[.55,.25,.4,.05],'string','Use abs motion','value',~art_use_diff_motion,'backgroundcolor',1*[1 1 1],'tooltipstring','Subject-motion threshold based on absolute motion parameter values');
        ht5=uicontrol('style','checkbox','units','norm','position',[.55,.2,.9,.05],'string','Use comp motion','value',art_use_norms,'backgroundcolor',1*[1 1 1],'tooltipstring','Subject-motion threshold based on composite motion measure');
        uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.10],'callback','uiresume');
        uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)');
        set([ht1a ht1 ht2a ht2 ht3a ht4a ht3b ht3c ht4b ht5],'enable','off');
        set(ht0,'callback','h=get(gcbo,''userdata''); switch get(gcbo,''value''), case 1, set(h.handles,''enable'',''off''); set(h.handles(2),''string'',num2str(h.default{1}(1))); set(h.handles(4),''string'',num2str(h.default{2}(1))); set(h.handles([5:6 10]),''value'',1); set(h.handles([7 9]),''value'',0); case 2, set(h.handles,''enable'',''off''); set(h.handles(2),''string'',num2str(h.default{1}(2))); set(h.handles(4),''string'',num2str(h.default{2}(2))); set(h.handles([5:6 10]),''value'',1); set(h.handles([7 9]),''value'',0); case 3, set(h.handles,''enable'',''off''); set(h.handles(2),''string'',num2str(h.default{1}(3))); set(h.handles(4),''string'',num2str(h.default{2}(3))); set(h.handles([5:6 10]),''value'',1); set(h.handles([7 9]),''value'',0); case 4, set(h.handles,''enable'',''on''); case 5, set(h.handles,''enable'',''off''); end;','userdata',struct('handles',[ht1a ht1 ht2a ht2 ht3a ht4a ht3b ht3c ht4b ht5],'default',{{art_global_thresholds, art_motion_thresholds}}));
        %@(varargin)set([ht1a ht1 ht2a ht2 ht3a ht4a ht3b ht4b ht5],'enable',subsref({'on','off'},struct('type','{}','subs',{{1+(get(gcbo,'value')~=3)}}))));
        set(ht5,'callback','h=get(gcbo,''userdata''); temp=str2num(get(h.handles(4),''string'')); if get(gcbo,''value''), set(h.handles(3),''string'',''Subject-motion mm threshold''); temp=temp(1); else, set(h.handles(3),''string'',''Subject-motion translation/rotation thresholds [mm, rad]''); if numel(temp)<2, temp=[temp .02]; end; end; set(h.handles(4),''string'',mat2str(temp));','userdata',struct('handles',[ht1a ht1 ht2a ht2 ht3a ht4a ht3b ht3c ht4b ht5],'default',{{art_global_thresholds, art_motion_thresholds}}));
        set(ht3a,'callback',@(varargin)set(ht3b,'value',~get(gcbo,'value')));
        set(ht3b,'callback',@(varargin)set(ht3a,'value',~get(gcbo,'value')));
        set(ht3c,'callback','v=get(gcbo,''value''); if v, v=str2double(conn_menu_inputdlg({''Number of initial scans to remove''},'''',1,{num2str(get(gcbo,''userdata''))})); if isempty(v), v=0; end; end; set(gcbo,''value'',v>0); if v>0, set(gcbo,''userdata'',v); end');
        set(ht4a,'callback',@(varargin)set(ht4b,'value',~get(gcbo,'value')));
        set(ht4b,'callback',@(varargin)set(ht4a,'value',~get(gcbo,'value')));
        uiwait(thfig);
        if ~ishandle(thfig), return; end
        art_global_threshold=str2num(get(ht1,'string'));
        temp=str2num(get(ht2,'string'));
        art_motion_threshold=temp;
        art_use_diff_global=get(ht3a,'value');
        art_use_diff_motion=get(ht4a,'value');
        art_use_norms=get(ht5,'value');
        if get(ht3c,'value'), art_drop_flag=get(ht3c,'userdata'); else art_drop_flag=0; end
        art_force_interactive=get(ht0,'value')==5;
        delete(thfig);
        drawnow;
        if numel(art_motion_threshold)<2, art_thresholds=[art_global_threshold(1) art_motion_threshold(1) art_use_diff_global(1) art_use_diff_motion(1) art_use_norms(1) art_force_interactive(1) nan art_drop_flag(1)];
        else art_thresholds=[art_global_threshold(1) art_motion_threshold(1) art_use_diff_global(1) art_use_diff_motion(1) art_use_norms(1) art_force_interactive(1) art_motion_threshold(2) art_drop_flag(1)];
        end
        %answ=conn_menu_inputdlg({'Enter scan-to-scan global signal z-value threshold','Enter scan-to-scan composite motion mm threshold'},'conn_setup_preproc',1,{num2str(art_global_threshold),num2str(art_motion_threshold)});
        %if isempty(answ), return; end
        %art_global_threshold=str2num(answ{1});
        %art_motion_threshold=str2num(answ{2});
    end
end
