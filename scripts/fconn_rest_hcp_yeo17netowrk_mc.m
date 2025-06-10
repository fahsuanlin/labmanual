close all; clear all;

data_path='/space_lin1/hcp';
%data_path='/Users/fhlin/workspace/hcp/';

fstem={
    '2_fsaverage_rfMRI_REST1_LR_hp2000_clean';
    '2_fsaverage_rfMRI_REST1_RL_hp2000_clean';
    '2_fsaverage_rfMRI_REST2_LR_hp2000_clean';
    '2_fsaverage_rfMRI_REST2_RL_hp2000_clean';
    };

fstem_reg={
    %100206_native_hippo_regressors_aseg_rest1_LR_hippo_left.mat
    'native_hippo_regressors_aseg_rest1_LR_hippo_left','native_hippo_regressors_aseg_rest1_LR_hippo_right';
    'native_hippo_regressors_aseg_rest1_RL_hippo_left','native_hippo_regressors_aseg_rest1_RL_hippo_right';
    'native_hippo_regressors_aseg_rest2_LR_hippo_left','native_hippo_regressors_aseg_rest2_LR_hippo_right';
    'native_hippo_regressors_aseg_rest2_RL_hippo_left','native_hippo_regressors_aseg_rest2_RL_hippo_right';
    };

fstem_confound={
    %	'rfMRI_REST1_LR_hp2000_clean_100206_regressors.mat';
    %	'rfMRI_REST1_LR_hp2000_clean';
    %        'rfMRI_REST1_RL_hp2000_clean';
    %        'rfMRI_REST2_LR_hp2000_clean';
    %        'rfMRI_REST2_RL_hp2000_clean';
    'native_regressors_rest1_LR';
    'native_regressors_rest1_RL';
    'native_regressors_rest2_LR';
    'native_regressors_rest2_RL';
    };

TR=0.72; %second

n_dummy=0;
flag_gavg=0;

subject='';
d=textread('subject_list_all.txt');
for d_idx=1:length(d)
    subject{d_idx}=num2str(d(d_idx));
end;
%subject=[];
%subject{1}='119732';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get annotation indices for left and right hemispheres....
file_annot=sprintf('%s/subjects/fsaverage/label/%s',getenv('FREESURFER_HOME'),'lh.Yeo2011_17Networks_N1000.annot');
[label_vertex label_value label_ctab] = read_annotation(file_annot);
for label_idx=2:length(label_ctab.struct_names)
    roi_lh(label_idx).label=find(label_ctab.table(label_idx,5)==label_value)-1; %0-based index
    [dummy,roi_lh(label_idx).stc_idx]=intersect([0:10241],roi_lh(label_idx).label);
    roi_lh(label_idx).name=label_ctab.struct_names{label_idx};
end;
file_annot=sprintf('%s/subjects/fsaverage/label/%s',getenv('FREESURFER_HOME'),'rh.Yeo2011_17Networks_N1000.annot');
[label_vertex label_value label_ctab] = read_annotation(file_annot);
for label_idx=2:length(label_ctab.struct_names)
    roi_rh(label_idx).label=find(label_ctab.table(label_idx,5)==label_value)-1; %0-based index
    [dummy,roi_rh(label_idx).stc_idx]=intersect([0:10241],roi_rh(label_idx).label);
    roi_rh(label_idx).name=label_ctab.struct_names{label_idx};
end;

%concatenate both hemispheres
for label_idx=2:length(label_ctab.struct_names)
    roi_bh(label_idx).stc_idx=cat(1,roi_lh(label_idx).stc_idx,roi_rh(label_idx).stc_idx+10242);
    roi_bh(label_idx).name=label_ctab.struct_names{label_idx};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for d_idx=1:length(subject)
    fprintf('[%s]...(%04d|%04d)....\r',subject{d_idx},d_idx,length(subject));

    fconn{d_idx}.subject=subject{d_idx};

    for f_idx=1:length(fstem)
        valid_subj_idx=[];

        roi=[];
        STC=[];

        ROI_avg=[];
        ROI_std=[];
        for hemi_idx=1:2
            switch hemi_idx
                case 1
                    hemi_str='lh';
                case 2
                    hemi_str='rh';
            end;

            fn=sprintf('%s/%s/analysis/%s_%s-%s.stc',data_path,subject{d_idx},subject{d_idx},fstem{f_idx},hemi_str);
            if(exist(fn))
                fprintf('data file [%s] found ...\n',fn);

                fconn{d_idx}.data_stem{f_idx}=fstem{f_idx};
                [stc{hemi_idx},v{hemi_idx},d0,d1,timeVec]=inverse_read_stc(fn);

                %remove dummy scans
                stc{hemi_idx}(:,1:n_dummy)=[];
                stc{hemi_idx}(:,end-n_dummy+1:end)=[];

                STC=cat(1,STC,stc{hemi_idx});
                flag_fe=1;
            else
                flag_fe=0;
            end;
        end;
        if(flag_fe)

            valid_subj_idx=cat(1,valid_subj_idx,d_idx);

            %fn=sprintf('%s/analysis/%s_regressors.mat',data_path,subject{d_idx});
            %fn=sprintf('%s/analysis/%s_%s_regressors.mat',data_path,fstem_confound{f_idx},subject{d_idx});
            fn=sprintf('%s/analysis/%s_%s.mat',data_path,subject{d_idx},fstem_confound{f_idx});

            if(exist(fn))
                D_reg=[];
                fprintf('confound file [%s] found ...\n',fn);

                fconn{d_idx}.file_confound{f_idx}=fn;

                load(fn);
                D_reg(:,1)=regressor_ventricle(1:end-1);
                D_reg(:,2)=regressor_wm(1:end-1);
                D_reg(1:n_dummy,:)=[];
                D_reg(end-n_dummy+1:end,:)=[];
            else
                D_reg=[];
            end;

            %motion regressors
            fn=sprintf('%s/%s/Preprocessed/Movement_Regressors.txt',data_path,subject{d_idx});

            if(exist(fn))
                fprintf('confound file [Movement_Regressors.txt] found ...\n');
                try
                    D_mc=textread(fn);
                    D_mc(:,end)=[];
                    D_mc(end,:)=[];
                catch
                    fprintf('error in reading motion regressor!\n');
                    D_mc=[];
                end;
            else
                D_mc=[];
            end;

            try
                D_reg=cat(2,D_reg,D_mc);
            catch
                fprintf('error in concatenating WM/ventricle regressors [%s] and motion regressors [%s]!\n',mat2str(size(D_reg)),mat2str(size(D_mc)));
            end;
            %remove global mean
            D=ones(size(STC,2),1);
            D(:,2)=[0:size(D,1)-1]./(size(D,1)-1);
            if(~isempty(D_reg))
                D=cat(2,D,D_reg);
            end;
            if(flag_gavg);
                D=cat(2,D,mean(STC,1)');
            end;
            STC=(STC'-D*(inv(D'*D)*D'*STC')).';

            stc_hemi{1}=STC(1:length(v{1}),:);
            stc_hemi{2}=STC(length(v{1})+1:end,:);


            hippo_reg=[];
            %for r_idx=1:length(fstem_reg(f_idx,:))
            for r_idx=2:length(label_ctab.struct_names)

                if(~isempty(roi_bh(r_idx).stc_idx))

                    data=STC(roi_bh(r_idx).stc_idx,:);

                    cc=corrcoef(data');
                    idx = triu(true(size(cc)), 1);  % logical matrix with 1s above diagonal
                    rr = cc(idx);
                    zz = 0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(STC,2)/(2/TR*2.34)-3));
                    %fconn{d_idx}.fconn{f_idx,r_idx}.fconn_corrcoef=cc(idx);
                    fconn{d_idx}.fconn{f_idx,r_idx}.fconn_z_avg(:,hemi_idx)=mean(zz);
                    fconn{d_idx}.fconn{f_idx,r_idx}.fconn_z_std(:,hemi_idx)=std(zz);

                    fconn{d_idx}.fconn{f_idx,r_idx}.name=label_ctab.struct_names{label_idx};
                end;
            end;


        end;
    end;

    fprintf('\n');

end;
save fconn_rest_hcp_yeo17network_aseg_mc.mat fconn roi_lh roi_rh roi_bh