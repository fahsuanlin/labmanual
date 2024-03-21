close all; clear all;

root_path='/Users/fhlin/workspace/seeg';

freesurfer_home=getenv('FREESURFER_HOME');

file_hippo_aseg={
    'aseg.mgz';
    };

setenv('SUBJECTS_DIR','/Users/fhlin/workspace/seeg/subjects');

subjects_dir=getenv('SUBJECTS_DIR');

fstem={
    'sfmcprstc';
    };

fmri_vol_data={
    '032';
    };

output_fstem={
    'native_hippo_regressors';
    }

error_subject={};

subject={
    's026';
    };

for stem_idx=1:length(fstem)
    for f_idx=1:length(subject)
        fprintf('[%s]...[%04d||%04d]...\r',subject{f_idx},f_idx,length(subject));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %create a registration matrix from the native subject to the target subject
        %"fsaverage". This registration matrix will be name as
        %"native2fsaverage.dat".
        %eval('!fslregister --s fsaverage --mov bold/004/f.nii --reg ./native2fsaverage.dat --initxfm --maxangle 70');
        %eval(sprintf('!fslregister --s %s --mov %s --reg %s --initxfm --maxangle 70',target_subject, file_register_source{f_idx}, file_register));

        %apply the "inverse" of the registration such that the aparc+aseg.mgz from
        %"fsaverage" will be transformed to the native subject's anatomical space.
        %The transformed apart+aseg file will be named as "aparc+aseg_f.nii".
        %eval('!mri_vol2vol --mov bold/004/fmc.nii --targ $SUBJECTS_DIR/fsaverage/mri/aparc+aseg.mgz --reg native2fsaverage.dat --o aparc+aseg_f.nii --inv --interp nearest');

        %registration has been done previously
        file_register=sprintf('%s/%s/fmri_data/unpack/register.dat',root_path,subject{f_idx});

        for v_idx=1:length(fmri_vol_data)
            file_mov=sprintf('%s/%s/fmri_data/unpack/bold/%s/%s.nii',root_path,subject{f_idx},fmri_vol_data{v_idx},fstem{stem_idx});

            if(exist(file_mov))

                %file_output=sprintf('/space_lin1/hcp/analysis/%s_%s.mat',subject{f_idx},output_fstem{stem_idx});
                %    eval(sprintf('!mri_vol2vol --mov %s --targ $SUBJECTS_DIR/%s/mri/aparc+aseg.mgz --reg %s --o %s --inv --interp nearest',file_mov,  target_subject, file_register, file_aseg));

                try

                    region_index={
                        [17]; %left hippo
                        [53]; %right hippo
                        %[203 233 235 237 239 241 243 245 211]; %head
                        %[234 236 238 240 242 244 246 212]; %body
                        %[226]; %tail
                        };

                    region_stem={
                        'hippo_left';
                        'hippo_right';
                        %'hippo_head';
                        %'hippo_body';
                        %'hippo_tail';
                        };

                    %     203: parasubiculum
                    %     211: HATA
                    %     212: fimbria
                    %     215: hippocampal fissure
                    %     226: HP_tail
                    %     233: presubiculum-head
                    %     234: presubiculum-body
                    %     235: subiculum-head
                    %     236: subiculum-body
                    %     237: CA1-head
                    %     238: CA1-head
                    %     239: CA3-head
                    %     240: CA3-body
                    %     241: CA4-head
                    %     242: CA4-body
                    %     243: GC-ML-DG-head
                    %     244: GC-ML-DG-body
                    %     245: molecular_layer_HP-head
                    %     246: molecular_layer_HP-body
                    %     7001: lateral nucleus
                    %     7003: basal nucleus
                    %     7005: central nucleus
                    %     7006: medial ncleaus
                    %     7007: cortical-nucleus
                    %     7008: accesory basal nucleus
                    %     7009: corticoamygdaloid transitio
                    %     7010: anterior amygdaloid area AAA
                    %     7015: paralaminar nucleus

                    for region_idx=1:length(region_index)

                        file_output=sprintf('%s/%s/fmri_analysis/%s_%s.mat',root_path,subject{f_idx},output_fstem{stem_idx},region_stem{region_idx});

                        tmp=[];

                        %                    if(~exist(file_output))

                        for aseg_idx=1:length(file_hippo_aseg)

                            %eval(sprintf('!mri_vol2vol --mov %s --targ $SUBJECTS_DIR/%s/mri/%s --reg %s --o /space_lin1/hcp/analysis/%s_%s --inv --interp nearest',file_mov,  subject{f_idx}, file_hippo_aseg{aseg_idx}, file_register, subject{f_idx}, file_hippo_aseg{aseg_idx}));

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %
                            %get white-matter and ventrical from "FreeSurferColorLUT.txt"

                            [dummy,fstm]=fileparts(file_hippo_aseg{aseg_idx});
                            %eval(sprintf('!mri_vol2vol --nearest --mov %s/%s/mri/%s --targ %s/average/mni305.cor.mgz  --xfm %s/%s/mri/transforms/talairach.xfm   --o %s/resting_analysis/%s_%s-mni305.mgz', subjects_dir, subject{f_idx}, file_hippo_aseg{aseg_idx}, freesurfer_home, subjects_dir, subject{f_idx}, root_path, subject{f_idx},fstm));
                            %eval(sprintf('!mri_vol2vol --mov %s/%s/mri/%s --reg register_id.dat   --tal --talres 2 --o %s/analysis/%s_%s-mni305-2mm.mgz', subjects_dir, subject{f_idx}, file_hippo_aseg{aseg_idx}, root_path, subject{f_idx},fstm));

                            %eval(sprintf('!mri_vol2vol --mov %s/%s/mri/%s --reg %s/average/mni305.cor.subfov2.reg   --tal --talres 2 --o %s/analysis/%s_%s-mni305-2mm.mgz', subjects_dir, subject{f_idx}, file_hippo_aseg{aseg_idx}, freesurfer_home, root_path, subject{f_idx},fstm));

                            %eval(sprintf('!mri_vol2vol --mov %s/%s/mri/%s --reg %s/average/mni305.cor.subfov2.reg   --tal --talres 2 --o %s/analysis/%s_%s-mni305-2mm.mgz', subjects_dir, subject{f_idx}, file_hippo_aseg{aseg_idx}, freesurfer_home, root_path, subject{f_idx},fstm));


                            eval(sprintf('!mri_vol2vol --nearest --mov %s --targ %s/%s/mri/%s --reg %s --o  %s/%s/fmri_analysis/%s-aseg.mgz --inv ',file_mov, subjects_dir,subject{f_idx}, file_hippo_aseg{v_idx}, file_register,root_path,subject{f_idx},fmri_vol_data{v_idx}));


                            d_aseg=MRIread(sprintf('%s/%s/fmri_analysis/%s-aseg.mgz',root_path, subject{f_idx},fmri_vol_data{v_idx}));

                            %d_regression=MRIread(file_regression_source);
                            %v_regression=d_regression.vol;
                            d=MRIread(file_mov);
                            acc=d.vol;
                            dim=size(acc);
                            %v_regression=reshape(acc,[dim(1)*dim(2)*dim(3), d_regression.nframes]);
                            %fprintf('functional data: [%d] voxels x [%d] time points\n',dim(1)*dim(2)*dim(3), d_regression.nframes);
                            v_regression=reshape(acc,[dim(1)*dim(2)*dim(3), dim(4)]);
                            fprintf('functional data: [%d] voxels x [%d] time points\n',dim(1)*dim(2)*dim(3), dim(4));

                            %vol = etc_MRIvol2vol(d_aseg,d, inv(d_aseg.tkrvox2ras*inv(d_aseg.vox2ras)*d.tkrvox2ras*inv(d.vox2ras)));
                            v_aseg=d_aseg.vol;
                            for roi_idx=1:length(region_index{region_idx})
                                idx=find(v_aseg(:)==region_index{region_idx}(roi_idx));
                                fprintf('%s: index [%d]: [%d] voxels...\n',region_stem{region_idx},region_index{region_idx}(roi_idx),length(idx));
                                if(isempty(tmp))
                                    tmp=v_regression(idx,:);
                                else
                                    tmp=cat(1,tmp,v_regression(idx,:));
                                end;
                            end;
                        end;
                        regressor_hippo=mean(tmp,1);

                        save(file_output,'regressor_hippo');

                        %eval(sprintf('!rm %s/analysis/%s_%s',root_path,subject{f_idx},file_aseg));
                        %                   end;

                    end;
                catch
                    error_subject{end+1}=subject{f_idx};
                end;
            else
                fprintf('regressor [%s] existed!\n',file_output);
            end;
        end;
    end;
    fprintf('\n');

    fprintf('error subject:\n');
    for ef_idx=1:length(error_subject)
        fprintf('[%s]\n',error_subject{ef_idx});
    end;
end;

fprintf('DONE!\n');