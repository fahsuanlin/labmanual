close all; clear all;

data_path={
    '/Users/fhlin/workspace/eegmri_memory/s006/resting_data/unpack/bold/005';
%	'/space_lin2/fhlin/eegmri_memory/s006/resting_data/unpack/bold/005';
    };

reg_path={
    '/Users/fhlin/workspace/eegmri_memory/s006/resting_analysis';
%        '/space_lin2/fhlin/eegmri_memory/s006/resting_analysis';
    };

mc_file='/Users/fhlin/workspace/eegmri_memory/s006/resting_analysis/mc_regressor.mat';
%mc_file='/space_lin2/fhlin/eegmri_memory/s006/resting_analysis/mc_regressor.mat';
    


fstem={
    'sfmcprstc';
    };

file_register='./register.dat';

file_aseg={
    'native_hippo_regressors_rest_hippo_left.mat';
    'native_hippo_regressors_rest_hippo_right.mat';
    };

var_aseg={
    'regressor_hippo';
    'regressor_hippo';
    };

output_aseg_stem={
    'hippo_left';
    'hippo_right';
    };


output_stem={
    'hippo_fconn_native_vol_aseg_120423';
    };

TR=2.0; %second

n_dummy=5;
flag_gavg=1;

confound_polynomial_order=1;
confound_period=[];


subject={
    's006';
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
roi_lh={};
roi_rh={};

for aseg_idx=1:length(file_aseg)

    for f_idx=1:length(fstem)
        valid_subj_idx=[];
        for d_idx=1:length(subject)
            fprintf('[%s]...(%04d|%04d)....\r',subject{d_idx},d_idx,length(subject));

            roi=[];
            STC=[];

            fn=sprintf('%s/%s.nii',data_path{d_idx},fstem{f_idx});
            if(exist(fn))
                vol=MRIread(fn);

                %remove dummy scans
                STC=vol.vol;
                STC=reshape(STC,[size(vol.vol,1)*size(vol.vol,2)*size(vol.vol,3),size(vol.vol,4)]);

                STC(:,1:n_dummy)=[];
                STC(:,end-n_dummy+1:end)=[];

                flag_fe=1;
            else
                flag_fe=0;
            end;

            n_confound=0;
            timepoints=size(STC,2);
            confound=[];
            for j=1:confound_polynomial_order
                n_confound=n_confound+1;
                confound(1:timepoints,n_confound)=([0:1/(timepoints-1):1].^(j))';
            end;
            timeVec=[0:TR:TR*(timepoints-1)];
            for j=1:length(confound_period)
                n_confound=n_confound+1;
                confound(1:timepoints,n_confound)=cos(2.*pi./confound_period(j).*timeVec)';
                n_confound=n_confound+1;
                confound(1:timepoints,n_confound)=sin(2.*pi./confound_period(j).*timeVec)';
            end;

            if(flag_fe)

                valid_subj_idx=cat(1,valid_subj_idx,d_idx);

                %fn=sprintf('%s/regressor_wm_ventrical_%s..mat',reg_path{d_idx},subject{d_idx});
                fn=sprintf('%s/regressor_wm_ventrical_005.mat',reg_path{d_idx});
                if(exist(fn))
                    D_reg=[];
                    load(fn);
                    D_reg(:,1)=regressor_ventricle(1:end);
                    D_reg(:,2)=regressor_wm(1:end);
                    D_reg(1:n_dummy,:)=[];
                    D_reg(end-n_dummy+1:end,:)=[];
                else
                    D_reg=[];
                end;

                D0=[];
                if(~isempty(mc_file))
                    load(mc_file);
                    D0=mc_regressor{f_idx};
                end;


                D=cat(2,D0);
                D(1:n_dummy,:)=[];
                D(end-n_dummy+1:end,:)=[];


                %remove global mean
                D=cat(2,D,ones(size(STC,2),1));
                if(~isempty(D_reg))
                    D=cat(2,D,D_reg);
                end;
                if(flag_gavg);
                    D=cat(2,D,mean(STC,1)');
                end;
                if(~isempty(confound))
                    D=cat(2,D,confound);
                end;
                STC=(STC'-D*(inv(D'*D)*D'*STC')).';

                load(file_aseg{aseg_idx});


                eval(sprintf('reg_aseg=%s;',var_aseg{aseg_idx}));

                %remove dummy scans
                reg_aseg(1:n_dummy)=[];
                reg_aseg(end-n_dummy+1:end)=[];

                reg_aseg=(reg_aseg(:)-D*(inv(D'*D)*D'*reg_aseg(:))).';


                rr=etc_corrcoef(reg_aseg(:),STC');
                z(d_idx,:)=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(STC,2)/(2/TR*2.34)-3));

            end;
        end;

        fprintf('\n');


        z_avg=mean(z(valid_subj_idx,:),1);
        z_std=std(z(valid_subj_idx,:),0,1);


        if(flag_gavg)
            output_stem_now{f_idx}=sprintf('%s_gavg',output_stem{f_idx});
        else
            output_stem_now{f_idx}=output_stem{f_idx};
        end;


        fconn=vol;
        fconn.vol=reshape(z,[size(vol.vol,1),size(vol.vol,2),size(vol.vol,3)]);
        fconn.nframes=1;

        fn=sprintf('%s_%s.nii',output_stem_now{f_idx},output_aseg_stem{aseg_idx});
        MRIwrite(fconn,fn);


        cmd=sprintf('!mri_vol2vol --reg %s --mov %s --fstarg  --o %s_%s-anat.nii',file_register,fn,output_stem_now{f_idx},output_aseg_stem{aseg_idx});
        eval(cmd);
    end;

end;
