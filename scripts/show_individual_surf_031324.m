close all; clear all;

subject={
    's025';
    's026';
    's031';
    's032';
    's034';
%    's041';
%    's045';
%    's047';
    };
    


cond_stem={
    'fmri_surf_soa_glm_h01_beta';
    'fmri_surf_soa_glm_h02_beta';
    'fmri_surf_soa_glm_h03_beta';
};

cond_output_stem={
    'h01';
    'h02';
    'h03';
};


cond_stem_str={
    'A';
    'V';
    'AV';
    };

root_dir='/Users/fhlin/workspace/seeg';

output_stem='fmri_surf_031324';

for subj_idx=1:length(subject)
    for cond_idx=1:length(cond_stem)
        for hemi_idx=1:2
            switch hemi_idx
                case 1
                    hemi_str='lh';
                case 2
                    hemi_str='rh';
            end;

            [dummy,v,a,b,timeVec]=inverse_read_stc(sprintf('%s/%s/fmri_analysis/%s-%s.stc',root_dir,subject{subj_idx},cond_stem{cond_idx},hemi_str));

            stc(:,cond_idx)=dummy(:,1);

            inverse_write_stc(repmat(stc(:,cond_idx),[1:5]),v,a,b,sprintf('%s_%s_%s-%s.stc',output_stem,subject{subj_idx},cond_output_stem{cond_idx},hemi_str));

        end;

         etc_render_fsbrain_stc({sprintf('%s_%s_%s',output_stem,subject{subj_idx},cond_output_stem{cond_idx})},[3 6],'flag_overlay_pos_only',1);

        print('-dpng',sprintf('%s_%s_%s.png',output_stem,subject{subj_idx},cond_stem_str{cond_idx}));

        close

    end;
end;
