close all; clear all;

file_mat='electrode_tfr_120121.mat';

subject={
    's025';
    's026';
    's027';
    's031';
    's032';
    's033';
    's034';
    's035';
    's036';
    's038';
    's041';
    's045';
    's046';
    's047';
    's048';
    's050';
    's051';
    's054';
    's055';
    's057';
    };

data_path='/Users/fhlin/workspace/seeg';

lh_index=[0 1 1 0 1 1 0 0 0 0 1 1 1 0 1 1 0 0 1 1];
rh_index=[1 0 0 1 0 0 1 1 1 1 1 0 0 1 0 1 1 1 1 1];

clim=[-8 8];

roi_name={
    'aud_early_lh';
    'aud_assoc_lh';
    'tpoj_lh';
    'aud_early_rh';
    'aud_assoc_rh';
    'tpoj_rh';
    };

output_stem='average_electrode_tfr_120121';


%%%%%%%%%
for subj_idx=1:length(subject)
    load(sprintf('%s/%s/analysis/%s',data_path,subject{subj_idx},file_mat));
    
    for roi_idx=1:size(tfr,1)
        for cond_idx=1:size(tfr,2)
            tfr_all(:,:,subj_idx,roi_idx,cond_idx)=tfr(roi_idx,cond_idx).log10tfr;
        end;
    end;
end;
    
timeVec=tfr(1,1).timeVec;
freqVec=tfr(1,1).freqVec;
        
base_idx=find(timeVec<0);

lh_subj=find(lh_index>0.5);
rh_subj=find(rh_index>0.5);


%colormap
a1=ones(1,64);
ag=linspace(0,1,64);
rr=[a1,a1,fliplr(ag)];
gg=[ag,a1,fliplr(ag)];
bb=[ag,a1,a1];
cmap=flipud([rr(:),gg(:),bb(:)]);


for roi_idx=1:size(tfr,1)
    for cond_idx=1:size(tfr,2)
        if(~isempty(findstr(tfr(roi_idx,1).name,'lh')))
            tfr_all_Z(:,:,roi_idx,cond_idx)=etc_z(squeeze(mean(tfr_all(:,:,lh_subj,roi_idx,cond_idx),3)),base_idx,'flag_baseline_correct',1,'dim',1);
        else
            tfr_all_Z(:,:,roi_idx,cond_idx)=etc_z(squeeze(mean(tfr_all(:,:,rh_subj,roi_idx,cond_idx),3)),base_idx,'flag_baseline_correct',1,'dim',1);
        end;
        
        figure;
        h=pcolor(timeVec,log10(freqVec),tfr_all_Z(:,:,roi_idx,cond_idx)'); shading flat; set(gca,'clim',clim);
        
        hold on;
        l=line([0 0],get(gca,'ylim')); set(l,'color',[1 1 1].*0.5,'LineWidth',2);
        
        cb=colorbar; set(cb,'ylim',clim,'ytick',[min(clim) max(clim)],'ticklength',3*get(cb,'ticklength'),'fontsize',16,'fontname','helvetica');
           
        set(gca,'ytick',log10([10 20 50 100 200 500]),'yticklabel',{'10','20','50','100','200','500'});
        colormap(cmap)
        xlabel('time (ms)');
        ylabel('frequency (Hz)');
        etc_plotstyle;

        hgexport(gcf,sprintf('%s_%s_%s',output_stem,roi_name{roi_idx},tfr(roi_idx,cond_idx).trig_str), hgexport('factorystyle'),'Format','png');

        fn{roi_idx,cond_idx}=sprintf('%s_%s_%s.png',output_stem,roi_name{roi_idx},tfr(roi_idx,cond_idx).trig_str);

    end;
end;



for roi_idx=1:size(tfr,1)
    for cond_idx=1:size(tfr,2)
        mont_file{cond_idx}=fn{roi_idx,cond_idx};
    end;
    [ccimg]=etc_montage([1 2 3],'mont_file',mont_file);
    
    fn_out=sprintf('%s_%s_composite.png',output_stem,roi_name{roi_idx});
    
    imwrite(double(ccimg./255),fn_out);
end;
