close all; clear all;

load electrodes_to_labels_061721.mat

roi_index=[1 4];

sfreq=2048; %Hz
freqVec=[6:10,12:2:20,25:5:60,70:10:120,140:20:300,350:50:500];

clim=[-8 8];

output_stem='electrode_tfr_061721';

%colormap
a1=ones(1,64);
ag=linspace(0,1,64);
rr=[a1,a1,fliplr(ag)];
gg=[ag,a1,fliplr(ag)];
bb=[ag,a1,a1];
cmap=flipud([rr(:),gg(:),bb(:)]);

for roi_idx=1:length(roi_index)
    for cond_idx=1:length(roi(roi_index(roi_idx)).erf_electrode_min_dist)
        tmp=squeeze(roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).data);
        
        tfr(roi_idx,cond_idx).data=tmp;
        
        timeVec=roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).timeVec;
        
        tfr(roi_idx,cond_idx).timeVec=timeVec;
        tfr(roi_idx,cond_idx).freqVec=freqVec;
        tfr(roi_idx,cond_idx).name=roi(roi_index(roi_idx)).name;
        tfr(roi_idx,cond_idx).trig_str=roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).trig_str;
        
        base_idx=find(timeVec<0);
        
        for freq_idx=1:length(freqVec)
            fprintf('condition [%s] f=%1.0f (Hz)...\r',roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).trig_str, freqVec(freq_idx));
            for trial_idx=1:size(tmp,2)
                tfr_abs(:,trial_idx)=abs(inverse_waveletcoef(freqVec(freq_idx),tmp(:,trial_idx)',sfreq,7));
            end;
            tfr_tmp=tfr_abs.^2;
            
            tfr_tmp=trimmean(tfr_tmp,5,'round',2);
            
            tfr_tmp=log10(tfr_tmp./repmat(mean(tfr_tmp(base_idx,:),1),[size(tfr_tmp,1),1]));
            
            roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).tfr(:,freq_idx)=tfr_tmp;
        end;
        tfr(roi_idx,cond_idx).log10tfr=roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).tfr;
        
        roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).tfr_Z=etc_z(roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).tfr,base_idx,'flag_baseline_correct',1,'dim',1);
        roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).freqVec=freqVec;
        
        tfr(roi_idx,cond_idx).log10tfr_Z=roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).tfr_Z;
        fprintf('\n');
        
        figure;
        h=pcolor(timeVec,log10(freqVec),(roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).tfr_Z)'); shading flat; set(gca,'clim',clim);
        
        hold on;
        l=line([0 0],get(gca,'ylim')); set(l,'color',[1 1 1].*0.5,'LineWidth',2);
        
        cb=colorbar; set(cb,'ylim',clim,'ytick',[min(clim) max(clim)],'ticklength',3*get(cb,'ticklength'),'fontsize',16,'fontname','helvetica');
           
        set(gca,'ytick',log10([10 20 50 100 200 500]),'yticklabel',{'10','20','50','100','200','500'});
        colormap(cmap)
        xlabel('time (ms)');
        ylabel('frequency (Hz)');
        etc_plotstyle;
        
        
        hgexport(gcf,sprintf('%s_%s_%s',output_stem,roi(roi_index(roi_idx)).name,roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).trig_str), hgexport('factorystyle'),'Format','png');

        fn{roi_idx,cond_idx}=sprintf('%s_%s_%s.png',output_stem,roi(roi_index(roi_idx)).name,roi(roi_index(roi_idx)).erf_electrode_min_dist(cond_idx).trig_str);
    end;
    
end;

for roi_idx=1:length(roi_index)
    for cond_idx=1:length(roi(roi_index(roi_idx)).erf_electrode_min_dist)
        mont_file{cond_idx}=fn{roi_idx,cond_idx};
    end;
    [ccimg]=etc_montage([1 2 3],'mont_file',mont_file);
    
    fn_out=sprintf('%s_%s_composite.png',output_stem,roi(roi_index(roi_idx)).name);
    
    imwrite(double(ccimg./255),fn_out);
end;

save electrode_tfr_061721.mat tfr
