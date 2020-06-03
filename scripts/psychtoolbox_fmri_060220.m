%
% a back-bone program for generating stimuli and collecting responses by psychtolbox in fMRI experiments
%
% fhlin@jun 2 2020
%

close all; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
Screen('Preference', 'SkipSyncTests', 1);

file_param='param_060220.txt';

output_stem='psychtoolbox_fmri_060220';
output_stem=''; %empty string disables the file saving

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR=1;                   %MRI TR (second), the time between two triggers
n_scan=5;               %total number of MRI scan (one trigger per scan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross hair setup
cross_length=20;
cross_width=1;

% static checkerboard setup
n_ring=16;       %nuber of rings
n_wedge=4*6;    %number of wedges in 2*PI

% dynamic checkerboard setup
flash_rate=8;               %Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize stimulus environment
whichScreen=0;
window=Screen('OpenWindow',0,0,[1 1 400 400]); %a smaller window; good for debugging
%window=Screen('OpenWindow',0,0); %a full window; good for experiment
windowRect=Screen(window,'Rect');

black=BlackIndex(window);
white=WhiteIndex(window);
GRAY=(white+black)/2;
inc=white-GRAY;

Screen(window,'FillRect',GRAY);
Screen(window,'TextSize',24);
Screen(window,'DrawText','initializing...',windowRect(3)/4,windowRect(4)/2,[256 256 256]);
Screen('Flip',window);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visual stimuli setup
xx=windowRect(3);
yy=windowRect(4);
[x,y]=meshgrid(-xx/2:xx/2-1,-yy/2:yy/2-1);
radius=sqrt(x.^2+y.^2);

%correct for magnification factor
mm=max(radius(:));
radius=(radius./mm).^(0.2).*mm;

r_diff=min([xx,yy])/2/n_ring;
R=mod(floor(radius/r_diff),2);
R(find(R==0))=-1;
rr=floor(radius/r_diff);

ang=atan2(y,x).*180./pi;
ang_diff=360/n_wedge;
A=mod(floor(ang/ang_diff),2);
A(find(A==0))=-1;

mask_righthemi=zeros(size(ang));
mask_righthemi(find(ang(:)>-30&ang(:)<30))=1;
mask_lefthemi=zeros(size(ang));
mask_lefthemi(find(ang(:)>150|ang(:)<-150))=1;

mask_bothhemi=mask_righthemi+mask_lefthemi;

mask_inner=zeros(size(ang));
mask_inner(find(rr(:)<=20))=1;
mask_outer=zeros(size(ang));
mask_outer(find(rr(:)>=23))=1;

NumImages=2;
%create enough offscreen windows for each picture in the experiment
%black and white checkerboards
for i = 1:NumImages
    pattern = GRAY+(inc*A.*R.*((-1).^(i))).*mask_bothhemi;
    pattern(round(windowRect(4)/2-cross_width):round(windowRect(4)/2+cross_width),round(windowRect(3)/2-cross_length):round(windowRect(3)/2+cross_length))=white;
    pattern(round(windowRect(4)/2-cross_length):round(windowRect(4)/2+cross_length),round(windowRect(3)/2-cross_width):round(windowRect(3)/2+cross_width))=white;
    offScrPtr(i)   = Screen('MakeTexture',window,pattern);
    
    visual_stim(i).ptr= offScrPtr(i);
    visual_stim(i).id= sprintf('v%d',i);
    
    %%% save images for documentation
    %imwrite(repmat(pattern,[1 1 3]),gray,sprintf('V_091215_%02d.tif',i),'tiff');
end

for i = 1:NumImages
    pattern = GRAY+(inc*A.*R.*((-1).^(i))).*mask_lefthemi;
    pattern(round(windowRect(4)/2-cross_width):round(windowRect(4)/2+cross_width),round(windowRect(3)/2-cross_length):round(windowRect(3)/2+cross_length))=white;
    pattern(round(windowRect(4)/2-cross_length):round(windowRect(4)/2+cross_length),round(windowRect(3)/2-cross_width):round(windowRect(3)/2+cross_width))=white;
    offScrPtr(i+2)   = Screen('MakeTexture',window,pattern);
    
    visual_stim(i+2).ptr= offScrPtr(i+2);
    visual_stim(i+2).id= sprintf('v%d',i+2);

    %%% save images for documentation
    %imwrite(repmat(pattern,[1 1 3]),gray,sprintf('VL_091215_%02d.tif',i),'tiff');
end

for i = 1:NumImages
    pattern = GRAY+(inc*A.*R.*((-1).^(i))).*mask_righthemi;
    pattern(round(windowRect(4)/2-cross_width):round(windowRect(4)/2+cross_width),round(windowRect(3)/2-cross_length):round(windowRect(3)/2+cross_length))=white;
    pattern(round(windowRect(4)/2-cross_length):round(windowRect(4)/2+cross_length),round(windowRect(3)/2-cross_width):round(windowRect(3)/2+cross_width))=white;
    offScrPtr(i+4)   = Screen('MakeTexture',window,pattern);
    
    visual_stim(i+4).ptr= offScrPtr(i+4);
    visual_stim(i+4).id= sprintf('v%d',i+4);

    %%% save images for documentation
    %imwrite(repmat(pattern,[1 1 3]),gray,sprintf('VR_091215_%02d.tif',i),'tiff');
end


cross=ones(size(pattern,1),size(pattern,2)).*GRAY;
cross(round(windowRect(4)/2-cross_width):round(windowRect(4)/2+cross_width),round(windowRect(3)/2-cross_length):round(windowRect(3)/2+cross_length))=white;
cross(round(windowRect(4)/2-cross_length):round(windowRect(4)/2+cross_length),round(windowRect(3)/2-cross_width):round(windowRect(3)/2+cross_width))=white;
offScrPtr(length(offScrPtr)+1)  = Screen('MakeTexture',window,cross);

visual_stim(7).ptr= offScrPtr(end);
visual_stim(7).id= sprintf('vbg');

%%% save images for documentation
%imwrite(repmat(cross,[1 1 3]),gray,sprintf('BK_091215.tif'),'tiff');

text_stim(1).text='a';
text_stim(1).textsize=24;
text_stim(1).textfont='arial';
text_stim(1).textcolor=[1 0 0].*256;
text_stim(1).id='t1';

text_stim(2).text='+';
text_stim(2).textsize=24;
text_stim(2).textfont='arial';
text_stim(2).textcolor=[1 0 1].*256;
text_stim(2).id='t2';

text_stim(3).text='#';
text_stim(3).textsize=24;
text_stim(3).textfont='arial';
text_stim(3).textcolor=[1 1 0].*256;
text_stim(3).id='t3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auditory stimuli setup
[s1,audio_fs]=audioread('sounds/500hz.wav');
[s2,audio_fs]=audioread('sounds/1000hz.wav');
[s3,audio_fs]=audioread('sounds/1500hz.wav');

s1=s1(:);
s2=s2(:);
s3=s3(:);
audio_stim(1).data=s1';
audio_stim(1).id='a1';
audio_stim(2).data=s2';
audio_stim(2).id='a2';
audio_stim(3).data=s3';
audio_stim(3).id='a3';


InitializePsychSound(1); %inidializes sound driver...the 1 pushes for low latency
pahandle = PsychPortAudio('Open', [], [], 1, audio_fs, 1, 0); % opens sound buffer at a different frequency


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finish stimulus setup. notify the user.
HideCursor;
Screen(window,'FillRect',GRAY);
Screen(window,'DrawText','hit any key to start the experiment ...',windowRect(3)/4,windowRect(4)/2,[256 256 256]);
Screen('Flip',window,0,0,2);
KbWait(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_start=GetSecs;

KbState=0;
KeyPressCnt=0;
KeyPressed_time=[];
KeyPressed=[];

Section_timebegin=GetSecs;

flag_cont=1;
time_begin=GetSecs;


[stim_time,stim_type,stim_id]=textread(file_param,'%f%s%s');
visual_stim_all={visual_stim(:).id};
text_stim_all={text_stim(:).id};
audio_stim_all={audio_stim(:).id};

for ii=1:length(stim_id)
    if(strcmp(stim_type(ii),'a'))
        IndexC = strfind(audio_stim_all,stim_id{ii});
        tmp=find(not(cellfun('isempty',IndexC)));    
        if(isempty(tmp))
            fprintf('error! no audio stimulus [%s] in the param prepared!\n',stim_id{ii});
            Screen('CloseAll');
            return;
        else
            stim_idx(ii)=tmp;
        end;
    elseif(strcmp(stim_type(ii),'v'))
        IndexC = strfind(visual_stim_all,stim_id{ii});
        tmp=find(not(cellfun('isempty',IndexC)));    
        if(isempty(tmp))
            fprintf('error! no  visual stimulus [%s] in the param prepared!\n',stim_id{ii});
            Screen('CloseAll');
            return;
        else
            stim_idx(ii)=tmp;
        end;
    elseif(strcmp(stim_type(ii),'t'))
        IndexC = strfind(text_stim_all,stim_id{ii});
        tmp=find(not(cellfun('isempty',IndexC)));    
        if(isempty(tmp))
            fprintf('error! no  visual stimulus [%s] in the param prepared!\n',stim_id{ii});
            Screen('CloseAll');
            return;
        else
            stim_idx(ii)=tmp;
        end;
    end;
end;
stim_count=1;

check_stim=1;


while(flag_cont)
    
    %check if A/V stim should be played
    if(check_stim)
        if(GetSecs-time_begin>stim_time(stim_count))
            
            switch lower(stim_type{stim_count})
                case 'a' %auditory stimulus
                    PsychPortAudio('FillBuffer', pahandle, audio_stim(stim_idx(stim_count)).data); % loads data into buffer
                    PsychPortAudio('Start', pahandle, 1,0); %starts sound immediatley
                     %PsychPortAudio('Stop', pahandle,1);% Stop sound playback
                case 'v' %visual stimulus
                    Screen('DrawTexture',window,offScrPtr(stim_idx(stim_count)));
                    [t0 StimulusOnsetTime FlipTimestamp]=Screen('Flip',window);
                case 't' %text stimulus
                    fprintf('t');
                    Screen(window,'TextSize',text_stim(stim_idx(stim_count)).textsize);
                    Screen(window,'TextColor',text_stim(stim_idx(stim_count)).textcolor);
                    Screen(window,'DrawText',text_stim(stim_idx(stim_count)).text,windowRect(3)/2,windowRect(4)/2);
                    %Screen(window,'TextSize',text_stim(stim_idx(stim_count)).textsize);
                    %Screen(window,'TextColor',text_stim(stim_idx(stim_count)).textcolor);
                    %
                    %                    Screen('DrawText',window,text_stim(stim_idx(stim_count)).text,windowRect(3)/4,windowRect(4)/2,text_stim(stim_idx(stim_count)).textcolor);
                    Screen('Flip',window,0,0,2);
            end;
            
            stim_count=stim_count+1;
            
            if(stim_count>length(stim_time))
                check_stim=0;
            end;
        end;
    end;
    
    %capture keys
    [keyIsDown, secs, keyCode] = KbCheck(-1);
    if xor(KbState,keyIsDown) && keyIsDown,
        KeyPressCnt=KeyPressCnt+1;
        KeyPressed_time(KeyPressCnt)= [secs-Section_timebegin];
        KeyPressed{KeyPressCnt}={KbName(keyCode)};
    end
    KbState=keyIsDown;
    
    %check if at the end of the experiment
    if((GetSecs-time_begin)>n_scan.*TR) 
        flag_cont=0; 
    end;
end;

PsychPortAudio('Close', pahandle);% Close the audio device:

time_end=GetSecs;

Screen('CloseAll');


if(~isempty(output_stem))
    c=clock;
    fn=sprintf('%s_[%4d_%02d_%02d]_[%02d_%02d_%02d].mat',output_stem,c(1),c(2),c(3),c(4),c(5),round(c(6)));
    fprintf('saving [%s]...\n',fn);
    save(fn,'KeyPressed_time','KeyPressed','stim_time','stim_type','stim_id','stim_idx');
end;
fprintf('DONE!\n');

