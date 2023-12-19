function trigger_monitor_readSerialData(src,~)
data = readline(src);
if(isempty(strfind(data,'Output')))

    fprintf("TTL read! %s", data);

    clock_now=datetime("now");
    laps_time=clock_now-src.UserData.PortObj.input_clock;

    src.UserData.LabelTTLInputCounter.Text=sprintf('%04d',str2num(src.UserData.LabelTTLInputCounter.Text)+1);
    src.UserData.PortObj.TTL_input_counter=src.UserData.PortObj.TTL_input_counter+1;

    src.UserData.LabelTTLInputLapseTime.Text=sprintf('%3.3f',seconds(laps_time));
end;
