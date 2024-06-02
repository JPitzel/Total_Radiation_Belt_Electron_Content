function download_data(PSD_directory,flux_directory,epoch_start,epoch_stop)

oldFolder = cd(flux_directory);

% Fetch MagEIS data files
sats=["rbspa","rbspb"];
for sat_i=1:2
    for epoch_t=epoch_start:epoch_stop
        [Y,M,D]=datevec(epoch_t);
        yyyy=sprintf('%04d',Y);
        mm=sprintf('%02d',M);
        dd=sprintf('%02d',D);
        url="https://rbsp-ect.newmexicoconsortium.org/data_pub/"+sats(sat_i)+"/mageis/level3/pitchangle/"+yyyy+"/";
        filename=sats(sat_i)+"_rel04_ect-mageis-L3_"+yyyy+mm+dd+"_v8.1.0.cdf";
        filename811=sats(sat_i)+"_rel04_ect-mageis-L3_"+yyyy+mm+dd+"_v8.1.1.cdf";
        if isfile(filename) || isfile(filename811)
            %disp(filename+" exists in the local directory.")
        else
            try
                disp("Fetching "+filename)
                websave(filename, url+filename);
            catch
                disp("Fetching "+filename811)
                websave(filename811, url+filename811);
            end
        end
    end
end

cd(oldFolder);
cd(PSD_directory);

% Fetch PSD data files
sats=["rbspa","rbspb"];
for sat_i=1:2
    for epoch_t=epoch_start:epoch_stop
        d=datevec(epoch_t);
        [Y,M,D]=datevec(epoch_t);
        yyyy=sprintf('%04d',Y);
        mm=sprintf('%02d',M);
        dd=sprintf('%02d',D);
        doy=epoch_t - datenum(Y,1,0);
        doy=sprintf('%03d',doy);
        t_secs = (epoch_t - 719529)*86400; % seconds since start epoch 1970-01-01
        url="https://rbspgway.jhuapl.edu/rTools/psdN/lib/php/getCDF.php?cli=%20-T%20";
        url=url+num2str(t_secs)+"%20-E%2086399%20-S%20"+sats(sat_i)+"%20-A%20TS04";
        filename="PSD_"+sats(sat_i)+"_mageis_"+yyyy+doy+".cdf";
        if isfile(filename)
            disp(filename+" exists in the local directory.")
        else
            disp("Fetching "+filename)
            disp(url)
            websave('tmp.tar.gz', url);
            try
                files=gunzip('tmp.tar.gz');
                untar(files{1});
                dir="PSD_"+yyyy+doy+"_"+yyyy+doy+"_"+sats(sat_i)+"_TS04";
                movefile(dir+"/"+filename,".")
                delete(dir+"/*")
                rmdir(dir)
            catch
                disp("ERROR: gzip file is corrupt.")
            end
        end
    end
end
warning('off')
delete("tmp.tar.gz")
delete("tmp.tar")
warning('on')
cd(oldFolder);

end
