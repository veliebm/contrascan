% Creates a MOVING moving window average for each subject.
%
% Driver code created by Benjamin Velie on 7/27/2021.
% veliebm@ufl.edu
%
% stead2singtrialsCont and regressionMAT created by Andreas Keil in May 2021.
% akeil@ufl.edu
% --------------------------------------------------------------------------

parameters_file = 'processed/movingmovingwindoweeg/parameters.json';
do_all(parameters_file)
delete_lock_file(mfilename('fullpath'))

function do_all(parameters_file)
    % Moving moving window all subjects whose metadata we've stored in a JSON file.
    all_parameters = read_json(parameters_file);

    for i = 1:numel(all_parameters)
        parameters = all_parameters(i)
        main(parameters.in_filename, parameters.in_dir, parameters.out_stem, parameters.out_tsv_name);
    end
end

function main(in_filename, in_dir, out_stem, out_tsv_name)
    % Run stead2singtrialsCont on a subject.
    eeglab;
    stead2singtrialsCont(in_filename, in_dir, 0, 1:1000, 1:1000, 12, 600, 500, 1, out_stem);
    convert_to_tsv(strcat(out_stem, '.amp.at'), out_tsv_name)
end
function [data] = read_json(in_path)
    % Read a JSON file.
    fname = in_path; 
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    data = jsondecode(str);
end
function convert_to_tsv(in_filename, out_filename)
    % Convert an emegs file to tsv.
    amplitude = ReadAvgFile(in_filename);
    dlmwrite(out_filename, amplitude, '\t');
end

function [winmat, fftamp] = stead2singtrialsCont(EEGsetfile, EEGsetfile_directory, plotflag, bslvec, ssvepvec, foi, sampnew, SampRate, saveflag, outname)
    % samprate is the NEW NEW sample rate  - user needs to calculate such that
    % integer number of sample points fits in a cycle  foi = freq of interest
    % this is for 2-D files from eeglab

    EEG = pop_loadset(EEGsetfile, EEGsetfile_directory);
     
    % determine settings for each moving window 
    % translate ssvepvec into new sample rate. the sample vector is indicated in old sample rate above 
    ssvepvec = floor(ssvepvec .* (sampnew./SampRate)); 
    % this to be done outside the loop to save time, needed for winshift proc
    sampcycle=1000/sampnew; %added code for the new samplerate
    tempvec = round((1000/foi)/sampcycle); % this makes sure that average win duration is exactly @@, which is the duration in sp of one cyle at @@ Hz = @@ ms, sampled at 250 Hz
    longvec = repmat(tempvec,1,100); % repeat this many times, at least for duration of entire epoch, subsegments are created later 
    winshiftvec_long = [ssvepvec(1) cumsum(longvec)+ ssvepvec(1)]; % use cumsum function to create growing vector of indices for start of the winshift
    tempindexvec = find(winshiftvec_long > ssvepvec(end)); 
    
    endindex = tempindexvec(1);  % this is the first index for which the winshiftvector exceeds the data segment 
    winshiftvec = winshiftvec_long(1:endindex-5);
   
    % need this stuff for the spectrum
    shiftcycle=round(tempvec*4);
    samp=1000/sampnew;
           freqres = 1000/(shiftcycle*samp); %added code to find the appropriate bin for the frequency of interest for each segment
           freqbins = 0:freqres:(sampnew/2); 
           min_diff_vec = freqbins-foi;  %revised part
           min_diff_vec = abs(min_diff_vec); %revised
           targetbin=find(min_diff_vec==min(min_diff_vec)); %revised 
           
   % now, determine the onsets of the segments of the 2-D file that will
   % replace the trials...
   
   % starting with finding the first S2 trigger
   searchindex = 1;
   while searchindex 
       if strcmp(EEG.event(searchindex).type, 'S  2')
       FirstS2index = searchindex; searchindex = 0; 
       else
       searchindex = searchindex + 1;
       end
   end
   
   %  finding the final S2 trigger
   
   for searchindex = 1:size(EEG.event, 2)
           if strcmp(EEG.event(searchindex).type, 'S  2')
           LastS2index = searchindex;
           end
   end
   
   %this is the 2-D data starting with the first S2 event, plus 1 second
   %baseline, ending after onset of last S2 event, plus 10 seconds
   mat = EEG.data(1:31, EEG.event(FirstS2index).latency-500:EEG.event(LastS2index).latency+5000); 
   
   % finally determine the onsets and duration of the overlapping 2-D segments that will stand in as
   % "trials"
   
   duration = 2000;
   
   segmentonsets = 1:1000:size(mat,2)-(duration+1); 
   
   % the above is still with the OLD sample rate, upsampling is done below,
   % whithin segments
   
   for fileindex = 1 : 1
      
       fftamp = []; 
       SNR = []; 
       phasestabmat = []; 
       fftcomplexphasevec = []; 
       
       NTrials = size(mat,3); 
   
     for trial = 1:length(segmentonsets)
           
       Data = mat(:,segmentonsets(trial):segmentonsets(trial)+(duration-1));
          
          fouriersum = []; 
      
       %============================================================
       % resample data
       %===========================================================   
       
           Data=double(Data');
               
           resampled=resample(Data,sampnew,SampRate);           
               
           Data = resampled';  
       
       %============================================================
       % 2. Baseline correction
       %===========================================================
       
       disp ('subtracting baseline')
       
       datamat = bslcorr(Data, bslvec);
   
   %===========================================================
       % 3. moving window procedure with 4 cycles  !!!
       % 
       %===========================================================
       disp('moving window procedure')
       
       winmatsum = zeros(size(datamat,1),shiftcycle); %4 cycles
       
        if plotflag, h = figure; end
      
        for winshiftstep = 1:length(winshiftvec)
           
           winmatsum = (winmatsum + regressionMAT(datamat(:,[winshiftvec(winshiftstep):winshiftvec(winshiftstep)+(shiftcycle-1)]))); % time domain averaging for win file
           fouriermat = fft(regressionMAT(datamat(:,[winshiftvec(winshiftstep):winshiftvec(winshiftstep)+(shiftcycle-1)]))');
           
           %for within trial phase locking
           fouriermat = fft(regressionMAT(datamat(:,[winshiftvec(winshiftstep):winshiftvec(winshiftstep)+(shiftcycle-1)]))');
           fouriercomp = fouriermat(targetbin,:)'; 
           
           if winshiftstep ==1
               fouriersum = fouriercomp./abs(fouriercomp);
           else
               fouriersum = fouriersum + fouriercomp./abs(fouriercomp);
           end
           
          if plotflag
              subplot(2,1,1), plot(1:sampcycle:shiftcycle*sampcycle, regressionMAT(datamat(:,[winshiftvec(winshiftstep):winshiftvec(winshiftstep)+(shiftcycle-1)]))'), title(['sliding window starting at ' num2str((winshiftvec(winshiftstep))*sampcycle)  ' ms ']), xlabel('time in milliseconds')
              subplot(2,1,2), plot(1:sampcycle:shiftcycle*sampcycle, winmatsum'), title(['sum of sliding windows; number of shifts:' num2str(winshiftstep) ]), ylabel('microvolts')
             %    subplot(3,1,3), hold on, circle([0,0],1,200,'-');
             %    plot([0;(imag(fouriercomp(120)./tenHZampfft(120)))], [0;(real(fouriercomp(120)./tenHZampfft(120)))]);title('phase angle of window') 
              pause(.4)
          end
           
     %    movmat(index) = getframe(h)
     
       end
       
       winmat = winmatsum./length(winshiftvec);
   
       %===========================================================
       % 5. determine amplitude and Phase using fft of the winmat (i.e. the
       % average
       %===========================================================
   
       % FFT of the average winmat
       % 
       
       NFFT = shiftcycle-1; % one cycle in sp of the desired frequency times 4 oscillations (-1)
       NumUniquePts = ceil((NFFT+1)/2); 
       fftMat = fft (winmat', (shiftcycle-1));  % transpose: channels as columns (fft columnwise)
       Mag = abs(fftMat);                                                   % Amplitude berechnen
       Mag = Mag*2;   
       
       Mag(1) = Mag(1)/2;                                                    % DC trat aber nicht doppelt auf
       if ~rem(NFFT,2)                                                       % Nyquist Frequenz (falls vorhanden) auch nicht doppelt
           Mag(length(Mag))=Mag(length(Mag))/2;
       end
       
       Mag=Mag/NFFT;                                                         % FFT so skalieren, dass sie keine Funktion von NFFT ist
       
       fftamp = [fftamp Mag((targetbin),:)'];
       
       SNR = [SNR log10(Mag((targetbin),:)'./ mean(Mag(([targetbin-2 targetbin+2]),:))').*10];
       
       % phase stability
   
       phasestabmat = [phasestabmat abs(fouriersum./winshiftstep)]; 
           
      end % trials
      
      outmat.fftamp = fftamp; 
      outmat.phasestabmat = phasestabmat; 
      outmat.winmat = winmat; 
      
      eval(['save ' outname '.fftamp.mat fftamp -mat'])
      
      if saveflag == 1 
          SaveAvgFile([outname '.amp.at'], fftamp)
      elseif saveflag == 2 
          SaveAvgFile([outname '.phastab.at'], phasestabmat)
       elseif saveflag == 3
            SaveAvgFile([outname '.twin.at'], winmat)
      elseif saveflag == 4
            SaveAvgFile([outname '.SNRdB.at'], SNR)
      end
          
   
      fclose('all');
   
    end % files
end
function[OutMat] = regressionMAT(InMat)
    % regression
    % compute linear regression and subtracts regression from input vector
    
    OutMat = zeros(size(InMat));
    
    % make sure InMat has electrodes as rows
        
    dummy = size(InMat);

    X = [ones(1,dummy(2));1:dummy(2)]';
    
    
    for channel = 1 : size(InMat,1)
        
        regkoefs = pinv(X' * X) * (X' * InMat(channel, :)');
        
        
        reggerade = [1:dummy(2)].*regkoefs(2,1) + regkoefs(1,1);
        OutMat(channel, :) = InMat(channel, :) - reggerade;
    end
end
