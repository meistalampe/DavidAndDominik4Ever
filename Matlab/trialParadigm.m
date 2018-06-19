% %% file header
% 
% % filename:     trialParadigm
% % author:       dominik limbach
% % date:         13.06.18
% 
% % description:
% %   - read stimuli files (audioread)
% %   - process stimuli, same length per vowel (a,i,u)
% %   - create a null sequence, for waiting periods (1second)
% %   - write a loop that randomly patches an audio file 
% %   - audio file: random stimuli, 10x null sequence, repeat 30 times (=
% %   5min)...60x null sequence, repeat all until every stimulus was
% %   presented 200 times 
% %   - implement subroutine that counts stimuli and eliminates a stimulus
% %   type as soon as 200 is reached ( while loop with break condition which
% %   becomes true if one variable is at 200,...., second and third loop that
% %   present remaining variables
% %   - create new audiofile(audiowrite)
% %   - create extra row for trigger which hold values according to the vowel
% %   while it is playing
% %   - create a time delay of 10% of the intervall time that is randomly
% %   weaven in and counters low frequency overlapping with the signal
% %   - normalize audiosignal to a amplitude of 1 and adjust loudness(dB)
% %   according to a sample with randomized phase (fft, randi phase,ifft -
% %   detect loudness)
%% file header

% filename:     trialParadigm
% author:       dominik limbach
% date:         13.06.18

% description:
%   - read stimuli files (audioread)
%   - process stimuli, same length per vowel (a,i,u)
%   - create a null sequence, for waiting periods (1second)
%   - write a loop that randomly patches an audio file 
%   - audio file: random stimuli, 10x null sequence, repeat 30 times (=
%   5min)...60x null sequence, repeat all until every stimulus was
%   presented 200 times 
%   - implement subroutine that counts stimuli and eliminates a stimulus
%   type as soon as 200 is reached ( while loop with break condition which
%   becomes true if one variable is at 200,...., second and third loop that
%   present remaining variables
%   - create new audiofile(audiowrite)
%   - create extra row for trigger which hold values according to the vowel
%   while it is playing


%% Clear
clc;
clear;
close all;

%% load stimuli files

[a,aFs] = audioread('a.wav');
[u,uFs] = audioread('u.wav');
[i,iFs] = audioread('i.wav');

% number of different stimuli
number_of_stim = 3;

% play command to test files
% player=audioplayer(a,aFs);
% play(player,[1 (get(player,'SampleRate')*3)]);

%% create filler sequence

% set sampling frequency if not uniform
% Fs = 44100;
Fs = aFs;

% make filler sequence of int_size seconds
int_size = 5;
interval = zeros(1,int_size*Fs);

% repetitions of every stimulus
trials_per_stim = 200;

% calculate the length of all stimuli for all trials including pauses
stim_length = trials_per_stim*length(interval); 
% how many stimuli are presented inbetween breaks
stim_per_interval = 30;
% calculate how many intervals are needed
number_intervals = ceil(trials_per_stim * number_of_stim)/stim_per_interval;
number_break = number_intervals-1;
% set break to 1 min
break_size = 30 / int_size;
interval_break = zeros(1,break_size*length(interval));
% set final length of the clip including breaks
clip_length = number_of_stim*stim_length + number_break * length(interval_break);

%% audio clip composition

% compose stimuli intervals
stimuli = zeros(number_of_stim*2, length(interval)); % stereo signal so 2 rows per stim
% a
stimuli(1,1:length(a)) = a(:,1);
stimuli(2,1:length(a)) = a(:,2);
% u
stimuli(3,1:length(u)) = u(:,1);
stimuli(4,1:length(u)) = u(:,2);
% i
stimuli(5,1:length(i)) = i(:,1);
stimuli(6,1:length(i)) = i(:,2);

% initiate audio clip and trigger signal
% 2 rows for audio signal 1 row for trigger
% trigger values( a = 2, u = 4, i=6, pause=1)
audio_clip_L = zeros(1,clip_length);
audio_clip_R = zeros(1,clip_length);
trigger_signal = zeros(1,clip_length);


loop_condition = true;

% counter and position variables
position = 1;
stim_count = 0;
count_a = 0;
count_u = 0;
count_i = 0;

interval_count = 0;

% create delay parameter
% set boundaries for delay
minimum = (0.05*length(interval));
maximum = (0.15*length(interval));
% create a vector with random delay times for each stimulus
N = (trials_per_stim*number_of_stim)+1;
delay_time = randi([minimum maximum],1,N);
% add all delays
delay_sum = sum(delay_time);
% create empty vector to compensate for delay
delay = zeros(1,delay_sum);

delay_count = 1;

% preallocation
audio_clip_L = [audio_clip_L delay];
audio_clip_R = [audio_clip_R delay];
trigger_signal = [trigger_signal delay];

while loop_condition == true
 
% random stimulus    
stim = randi(number_of_stim);
position = position + delay_time(delay_count);


switch stim
    case 1
        if(count_a <= trials_per_stim)
            audio_clip_L(1,position:(position+length(interval)-1)) = stimuli(1,:);
            audio_clip_R(1,position:(position+length(interval)-1)) = stimuli(2,:);
            %trigger_signal(1,position:(position+length(a)-1)) = 1;
            trigger_signal(1,position:(position + Fs*0.02)) =1;

            position = position + length(interval);
            count_a = count_a + 1;
            stim_count = stim_count + 1;
            delay_count = delay_count + 1;
        end
        
    case 2
        if(count_u <= trials_per_stim)
            audio_clip_L(1,position:(position+length(interval)-1)) = stimuli(3,:);
            audio_clip_R(1,position:(position+length(interval)-1)) = stimuli(4,:);
            %trigger_signal(1,position:(position+length(u)-1)) = 3;
            trigger_signal(1,position:(position + Fs*0.04)) =1;

            position = position + length(interval);
            count_u = count_u + 1;
            stim_count = stim_count + 1;
            delay_count = delay_count + 1;
        end
    case 3
        if(count_i <= trials_per_stim)
            audio_clip_L(1,position:(position+length(interval)-1)) = stimuli(5,:);
            audio_clip_R(1,position:(position+length(interval)-1)) = stimuli(6,:);
            %trigger_signal(1,position:(position+length(i)-1)) = 5;
            trigger_signal(1,position:(position + Fs*0.06)) =1;

            position = position + length(interval);
            count_i = count_i + 1;
            stim_count = stim_count + 1;
            delay_count = delay_count + 1;
        end
end


if(stim_count == stim_per_interval)
    audio_clip_L(1,position:(position+length(interval_break)-1)) = 0;
    audio_clip_R(1,position:(position+length(interval_break)-1)) = 0;
    trigger_signal(1,position:(position+length(interval_break)-1)) = 0;
    trigger_signal(1,position:(position + Fs*0.1)) = 1;
    
    position = position + length(interval_break);
    stim_count = 0;
    interval_count = interval_count + 1;
end

if(interval_count == number_intervals)
    loop_condition = false;
    
end

end

% normalize amplitude
audio_clip_L = audio_clip_L ./ max(abs(audio_clip_L));
audio_clip_R = audio_clip_R ./ max(abs(audio_clip_R));
% 
% % fft
% fft_L = fft(audio_clip_L);
% amp = 2*abs(fft_L)/length(fft_L);
% phs = angle(fft_L);
% 
% for j = length(phs)
% phs(j) = (pi+pi).*rand(1,1) + -pi;
% end
% ifft_L = zeros(1,length(fft_L));
%  
% for i = length(fft_L)
%    real = amp(i);
%    imag = phs(i);
%    ifft_L(i) = real + imag;
%     
% end
% 
% ifft_L = ifft(ifft_L);
% 
% audiowrite('test.wav',ifft_L(1,:),aFs);
% % 
% % save 
% save('vowelStimulus.mat','-v7.3');
% 
% % write audio file
% filename = 'clip.wav';
% audiowrite(filename,audio_clip(1,:),aFs);
% 
% 
% % trigger signal 
% 
% filenamet = 'trigger.wav';
% audiowrite(filenamet,trigger_signal(1,:),aFs);

time = 1:(length(trigger_signal));
time = time./aFs;

figure;
hold on;
plot(time(158000000:end),trigger_signal(158000000:end));
plot(time(158000000:end),audio_clip_L(158000000:end));
ylabel('amplitude');
xlabel('time [s]');
legend('trigger signal','audio signal')



% 
% (158000000:end)

% %% Clear
% clc;
% clear;
% close all;
% 
% %% load stimuli files
% 
% [a,aFs] = audioread('a.wav');
% [u,uFs] = audioread('u.wav');
% [i,iFs] = audioread('i.wav');
% 
% % number of different stimuli
% number_of_stim = 3;
% 
% % play command to test files
% % player=audioplayer(a,aFs);
% % play(player,[1 (get(player,'SampleRate')*3)]);
% 
% %% create filler sequence
% 
% % set sampling frequency if not uniform
% % Fs = 44100;
% Fs = aFs;
% 
% % make filler sequence of 10 seconds
% interval = zeros(1,5*Fs);
% 
% % repetitions of every stimulus
% trial_per_stim = 300;
% 
% % calculate the length of all stimuli for all trials including pauses
% stim_length = trial_per_stim*length(interval); 
% % how many stimuli are presented inbetween breaks
% stim_interval_count = 30;
% % calculate how many intervals are needed
% number_interval = ceil(trial_per_stim * number_of_stim)/stim_interval_count;
% number_break = number_interval;
% % set break to 1 min
% interval_break = zeros(1,6*length(interval));
% % set final length of the clip including breaks
% clip_length = number_of_stim*stim_length + number_break * length(interval_break);
% 
% %% audio clip composition
% 
% % compose stimuli intervals
% stimuli = zeros(number_of_stim*2, length(interval)); % stereo signal so 2 rows per stim
% % a
% stimuli(1,1:length(a)) = a(:,1);
% stimuli(2,1:length(a)) = a(:,2);
% % u
% stimuli(3,1:length(u)) = u(:,1);
% stimuli(4,1:length(u)) = u(:,2);
% % i
% stimuli(5,1:length(i)) = i(:,1);
% stimuli(6,1:length(i)) = i(:,2);
% 
% % initiate audio clip and trigger signal
% % 2 rows for audio signal 1 row for trigger
% % trigger values( a = 1, u = 3, i=5)
% audio_clip = zeros(2,clip_length);
% trigger_signal = zeros(1,clip_length);
% 
% 
% loop_condition = true;
% 
% % counter and position variables
% position = 1;
% stim_count = 0;
% count_a = 0;
% count_u = 0;
% count_i = 0;
% 
% interval_count = 0;
% 
% while loop_condition == true
%  
% % random stimulus    
% stim = randi(number_of_stim);
% 
% switch stim
%     case 1
%         if(count_a <= trial_per_stim)
%             audio_clip(1,position:(position+length(interval)-1)) = stimuli(1,:);
%             audio_clip(2,position:(position+length(interval)-1)) = stimuli(2,:);
%             %trigger_signal(1,position:(position+length(a)-1)) = 1;
%             trigger_signal(1,position) = 1;
% 
%             position = position + length(interval);
%             count_a = count_a + 1;
%             stim_count = stim_count + 1;
%         end
%         
%     case 2
%         if(count_u <= trial_per_stim)
%             audio_clip(1,position:(position+length(interval)-1)) = stimuli(3,:);
%             audio_clip(2,position:(position+length(interval)-1)) = stimuli(4,:);
%             %trigger_signal(1,position:(position+length(u)-1)) = 3;
%             trigger_signal(1,position) = 3;
% 
%             position = position + length(interval);
%             count_u = count_u + 1;
%             stim_count = stim_count + 1;
%         end
%     case 3
%         if(count_i <= trial_per_stim)
%             audio_clip(1,position:(position+length(interval)-1)) = stimuli(5,:);
%             audio_clip(2,position:(position+length(interval)-1)) = stimuli(6,:);
%             %trigger_signal(1,position:(position+length(i)-1)) = 5;
%             trigger_signal(1,position) = 5;
% 
%             position = position + length(interval);
%             count_i = count_i + 1;
%             stim_count = stim_count + 1;
%         end
% end
% 
% 
% if(stim_count == stim_interval_count)
%     audio_clip(1,position:(position+length(pause)-1)) = 0;
%     audio_clip(2,position:(position+length(pause)-1)) = 0;
%     trigger_signal(1,position:(position+length(pause)-1)) = 0;
%     
%     position = position + length(pause);
%     stim_count = 0;
%     interval_count = interval_count + 1;
% end
% 
% if(interval_count == number_interval)
%     loop_condition = false;
% end
% 
% end
% 
% % write audio file
% filename = 'clip.wav';
% audiowrite(filename,audio_clip(1,:),aFs);
% 
% % trigger signal 
% time = 1:length(trigger_signal);
% time = time./aFs;
% 
% figure;
% plot(time,trigger_signal);
% label