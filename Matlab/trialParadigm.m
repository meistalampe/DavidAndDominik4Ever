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
clear all;
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

% make filler sequence of 10 seconds
interval = zeros(1,10*Fs);

% repetitions of every stimulus
trial_per_stim = 200;

% calculate the length of all stimuli for all trials including pauses
stim_length = trial_per_stim*length(interval); 
% how many stimuli are presented inbetween breaks
stim_interval_count = 40;
% calculate how many intervals are needed
number_interval = ceil(trial_per_stim * number_of_stim)/stim_interval_count;
number_break = number_interval;
% set break to 1 min
interval_break = zeros(1,6*length(interval));
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
% trigger values( a = 1, u = 3, i=5)
audio_clip = zeros(2,clip_length);
trigger_signal = zeros(1,clip_length);


loop_condition = true;

% counter and position variables
position = 1;
stim_count = 0;
a_count = 0;
u_count = 0;
i_count = 0;

interval_count = 0;

while loop_condition == true
 
% random stimulus    
stim = randi(number_of_stim);

switch stim
    case 1
        if(a_count <= trial_per_stim)
            audio_clip(1,position:(position+length(interval)-1)) = stimuli(1,:);
            audio_clip(2,position:(position+length(interval)-1)) = stimuli(2,:);
            trigger_signal(1,position:(position+length(a)-1)) = 1;

            position = position + length(interval);
            a_count = a_count + 1;
            stim_count = stim_count + 1;
        end
        
    case 2
        if(u_count <= trial_per_stim)
            audio_clip(1,position:(position+length(interval)-1)) = stimuli(3,:);
            audio_clip(2,position:(position+length(interval)-1)) = stimuli(4,:);
            trigger_signal(1,position:(position+length(u)-1)) = 3;

            position = position + length(interval);
            u_count = u_count + 1;
            stim_count = stim_count + 1;
        end
    case 3
        if(i_count <= trial_per_stim)
            audio_clip(1,position:(position+length(interval)-1)) = stimuli(5,:);
            audio_clip(2,position:(position+length(interval)-1)) = stimuli(6,:);
            trigger_signal(1,position:(position+length(i)-1)) = 5;

            position = position + length(interval);
            i_count = i_count + 1;
            stim_count = stim_count + 1;
        end
end


if(stim_count == stim_interval_count)
    audio_clip(1,position:(position+length(pause)-1)) = 0;
    audio_clip(2,position:(position+length(pause)-1)) = 0;
    trigger_signal(1,position:(position+length(pause)-1)) = 0;
    
    position = position + length(pause);
    stim_count = 0;
    interval_count = interval_count + 1;
end

if(interval_count == number_interval)
    loop_condition = false;
end

end

% write audio file
filename = 'clip.wav';
audiowrite(filename,audio_clip(1,:),aFs);

