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

load handel;
player=audioplayer(y,Fs);
play(player,[1 (get(player,'SampleRate')*3)]);

