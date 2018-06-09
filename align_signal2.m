% align_signals2   Takes 2 structs for OptiTrack (OT) and OnBoard (OB) and
% specifies the method used (depends on the onboard data structure). It
% crops, aligns the data sets based on the phi or theta state and uses the
% same translations on the other state. This function is specifically built
% for data alignment for the Bebop 2 quadrotor and the OptiTrack as used by
% TU Delft C&S.
%
% INPUT: 
%   OT        OptiTrack data structure  
%   OB        OnBoard data structure
%   method    switch case specifying the state which shall be used

% OUTPUT: 
%   OT_a                Shifted OptiTrack data structure
%   OB_a                Shifted OnBoard data structure
%   s1                  Signal used for OptiTrack
%   s2                  Signal used for OnBoard
%   shifting_vector     the shifting vector used to shift other states

% By Leon Sijbers 2017

function [OT_shift,OB_shift,s1,s2,shifting_vector] = align_signal2(OT_a,OB_a,method)

% In order to show the operations the script performs
show_animation  = 0;
show_crop       = 0;
show_while_loop = 0;
% To cross check another state to verify the shift went okay
VERIFICATION    = 1;
% Optionally output the .mat file with shifted data structs
saving          = 0;

% Set maxlag in order to indicate the maximum shift required 
% Depending on the data set it might be give different results
% Try tweaking with the maxlag to get better results

% If the whole data set is slightly misaligned: alter percent_firstrun as
% this aligns the first # percent of the data

% if parts of the data set is slightly misaligned: alter maxlag_otherloops
switch method
    case 'phi_ot'
        s1 = OB_a.phi_ot*57.3;
        s2 = OT_a.PHI;
        maxlag_firstloop = 500;
        maxlag_otherloops = 300;
        percent_firstrun = 0.15;
    case 'theta_ot'
        s1 = OB_a.theta_ot*57.3;
        s2 = OT_a.THETA;
        maxlag_firstloop = 300;
        maxlag_otherloops = 300;
        percent_firstrun = 0.15;
    case 'phi'
        s1 = OB_a.phi*57.3;
        s2 = OT_a.PHI;
        maxlag_firstloop = 3;
        maxlag_otherloops = 1;
        percent_firstrun = 0.15;
    case 'theta'
        s1 = OB_a.theta*57.3;
        s2 = OT_a.THETA;
        maxlag_firstloop = 3;
        maxlag_otherloops = 1;
        percent_firstrun = 0.15;
end
        
% append shortest of two signals with zeros to get equal size
if length(s1) > length(s2) 
s2(numel(s1))=0;
else 
s1(numel(s2))=0;
end

% Crop out bad data points i.e. the leading and trailing zeros
[crop_signal_1,crop_signal_2,crop_vector] = CropSignal(s1,s2,0.2,100,show_crop);

if (show_crop == 1) % Plot if required
figure; hold on; plot(s1); plot(s2);
h=line([crop_vector(1,1); crop_vector(1,1)],[80;-80]); set(h,'Color','black');
h=line([crop_vector(2,1); crop_vector(2,1)],[80;-80]); set(h,'Color','black');
hold off; title('Original signals with cutoff lines')
end
% Use cropped signals from now on
cs1 = crop_signal_1;
cs2 = crop_signal_2;

shifting_vector = zeros(2,1);

if (show_crop == 1) % Plot if required
figure; hold on; plot(cs1); plot(cs2); 
hold off; title('Cropped signal')
end

%% Next part
% First make sure the first #% of the samples is aligned well
K = round(percent_firstrun*size(cs1,1));
% % K = 5000;

% The first time shift (for the entire sample)

[ds1,ds2,~,~,~,LAGS_max] = TimeShift2(cs1,cs2,maxlag_firstloop,K,show_animation);

% First shift
shifting_vector(:,1) = [1; LAGS_max];

[~,sawtooth_K] = Outliers(ds1,ds2,1000,5,20);

% Find first element in the signal that has more than K elements larger
% than L*std.dev.
first_nonzero_element = 1; %find(sawtooth_K,1,'first');
% sawtooth_reset = ones(size(sawtooth_K,1),1);
mini = 1; maxi = size(sawtooth_K,1);
counting = 0;
show_animation = 0;
while (sum(sawtooth_K) > 1) && (counting < 100)
counting = counting + 1;
% Locate where to shift signal
[sawtooth_signal,sawtooth_K] = Outliers(ds1,ds2,150,3.5,20);
% And eliminate the last shift such that it isnt repeated
sawtooth_K(1:(first_nonzero_element+sawtooth_signal(first_nonzero_element)+750)) = 0;
% Find first detected outlier
first_nonzero_element = find(sawtooth_K,1,'first');
if isempty(first_nonzero_element)
    
   break 
end

% Determine the lead/lag shift for the next part
K = min(50,(size(sawtooth_K,1)-first_nonzero_element));
% Calculate shift for first K parts after the first detected outlier
a1 = ds1(first_nonzero_element:(first_nonzero_element+K));
a2 = ds2(first_nonzero_element:(first_nonzero_element+K));
% And find the shift required as LAGS_max
[~,~,~,~,~,LAGS_max] = TimeShift2(a1,a2,maxlag_otherloops,K,show_animation);

% Insert the lead/lag shift into the cropped signal ds1
if ~(LAGS_max == 0)
% First entry already taken so next element 
shifting_vector(:,(counting+1)) = [first_nonzero_element; LAGS_max];

ds1 = [ds1(1:first_nonzero_element-1); ds1(first_nonzero_element-1)*ones(LAGS_max*-1,1); ds1(first_nonzero_element:end)];
ds2(numel(ds1)) = 0;
end
if show_while_loop == 1
    delete(gco)
    figure(3); set(gcf, 'Position', [800, 50, 700, 650]);   plot(sawtooth_signal(mini:maxi),'c');hold on; plot(ds2(mini:maxi),'k'); plot(ds1(mini:maxi),'r'); plot(sawtooth_K(mini:maxi)*10,'y'); 
    hold off;
    axis([0 inf -25 250 ]);
    txt2 = sprintf('%d',round(counting));
    text(1000,-20,txt2);
% axis([first_nonzero_element-1000 first_nonzero_element+4000 -80 10 ])
end
% pause(1/5)
end

% Script to check shifted signal vs other state in the same dataset
if (VERIFICATION == 1)
VERIFICATION_script;
end

%% Translation of other states

clear OT_shift OB_shift
[~, ~, OT_shift, OB_shift] = AlignStates(OT_a,OB_a,shifting_vector,crop_vector);

if (saving == 1)
filename = sprintf('Shifted_data_set_index_%d.mat',index);
save(['D:\001 - Studie\01 - Master Thesis\02_Data\OJF_november\_shifted_data_OJF_november_2017\' filename],'OT_shift','OB_shift','OT_a','OB_a','WIND','PARA','take')
end
end

function [crop_signal_1,crop_signal_2,crop_vector] = CropSignal(signal_1,signal_2,std_limit,window_size,show_crop)
% CropSignal    Takes two signal and calculates the effective range of the data based on
%   the zeros in the leading or trailing part of the set.
%   Data is cropped to exclude the parts where both derivatives of signal 1
%   and 2 are zero
%   based on the original signal; find the first good data point and exclude
%   data before this element
%
% INPUT: 
%   signal_1        first data set 
%   signal_2        second data set
%   std_limit       signal is cropped if movstd <  0.2*std_limit
%   window_size     moving std.dev window size 
%   show_crop       boolean to generate the plots
% OUTPUT: 
%   crop_signal_1   Cropped signal 1 
%   crop_signal_2   Cropped signal 2 with size(crop_signal_1)
%   crop_vector     Contains first and last element to be included of data set
% 
% By Leon Sijbers 2017
if ~exist('show_crop','var')
    show_crop = 0;
end
if ~exist('std_limit','var')
    std_limit = 0.2;
end
if ~exist('window_size','var')
    window_size = 100;
end

% set_axis = [size(signal_1,1)-5000 size(signal_1,1) -inf inf];
set_axis = [0 inf -inf inf];
% figure; axis(set_axis); hold on;refline([0 std(diff(signal_1))]); plot(abs(signal_1 - mean(signal_1)),'k'); plot(statement_1,'r'); plot(signal_1); hold off;

% Check if the signal 1 or 2 has deviations that are insignificant 
% Calculate the std deviation of the data segment based on -0.5*window_size
% and 0.5*window_size and compare it with some low value that is arbitrary
% for now.
statement_1 = ~(movstd(signal_1,window_size) < 0.2*std_limit);
statement_2 =  ~(movstd(signal_2,window_size) < 0.2*std_limit);

statement_max = and(statement_1,statement_2);
figure; axis(set_axis); hold on; plot(statement_1,'g'); plot(statement_2,'r'); plot(statement_max); plot(signal_1); plot(signal_2); hold off;
cutoff_vector = ( statement_max == 1 );


% Find first nonzero element starting from first element
% Add window_size/2 as the std.dev also based on some future elements
first_good_datapoint = find(cutoff_vector,1,'first')+window_size/2-1;

% And crop the signal based on the cutoff_vector
crop_signal_1 = signal_1(first_good_datapoint:end);
crop_signal_2 = signal_2(first_good_datapoint:end);

% based on the resized signal; find the first bad data point and exclude
% data after this element

% Find first nonzero element starting from last element
first_bad_datapoint = find(cutoff_vector((first_good_datapoint+1):end),1,'last')-window_size/2;
crop_signal_1 = crop_signal_1(1:first_bad_datapoint);
crop_signal_2 = crop_signal_2(1:first_bad_datapoint);

crop_vector   = [first_good_datapoint; first_good_datapoint+first_bad_datapoint];
if (show_crop == 1)
    leading_zeros = 10*ones(first_good_datapoint,1);
    trailing_zeros = 10*ones( (size(signal_1,1)-first_bad_datapoint),1);
    show_crop_signal_1 = [leading_zeros+1; crop_signal_1; trailing_zeros+1];
    show_crop_signal_2 = [leading_zeros; crop_signal_2; trailing_zeros];
    figure(5); set(gcf, 'Position', [800, 50, 700, 350]);  %axis([0 inf -10 1]);
    hold on; 
    plot(signal_1+1); plot(signal_2+1);  plot(show_crop_signal_1); plot(show_crop_signal_2);
    hold off; 
    title('Show cropping')
    legend('signal_1','signal_2','crop_1','crop_2');
end

end
function [y] = ConsequentNumber(x)
% ConsequentNumber Builds a vector that contains a cumulative of consequetive ones such that
% a step size becomes a ramp
% This creates 'spikes' with height as the width of the step block
 
% By Leon Sijbers 2017
n = size(x,1);
f = find([true;not(x)~=0;true]);
y = zeros(1,n);
y(f(1:end-1)) = diff(f);
y = cumsum(y(1:n))-(1:n);
y = y';
end
function [sawtooth_signal,sawtooth_K] = Outliers(signal_1,signal_2,K,L,M)
% Outliers  Calculates the 'outliers' of the root mean square for two signals and
% calculates their significance based on the amount of consequent outliers
%
% INPUT: 
%   signal_1 and signal_2 a (n x 1) vector
%  K, L and M scalar integers: 
%   K is samples to calculate std.dev
%   L to calculate the amount of std.devs to indicate an outlier
%   M the amount of consequent samples required for a reshift of signal
%
% OUTPUT: 
%   sawtooth_signal vector that contains cumsum of occurances of outliers
%   sawtooth_K      logical vector where sawtooth_signal exceeds M 
%
% By Leon Sijbers 2017

% Take the difference
diff_ds = abs(signal_2-signal_1);

% Which elements are larger than L times the standard deviation
deviation = diff_ds > (L * std(diff_ds(1:K)));

% Get the number of repeated outliers
% Called sawtooth as this produced a sawtooth-like signal for each repeated
% segment
[sawtooth_signal] = ConsequentNumber(deviation);

% See which of the repeated outliers are more than M times in a row
sawtooth_K = sawtooth_signal > M;

% figure(8); hold on; plot(signal_1); plot(signal_2); plot(sawtooth_signal); refline([0 (L * std(diff_ds(1:K)))]); plot(diff_ds); hold off;
end
function [shifted_signal_1,shifted_signal_2,C,LAGS,C_max,LAGS_max] = TimeShift2(signal_1,signal_2,MAXLAG,K,show_animation)
% TimeShift2   Calculates the delay between signal_1 and signal_2 for 
%   the first K samples based on a least mean squared error convolution of 
%   range -MAXLAG:MAXLAG (Based on the xcorr()-function which is
%   unavailable without the econometrics toolbox and also less robust to y-axis offsets)
% 
% INPUT: 
%   signal_1        nx1 vector that is kept stationary
%   signal_2        nx1 vector that is shifted 
%   K               scaler that denotes the 1:K part of signal 1 and 2
%   MAXLAG          scalar that denotes the limit range over which to shift
%   show_animation  boolean to allow the animation to generate

% OUTPUT:
%   shifted_signal_1    vector of a shifted signal 1
%   shifted_signal_2    vector of a shifted signal 1
%   C                   vector containing indices for the lags
%   LAGS                Corresponding lags/delays
%   C_max               Index correponding to the maximum delay
%   LAGS_max            delay required for minimal least squared error
% 
% By Leon Sijbers, 2018
%
% See also XCORR
if ~exist('show_animation','var')
    show_animation = 0;
end
signal_3 = signal_1(1:K);
signal_4 = signal_2(1:K);
C1 = ones(MAXLAG,1).*2;
C2 = ones(MAXLAG,1).*2;
LAGS1 = zeros(MAXLAG,1);
LAGS2 = zeros(MAXLAG,1);

for i = 1:MAXLAG
% Calculate squared differences for each element
dif_sig = (signal_3((1):(end-i))  - signal_4((1):(end-i))).^ 2 / max((signal_3((1):(end-i))  - signal_4((1):(end-i))).^ 2) ;
% count the lagged element number
LAGS1(i) = i;
% Find mean squared dif
C1(i) = mean(dif_sig); 
if show_animation == 1
delete(gco)
figure(6); set(gcf, 'Position', [50, 50, 1500, 0700]);  axis([0 inf -10 1]);
plot(signal_3,'k'); hold on; plot(signal_4); plot(dif_sig*10); plot(C1,'r'); hold off;
pause(1/100)
end
signal_3 = [signal_1( (i+1):K); zeros((i),1)];
end
signal_3 = signal_1(1:K);
signal_4 = signal_2(1:K);
for i = 1:MAXLAG
% Calculate squared differences for each element
dif_sig = (signal_3((1+i):(end))  - signal_4((1+i):(end))).^ 2 / max( (signal_3((1+i):(end))  - signal_4((1+i):(end))).^ 2);
% count the lagged element number
LAGS2(i) = -i;
% Find mean squared dif
C2(i) = mean(dif_sig); 
if show_animation == 1
    dif_sig = [zeros(i,1); dif_sig];
delete(gco)
figure(3); set(gcf, 'Position', [50, 50, 1500, 0700]);  axis([0 inf -10 1]);
plot(signal_3,'k'); hold on; plot(signal_4); plot(dif_sig*10); hold off;
pause(1/100)
end
signal_3 = [zeros((i),1); signal_1( (1):(K-i))];
end
% Find least mean squared value and the corresponding lag value
if min(C1) > min(C2)
  [C_max,I] = min(C2);
LAGS_max = LAGS2(I);  
else
   [C_max,I] = min(C1);
LAGS_max = LAGS1(I); 
end
LAGS = [LAGS1; LAGS2];
C = [C1; C2];

if LAGS_max > 0
ds1 = signal_1(abs(LAGS_max):end);
ds1(numel(signal_2))=0;
ds2 = signal_2;
elseif LAGS_max < 0 
ds2 = signal_2(abs(LAGS_max):end);
ds2(numel(signal_1))=0;
ds1 = signal_1;
else 
    ds2 = signal_2;
    ds1 = signal_1;
end

shifted_signal_1 = ds1;
shifted_signal_2 = ds2; % UPDATE

end
function [OT_a, OB_a, OT_shift, OB_shift] = AlignStates(OT_a,OB_a,shifting_vector, crop_vector)
% ALIGNSTATES   Shifts vector OB_a such that it aligns with OT_a
%   First, OT_a and OB_a are cropped acording to crop_vector 
%   Second, OB_a is shifted using the shifting_vector
%   Third, OT_a size is matched to OB_a size by appending zeros
%   
% INPUT: 
%   OT_a the structure of states from Optitrack (~48 fields)
%   OB_a the structure of states logged onboard Bebop (~47 fields)
%   note: OB_a and OB_t indices do not typically correspond to same states
% 
% OUTPUT: 
%   OT_a        Original Optitrack structure
%   OB_a        Original Onboard structure
%   OT_shift    Shifted Optitrack structure
%   OB_shift    Shifted Onboard structure
% 
% By Leon Sijbers, 2018

OT_shift = OT_a; 
OB_shift = OB_a;
OT_names = fieldnames(OT_shift);
OB_names = fieldnames(OB_shift);
a2 = OT_shift.(OT_names{2});

for j = 1:size(OB_names,1)
% Try first for Phi states
b = OB_shift.(OB_names{j});
a = a2((crop_vector(1,1):crop_vector(2,1)));
% Cropping all states of both the Optitrack as Onboard logging according to
% the cropping vector created by function: CropSignal()
% crop_vector could be made to manually crop by inserting 
% crop_vector = [start_element; end_element]
b = b(crop_vector(1,1):crop_vector(2,1));

% Do first shift to align first part
b = [b(1:shifting_vector(1,1)); b(shifting_vector(1,1)+shifting_vector(2,1):end)];
% Append size of the shifted Onboard signal with zeros to match Optitrack 
b(numel(a)) = 0;
if size(shifting_vector,2) > 1 
n = length(shifting_vector);
% Loop over the next n-1 shifts in the vector
for i = 2:n 
    element = shifting_vector(1,i);
    delay = shifting_vector(2,i);  
    if ~(delay == 0)
        b = [b(1:element-1); b(element-1)*ones(delay*-1,1); b(element:end)];
    end
end


OB_shift.(OB_names{j}) = b;    
end
end
for i = 1:size(OT_names,1)
a = OT_shift.(OT_names{i});
a = a(crop_vector(1,1):crop_vector(2,1),:);
% Adjust size of Optitrack state to match Onboard state
   a(numel(b),:) = 0;
OT_shift.(OT_names{i}) = a;   
end

end

