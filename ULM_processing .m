%%
clear all
close all
clc

addpath('C:\Users\vinco\Documents\MATLAB');
addpath('C:\Users\vinco\Documents\Magistrale\Ultrasound Technologies\UTMA Project');

%% Import of the RFB data and set of parameters

RFB_DATA_DIRECTORY = 'D:\UNIVERSITA\IE\Ultrasound Technologies for Medical Applications\Acquisitions\A3';
RFB_DATA_FILE = '\A3_1.rfb';
Nlines=129; % number of lines in the image
PRF = 4000; % Hz

% sampling frequency
fs=50e6;
% propagating medium S.O.S.
c0=1540;
% transmit offset in m
OffSet=0.005;
% Pitch in mm
Pitch=0.245;
LastGate = 2944;

% Reading the data
DataObj=DataUlaopPostBeamforming(strcat(RFB_DATA_DIRECTORY, RFB_DATA_FILE)); 
Read(DataObj,'firstPri',1,'npri',GetTotalPri(DataObj)); %Read all the avaialble PRIs
[nGate,nPri]=size(DataObj.LastReadData); %extract the number of gate and the number of pris
time=DataObj.LastReadTime(1)+(0:nGate-1)/DataObj.fs; %time axis
rf=DataObj.LastReadData;

% Total data
figure
imagesc(20*log10(abs(rf(1:LastGate, :)))); colorbar;


%% From first raw frame, compute the actual depth and axial resolution
firstColumn = 56;
firstFrame = rf(1:LastGate, firstColumn:firstColumn+Nlines-1);
LateralDimension=[1:Nlines]*Pitch;
LateralDimension=LateralDimension-(Pitch*floor(Nlines/2));
Depth=OffSet+[1:LastGate]*c0*1/(2*fs);
f = figure;
set(gca,'Color', [242, 244, 245]./255);
imagesc(LateralDimension, Depth*1000, 20*log10(abs(firstFrame)./max(max(abs(firstFrame))))); colorbar;
ylabel('Depth [mm]'); xlabel('Lateral Dimension [mm]');
caxis([-50 0])
%colormap('gray');
pbaspect([30 50 1]);


%% Build a 3D matrix with fast direction 
columnIndex = firstColumn;   
nFrames = floor((nPri-columnIndex)/ Nlines); % Find frames number -  flooring
matrix3D = ones(LastGate, Nlines, nFrames);  % Inizialization

for frameIndex = 1:nFrames
        matrix3D(:, :, frameIndex) =  rf(1:LastGate, columnIndex:columnIndex+Nlines-1); 
        columnIndex = columnIndex + Nlines;
end

%showVideoFrom3DMatrix(matrix3D, Nlines, LateralDimension, Depth, PRF, fs)

%% Build the filter
tx_freq = 8e6
cycles = 4
first_value = tx_freq - tx_freq/cycles % Min cut-off frequency
second_value = tx_freq + tx_freq/cycles % Max cut-off frequency
% filtering
[A,B] = butter(8,[first_value second_value]/25e6);
%h = fvtool(A,B)
w = linspace(0, fs/2, LastGate/2);
H = freqz(A,B, LastGate/2);
plot(w, 20*log10(abs(H)));

%% Filtering and envelope detection on 3DMatrix 
filteredMatrix3D = ones(LastGate, Nlines, nFrames);  % Inizialization
for frameIndex = 1:nFrames
    for i = 1:Nlines
        filteredsignal = filter(A, B, matrix3D(:, i, frameIndex)); 
        filteredMatrix3D(:, i, frameIndex) = abs(hilbert(filteredsignal));
    end
end

%showVideoFrom3DMatrix(filteredMatrix3D, Nlines, LateralDimension, Depth, PRF, fs, "A14");

%% Background subtraction
noBackMatrix3D = ones(LastGate, Nlines, nFrames-1);
for frameIndex = 1:nFrames-1
    noBackMatrix3D(:, :, frameIndex) = filteredMatrix3D(:, :, frameIndex+1) - filteredMatrix3D(:, :, 1);
end

%showVideoFrom3DMatrix(noBackMatrix3D, Nlines, LateralDimension, Depth, PRF, fs, "A7");
%createVideoFrom3DMatrix(noBackMatrix3D, Nlines, LateralDimension, Depth, PRF, "A7");

%% FILTERING: Single Value Decomposition filtering
% From filteredMatrix3D build the Casorati Matrix (coloum=1 frame, row=pixel over time)
casorati_matrix = reshape(filteredMatrix3D, [numel(filteredMatrix3D(:,:,1)), size(filteredMatrix3D, 3)]);

% Plot casorati matrix
%imagesc(20*log10(abs(casorati_matrix)./max(max(abs(casorati_matrix))))); colorbar;

% Apply SVD on Casorati matrix
[U,S,V] = svd(casorati_matrix, 'econ');

% Plot singular values curve
plot([1:nFrames], diag(S));


S_filtered = S;
clutter_threshold = 6;
noise_threshold = 80;

S_filtered(1:clutter_threshold, 1:clutter_threshold) = 0;
S_filtered(noise_threshold:end, noise_threshold:end) = 0;

casorati_filtered = U * S_filtered * V';
svdFilteredMatrix3D = reshape(casorati_filtered, size(filteredMatrix3D));

showVideoFrom3DMatrix(svdFilteredMatrix3D, Nlines, LateralDimension, Depth, PRF, fs)

%% DETECTION: cross-correlation - template vs frame
% Suppongo che questa elaborazione debba essere fatta per ogni frame
detection_frame = svdFilteredMatrix3D(:, :, 114);


%1. Identifico il template
f = figure;
set(gca,'Color', [242, 244, 245]./255);
imagesc(20*log10(abs(detection_frame)./max(max(abs(detection_frame))))); colorbar;
ylabel('Depth [mm]'); xlabel('Lateral Dimension [mm]');
caxis([-50 0])
%colormap('gray');
pbaspect([30 50 1]);
hold on;

% Disegna un rettangolo per selezionare la ROI
h = drawrectangle;  % Disegna il rettangolo interattivo

% 3. Ottieni le coordinate e le dimensioni della ROI
roi_position = round(h.Position);  % Ottieni posizione e dimensioni come [x, y, width, height]

% 4. Estrai la ROI dall'immagine come template
x = roi_position(1);
y = roi_position(2);
width = roi_position(3);
height = roi_position(4);
template = detection_frame(y:y+height-1, x:x+width-1);  % Ritaglia la ROI dall'immagine

% 5. Visualizza il template per confermare
f = figure;
set(gca,'Color', [242, 244, 245]./255);
imagesc(20*log10(abs(template)./max(max(abs(detection_frame))))); colorbar;
ylabel('Depth [mm]'); xlabel('Lateral Dimension [mm]');
caxis([-50 0])

% 6. Normalizziamo il template in modo da osservare solo la struttura,
% senza preoccuparci delle differenze di intensità
template_normalized = (template - mean(template(:))) / std(template(:));


%Normalizziamo anche il frame per cui facciamo la cross-correlation
detection_frame_normalized = (detection_frame - mean(detection_frame(:))) / std(detection_frame(:));

%% Cross correlation
J = normxcorr2(template_normalized, detection_frame_normalized); 

f = figure;
set(gca,'Color', [242, 244, 245]./255);
imagesc(J); colorbar;
%caxis([0.2 1]);


%% DETECTION: max intensity
%% Find local maxima
local_maxima_matrix = ones(size(noBackMatrix3D));

%  local_maxima_matrix is binary
for f = 1:size(noBackMatrix3D,3)
    local_maxima_matrix(:, :, f) = imregionalmax(noBackMatrix3D(:,:,f));
end

% take the original values
for f = 1:size(noBackMatrix3D,3)
    local_maxima_matrix(:, :, f) = abs(noBackMatrix3D(:,:,f) .* local_maxima_matrix(:, :, f));
end

showVideoFrom3DMatrix(local_maxima_matrix, Nlines, LateralDimension, Depth, PRF, fs, "A3");

%% Apply a threshold
max_intensity_matrix = ones(size(local_maxima_matrix));

threshold = 0.05;
for f = 1:size(local_maxima_matrix, 3)
    max_intensity_frame = ones(size(normalized_frame, 1), size(normalized_frame, 2));
    matrix_frame = local_maxima_matrix(:, :, f);
    matrix_frame = abs(matrix_frame)./max(max(abs(matrix_frame)));
    for i = 1:size(normalized_frame, 1)   
        row = matrix_frame(i, :); 
        for c = 1:Nlines
            if row(c) < threshold
               max_intensity_frame(i, c) = 0; 
            else
              max_intensity_frame(i, c) = row(c); 
            end 
        end
    end
    max_intensity_matrix(:, :, f) = max_intensity_frame;
end


showVideoFrom3DMatrix(max_intensity_matrix, Nlines, LateralDimension, Depth, PRF, fs, "A3");

%% Have arrays for each frame with local maxima
max_row = 0;
for f = 1: size(max_intensity_matrix, 3)
  [rows cols] = find(max_intensity_matrix(:,:,f)); 
  max_row = max(max_row, size(rows,1));
end

local_maxima_arrays = zeros(max_row, 2, size(max_intensity_matrix, 3));

for f = 1: size(max_intensity_matrix, 3)
  [rows cols] = find(max_intensity_matrix(:,:,f)); 
  local_maxima_arrays(1:size(rows,1), 1, f) = rows;
  local_maxima_arrays(1:size(cols,1), 2, f) = cols;
end


%% old 
max_intensity_matrix = ones(size(noBackMatrix3D));

dectectionFrame = noBackMatrix3D(:, :, 50);
normalized_frame = abs(dectectionFrame)./max(max(abs(dectectionFrame)));
threasholded_frame = ones(size(normalized_frame, 1), size(normalized_frame, 2));

%row = normalized_frame(950, :);
%plot(row);

threshold = 0.1;
for f = 1:size(noBackMatrix3D, 3)
    max_intensity_frame = ones(size(normalized_frame, 1), size(normalized_frame, 2));
    matrix_frame = noBackMatrix3D(:, :, f);
    matrix_frame = abs(matrix_frame)./max(max(abs(matrix_frame)));
    for i = 1:size(normalized_frame, 1)   
        row = matrix_frame(i, :); 
        for c = 1:Nlines
            if row(c) < threshold
               max_intensity_frame(i, c) = 0; 
            else
              max_intensity_frame(i, c) = row(c); 
            end 
        end
    end
    max_intensity_matrix(:, :, f) = max_intensity_frame;
end

showVideoFrom3DMatrix(max_intensity_matrix, Nlines, LateralDimension, Depth, PRF, fs, "A3");


for f = 1:size(max_intensity_matrix,3)
    
    max_intensity_matrix(:, :, f) = imgaussfilt(max_intensity_matrix(:,:,f));
end
%showVideoFrom3DMatrix(max_intensity_matrix, Nlines, LateralDimension, Depth, PRF, fs, "A3");

%%
f = figure;
set(gca,'Color', [242, 244, 245]./255);
imagesc(20*log10(normalized_frame)); colorbar;
ylabel('Depth [mm]'); xlabel('Lateral Dimension [mm]');
%caxis([-50 0])


f = figure;
set(gca,'Color', [242, 244, 245]./255);
imagesc(20*log10(threasholded_frame)); colorbar;
ylabel('Depth [mm]'); xlabel('Lateral Dimension [mm]');
%caxis([-50 0])

%%
%J = medfilt2(dectectionFrame);
%f = figure;
%set(gca,'Color', [242, 244, 245]./255);
%imagesc(20*log10(abs(J)./max(max(abs(J))))); colorbar;
%ylabel('Depth [mm]'); xlabel('Lateral Dimension [mm]');
%caxis([-50 0])

local_maxima_matrix = ones(size(max_intensity_matrix));

for f = 1:size(max_intensity_matrix,3)
    local_maxima_matrix(:, :, f) = imregionalmax(max_intensity_matrix(:,:,f));
end

showVideoFrom3DMatrix(local_maxima_matrix, Nlines, LateralDimension, Depth, PRF, fs, "A3");


%ROW 1258
row = abs(noBackMatrix3D(:, 74, 50));
%plot(row);


%% TRACKING

frame1 = size(find(local_maxima_arrays(:,1,1)),1);
frame2 = size(find(local_maxima_arrays(:,1,2)),1);

cost_matrix = zeros(frame1, frame2);

% 1. Build cost matrix
for f = 1:size(max_intensity_matrix,3)
    %To finish
    cost_matrix = sqrt(sum((local_maxima_array(1:frame1, :, 2) - local_maxima_array(1:frame1, :, 1) .^ 2));
end




