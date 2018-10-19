%%%%%%%%%%%%% part_two.m file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      Filtering a corrupted image
%
%
% Code flow: 
%      1.  Load input image 'lake.tiff'.
%      2.  Plot the necessary functions before and after modulating and using the
%      logarithmic scale
%      3.  Add noise to the input image
%      4. Design a notch filter
%      5. Perform spatial convolution and observe the output image and
%      compare it with the original image. 
% 
%       
%  Author:      Pavan Gurudath
%  Date:        10/21/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2

% Clear out all memory and close all open MATLAB windows
close all; clear all; clc;
imtool close all;
warning off;
delete *.jpeg;

% Read image lake.tif 
lake = imread('lake.tif');

% Obtain the first layer since it contains the necessary information in the .tif image
f = lake(:,:,1); 
[rows,col] = size(f);  

% Display and store f(x,y)
hFigure = imtool(f);
set(hFigure,'NumberTitle','off','Name','f(x,y)');
lake_new = getimage(hFigure);
imwrite(uint8(lake_new),'f(x,y).jpeg');


%Generate noisy lake image
     
c = double(zeros(rows,col));                   %the corrupted version of the original image


for i = 1:rows
    for j=1:col
        c(i,j) = double(f(i,j))+(32).*cos((2*pi*32*i)/512);
    end
end    

%Display and store c(x,y)
hFigure = imtool(c,[]);    
set(hFigure,'NumberTitle','off','Name','lake.tif with noise');
lake_corrupted = getimage(hFigure);
imwrite(uint8(lake_corrupted),'lake_corrupted.jpeg');

%Display and store |C(u,v)|
C = fft2(c,512,512);

    % Modulation of input image
    c_modulated = double(zeros(rows,col));
    c = double(c);
    for i = 1:rows
        for j = 1:col
            c_modulated(i,j) = c(i,j)*((-1)^(i+j));
        end
    end
    
Cm = fft2(c_modulated,rows,col);
hFigure = imtool(log(1+abs(Cm)),[]);    
set(hFigure,'NumberTitle','off','Name','|C(u,v)|');
lake_noise = getimage(hFigure);
lake_noise = lake_noise.*255./(max(max(lake_noise))); 
imwrite(uint8(lake_noise),'[C(u,v)].jpeg');    


% Modulation of input image
f_modulated = double(zeros(rows,col));
f= double(f);

for i = 1:rows
    for j = 1:col
        f_modulated(i,j) = f(i,j)*((-1)^(i+j));
    end
end

%Display and store |F(u,v)|
F_modulated = fft2(double(f_modulated));
hFigure = imtool(log(1+abs(F_modulated)),[]);
set(hFigure,'NumberTitle','off','Name','log(1+abs[F(u,v)])');
Fm = getimage(hFigure);
Fm = Fm.*255./(max(max(Fm)));
imwrite(uint8(Fm),'log(1+abs[F(u,v)]).jpeg');


% design of the notch filter
H = double(ones(rows,col));
H(33,1) = 0;
H(481,1) = 0;
hFigure = imtool(H,[]);
set(hFigure,'NumberTitle','off','Name',' |H(u,v)| ');
Mag_NotchFilter = getimage(hFigure);
Mag_NotchFilter = Mag_NotchFilter.*255./(max(max(Mag_NotchFilter)));
imwrite(uint8(Mag_NotchFilter),'mod[H(u,v)].jpeg');



% filter the noisy image with the notch filter
G = double(zeros(rows,col));
for i=1:rows
    for j=1:col
        G(i,j)=H(i,j).*C(i,j);
    end
end

% the output image after filtering
g = real(ifft2(G));

hFigure = imtool(g,[]);
set(hFigure,'NumberTitle','off','Name','Filtered Output image');
lake_filtered = getimage(hFigure);
imwrite(uint8(lake_filtered),'g(x,y).jpeg');

g_modulated = double(zeros(rows,col));
for i = 1:rows
    for j = 1:col
        g_modulated(i,j) = g(i,j)*((-1)^(i+j));
    end
end

G_modulated = fft2(g_modulated,512,512);
hFigure = imtool(log(1+abs(G_modulated)),[]);
set(hFigure,'NumberTitle','off','Name','|G(u,v)|');
Mag_res_lake_filtered = getimage(hFigure);
Mag_res_lake_filtered = Mag_res_lake_filtered.*255./(max(max(Mag_res_lake_filtered)));
imwrite(uint8(Mag_res_lake_filtered),'mod[G(u,v)].jpeg');


% difference between the original image and output image
diffImage = double(zeros(rows,col));

for i = 1:rows
    for j = 1:col
        diffImage(i,j) = f(i,j)-g(i,j);
    end
end

hFigure = imtool(diffImage);
set(hFigure,'NumberTitle','off','Name','Difference blw original and output image');
difference = getimage(hFigure);
diffImage = diffImage.*255./(max(max(diffImage))-min(min(diffImage)));
imwrite(uint8(diffImage),'diffImage(x,y).jpeg');

d_modulated = double(zeros(rows,col));
for i = 1:rows
    for j = 1:col
        d_modulated(i,j) = (diffImage(i,j))*((-1)^(i+j));
    end
end

D_modulated = fft2(double(d_modulated));
hFigure = imtool(log(1+abs(D_modulated)),[]);
set(hFigure,'NumberTitle','off','Name','DFT magnitude response of difference b/w original and output image');
Mag_res_difference = getimage(hFigure);
imwrite(Mag_res_difference,'log(1+mod[diffImage(u,v)]).jpeg');