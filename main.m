%%%%%%%%%%%%% main.m file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      Basic DFT-based Frequency Analysis:
%
%
% Code flow: 
%      1.  Load input image 'checker.gif'.
%      2.  Plot the necessary functions before and after modulating and using the
%      logarithmic scale
%      3.  Set DC component to zero
%      4. Plot the necessary functions after modulating and using the
%      logarithmic scale
%      5. Create a gaussian filter and pass it through the input image and plot appropriate figures
%      6. Perform the above step after zero padding the input image

%  The following functions are called:
%       lpfilter
%       dftuv
%       
%  Author:      Pavan Gurudath
%  Date:        10/21/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1a: 
% Clear out all memory and close all open MATLAB windows
close all; clear all; clc;
imtool close all;
warning off;
delete *.jpeg;

% Read image
f=imread('checker.gif');
[M, N] = size(f);

hFigure = imtool(double(f));
set(hFigure,'NumberTitle','off','Name','input checker image');

f=double(f);
fpad=f;                         %Store temporary f(x,y) which is used later    

F   = fft2(f,M,N);              %FFT of f(x,y)

F_temp = F;                     %Store temporary F(u,v) which is used later

%Display and store |F(u,v)|
hFigure = imtool(abs(F));    
set(hFigure,'NumberTitle','off','Name','|F(u,v)|');
modF = getimage(hFigure);
% modF = modF./(max(max(modF))).*255;
imwrite(uint8(modF),'mod[F(u,v)].jpeg');

%Display and store log(1+|F(u,v)|)
hFigure = imtool(log(1+abs(F)),[]);     
set(hFigure,'NumberTitle','off','Name','log(1+|F(u,v)|)');
logmodF= getimage(hFigure);
logmodF=logmodF./(max(max(logmodF))).*255;
imwrite(uint8(logmodF),'log(1+mod[F(u,v)]).jpeg');


%Modulate f(x,y) and compute F(u,v)
f_modulation = f;
for i=1:M
    for j=1:N
        f_modulation(i,j) = f(i,j).*((-1)^(i+j));
    end
end

F_modulation = fft2(f_modulation); 
    
%Display and store |Fm(u,v)|
hFigure = imtool(abs(F_modulation),[]);    
set(hFigure,'NumberTitle','off','Name','|Fm(u,v)|');
modFm = getimage(hFigure);
modFm = modFm.*255./(max(max(modFm)));
imwrite(uint8(modFm),'mod(Fm(u,v)).jpeg');

%Display and store log(1+|Fm(u,v)|)
hFigure = imtool(log(1+abs(F_modulation)),[]);
set(hFigure,'NumberTitle','off','Name','log(1+|Fm(u,v)|)');
logmodFm = getimage(hFigure);
logmodFm = logmodFm.*255./(max(max(logmodFm)));
imwrite(uint8(logmodFm),'log(1+mod[Fm(u,v)]).jpeg');

F(1,1) =0.0;        %Set F(0,0) to zero
%Display and store the new image g(x,y)
f_newImage = real(ifft2(F));    %Inverse dft after F(0,0) is set to 0

hFigure = imtool(f_newImage);
set(hFigure,'NumberTitle','off','Name','g(x,y)');
g_newImage = getimage(hFigure);
for i=1:M
    for j=1:N
        if g_newImage(i,j) <0
            g_newImage(i,j) =0;
        else
            g_newImage(i,j) =g_newImage(i,j)*255;
        end
    end
end

imwrite(uint8(g_newImage),'g(x,y)_prob1.jpeg');

%% PART 1b: 
%Obtaining the Gaussian filter
sig = 15;           %Value of sigma
figure(1);
H  = (lpfilter('gaussian', M, N, sig));
mesh1=mesh(H(1:5:256, 1:5:256));      %3D plot of H(u,v)
saveas(mesh1,'mesh1.png');
%Compute the convolution of f(x,y) and h(x,y) 
 G = double(zeros(M,N));
F = F_temp;
G = F.*H;                       %Lowpass gaussian filtered image in frequency domain 


%Display and store log(1+|G(u,v)|) without modulation
hFigure = imtool(log(1+abs(G)));
set(hFigure,'NumberTitle','off','Name','log(1+|G(u,v)|) w/o modulation');
Glpf = getimage(hFigure);
Glpf = Glpf.*255/(max(max(Glpf)));
imwrite(Glpf,'log(1+mod[G(u,v)]_wo_mod.jpeg');


g= real(ifft2(G));              %Lowpass gaussian filtered image in spatial domain
%Display and store g(x,y) 
hFigure = imtool(g);
set(hFigure,'NumberTitle','off','Name','g(x,y) w/o modulation');
glpf = getimage(hFigure);
glpf = glpf.*255./(max(max(glpf)));
imwrite(uint8(glpf),'g(x,y)_wo_mod.jpeg');


%Modulation of g(x,y)
g_modulation = g;
for i=1:M
    for j=1:N
        g_modulation(i,j) = g(i,j).*((-1)^(i+j));
    end
end

%Display and store g_modulation(x,y) 
hFigure = imtool(abs(g_modulation));
set(hFigure,'NumberTitle','off','Name','g_modulation(x,y)');
gmlpf = getimage(hFigure);
gmlpf = gmlpf.*255./(max(max(gmlpf)));
imwrite(uint8(gmlpf),'g_modulation(x,y).jpeg');


%Display and store G_modulation(u,v) appropriately
G_modulation = fft2(g_modulation);

hFigure = imtool(log(1+abs(G_modulation)),[]);
set(hFigure,'NumberTitle','off','Name','G_modulation(u,v)');
Gmlpf = getimage(hFigure);
Gmlpf = Gmlpf.*255./(max(max(Gmlpf)));
imwrite(uint8(Gmlpf),'log(1+mod[Gmodulation(u,v)]).jpeg');


%% Part 1c:

% Zero Padding
P=2*M;
Q=2*N;

% Gaussian LPF of size 512x512
sig=15;
Hpad   = lpfilter('gaussian', P, Q, sig);

figure(2);
mesh2=mesh(Hpad(1:7:P, 1:7:Q));               %3D plot of H(u,v)
saveas(mesh2,'Mesh2.png');
fpad = double(fpad);        %Recalling the input image f
Fpad=fft2(fpad,P,Q);        %FFT of input image with zero padding

%Gaussian low pass filtered image in frequency domain for zero padded
%images

Gpad = double(zeros(P,Q));  
for i=1:P
    for j=1:Q
        Gpad(i,j) = Fpad(i,j)*Hpad(i,j);
    end
end

% imtool(log(1+abs(Gpad)),[]);        Not displayed because g(x,y) is not
%                                     modulated

gpad= real(ifft2(Gpad));            %IFFT of G(u,v)
%Modulation of g(x,y) with zero padding
for i=1:P
    for j=1:Q
        gpad_modulation(i,j) = gpad(i,j)*((-1).^(i+j));
    end
end

%Display and store |gpad_modulation(x,y)|
hFigure = imtool(abs(gpad_modulation));
set(hFigure,'NumberTitle','off','Name','gpad_modulation(x,y)');
gpadm = getimage(hFigure);
gpadm = gpadm.*255./(max(max(gpadm)));
imwrite(uint8(gpadm),'gpad_modulation(x,y).jpeg');

%Display and store log(1+|Gpad_modulation(u,v)|)
Gpad_modulation = fft2(gpad_modulation, P,Q);
    
hFigure = imtool(log(1+abs(Gpad_modulation)),[]);
set(hFigure,'NumberTitle','off','Name','Gpad_modulation(u,v)');
Gpadm = getimage(hFigure);
Gpadm = Gpadm.*255./(max(max(Gpadm)));
imwrite(uint8(Gpadm),'Gpad_modulation(u,v).jpeg');

%% PART 2
%   Please run part_two.m for the second question
part_two();