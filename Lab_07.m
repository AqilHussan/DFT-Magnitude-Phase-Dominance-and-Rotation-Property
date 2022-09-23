close all;
Fourier = imread("fourier.png");                                 %%Image Read
figure('Name','Input:fourier');
imshow(Fourier);                                                 %%Image Display
Fourier_transform = imread("fourier_transform.png");             %%Image Read
figure('Name','Input:fourier_transform');                        %%Image Display
imshow(Fourier_transform);
peppers = imread("peppers_small.png");                           %%Image Read
figure('Name','Input:peppers_small');                                  %%Image Display
imshow(peppers);  
dftFourier=findDFT(Fourier);                                     %%Finding the 2D Dft 
FourierMag=abs(dftFourier);                                      %%Taking the magnitude of the 2D Dft 
FourierPhase=angle(dftFourier);                                  %%Taking the Phase the 2D Dft 
dftFourier_transform=findDFT(Fourier_transform);                 %%Finding the 2D Dft 
Fourier_transformMag=abs(dftFourier_transform);                  %%Taking the magnitude of the 2D Dft
Fourier_transformPhase=angle(dftFourier_transform);              %%Taking the Phase the 2D Dft 
figure('Name','Fourier_transform:fourier');
imshow(log(abs(fftshift(dftFourier))), []);                      %%Displaying the 2d fft
figure('Name','Fourier_transform:fourier_transform');
imshow(log(abs(fftshift(dftFourier_transform))), []);            %%Displaying the 2d fft
DFTImage3=FourierMag.*exp(1j*Fourier_transformPhase);            %%F3(k, l) = |F1(k, l)|exp(jφ2(k,l))
DFTImage4=Fourier_transformMag.*exp(1j*FourierPhase);            %%F4(k, l) = |F2(k, l)|e(jφ1(k,l))
Image3=findIDFT(DFTImage3);                                      %%Finding IDFT 
figure('Name','Image3');
imshow(uint8(Image3));                                           %%Display IDFT
Image4=findIDFT(DFTImage4);                                      %%Finding IDFT
figure('Name','Image4');
imshow(uint8(Image4));                                           %%Display IDFT
RotatedDFT=complex(zeros(size(peppers)));                        %%Initializing Rotated DFT 
RotatedPeppers=(zeros(size(peppers)));                           %%Initializing Rotated Image
[rows,colms]=size(peppers);                                      %%Finding the size of Image
theta=pi/2;                                                      %%Angle of rotation
R=[cos(theta),-sin(theta);sin(theta),cos(theta)];                %%Defining the rotation matrix
Comp=double(peppers);                                            %%Copying image as double
centre_height=floor(((rows)/2));                                 %%Finding the center of image
for k=1:rows
    for l=1:colms
        for m=1:rows
            for n=1:colms
                %%%Finding Rotated dft with origin is at the center of the image.
                k1=k-centre_height  ;                 
                l1=l-centre_height ;
                P=[m n]*R*([k1 l1]');
                RotatedDFT(k,l)=RotatedDFT(k,l)+((Comp(m,n))*exp(-1i*2*pi*P/(rows)));
            end
        end
    end
end
%%%Rotating the input image by same angle
for i=1:rows
    for j=1:colms
        x=rows-1-i-centre_height ;                   %Co-ordinates of pixel with respect to the centre 
        y=colms-1-j-centre_height; 
        [new_x,new_y]=rotate(theta,x,y) ;              %Applying the transform       
        new_x=centre_height-new_x  ;                   %Change new_x and new_y with respect to the new centre   
        new_y=centre_height-new_y;
        X=floor(new_x);                                %Find X' Y'
        Y=floor(new_y);
        a=new_x-X;                                     %Finding the deviation from actual point and X' Y'
        b=new_y-Y;
        %Applying bilinear interpolation
        if(1<X+1 && X+1<rows+1 && 1<Y+1 && Y+1<colms+1)      
            RotatedPeppers(i,j)=(1-a)*(1-b)*peppers(X,Y)+(1-a)*(b)*peppers(X,Y+1)+(a)*(1-b)*peppers(X+1,Y)+(a)*(b)*peppers(X+1,Y+1);
            
        end
    end
end
%Display Rotated DFT
figure('Name','Rotated Fourier_transform:peppers_small');
imshow(log(abs((RotatedDFT))), []);
ImageRotated=findIDFT(RotatedDFT);                     %Finding the IDFT of rotated dft
figure('Name','Inverse Fourier_transform:Rotated');
imshow(uint8(abs(ImageRotated)));                      %%Display  IDFT
figure('Name','Rotated peppers_small');
imshow(uint8(RotatedPeppers));                         %%Display the rotated image



%%%Function to find  2D DFT using row-column decomposition.
function DFT=findDFT(Image)
[rows,colms]=size(Image);               %Finding the rows and column size
Intermediate=zeros(size(Image));        %Initializing the Intermediate matrix 
DFT=zeros(size(Image));                 %Initializing the DFT matrix 
%Taking Rows and finding DFT of each row
for i=1:rows
    Intermediate(i,:)=fft(Image(i,:));
end
%Taking Columns of intermediate and finding DFT of each column 
for i=1:colms
    DFT(:,i)=fft(Intermediate(:,i));
end
end
%%%Function to find  2D IDFT using row-column decomposition.
function IDFT=findIDFT(Image)
[rows,colms]=size(Image);                 %Finding the rows and column size
Intermediate=zeros(size(Image));          %Initializing the Intermediate matrix
IDFT=zeros(size(Image));                  %Initializing the IDFT matrix
%Taking Rows and finding IDFT of each row
for i=1:rows
    Intermediate(i,:)=ifft(Image(i,:));
end
%Taking Columns of intermediate and finding IDFT of each column 
for i=1:colms
    IDFT(:,i)=ifft(Intermediate(:,i));
end
end

function [newX,newY]=rotate(angle,x,y)
    

%      |newX|   |cos(angle)   sine(angle)| | x |
%      |    | = |                        | |   |
%      |newY|   |-sine(angle) cos(angle) | | y |
%      


   
    cosine=cos(angle);                  %Finding the cosine               
    sine=sin(angle);                    %Finding the sine        
    newX=x*cosine-y*sine ;              %Applying the transform                  
    newY=x*sine+y*cosine;
   
end


    
    