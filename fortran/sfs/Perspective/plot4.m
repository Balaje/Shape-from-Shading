%%
% Script to generate the grayscale image

N = 400;
A = imread('thing1.png');
gray = mat2gray(imresize(A,[N,N]));

GRAY = gray(:,:,1);

for i=1:N
    for j=1:N
        if(GRAY(i,j) < 0.1)
            GRAY(i,j) = 0.1;
        end
    end
end

imshow(GRAY)
dlmwrite('peaks.txt',GRAY);