% Read in image
I2 = imread("ece6560_img4.jpg");
figure(1)
image(I2)
% find dimensions
[m, n, o] = size(I2);
% add noise to image if necessary, else uncomment line 9
 J = imnoise(I2,'gaussian');
%  J = I2;
% Output noisy image
 figure(2)
 image(J)

J = double(J);
% determine step size for x and y direction
dx = 1/(m-1);
dy = 1/(n-1);

% set number of iterations and time step
numt = 120;
dt = 6e-14;
% Set constant K
K = 4.7e7;
% Run the PDE implementation
for l=1:numt
       for k =1:3
       for j = 3:m-2
           for i=3:n-2
                 I_mag_grad = abs((J(j,i+1,k)-2*J(j,i,k)+J(j,i-1,k))/dy^2+(J(j+1,i,k)-2*J(j,i,k)+J(j-1,i,k))/dx^2);
                 alpha = 1/sqrt((1+(I_mag_grad/K)^2));
                  J(j,i,k) = J(j,i,k) - dt/(dx^4)*(alpha)*(J(j-2,i,k)-4*J(j-1,i,k)+6*J(j,i,k)-4*J(j+1,i,k)+J(j+2,i,k)) - ... 
                  dt/(dy^4)*(alpha)*(J(j,i-2,k)-4*J(j,i-1,k)+6*J(j,i,k)-4*J(j,i+1,k)+J(j,i+2,k)); 
           end
       end
       end 
end

% Calculate the max value of the original image
max = 0;
for k = 1:3
    for j = 1:n
        for i = 1:m
            if (max < I2(i,j,k))
                max = I2(i,j,k);
            end
        end
    end
end

% Calculate the PSNR
Temp = double(I2);
sum = 0.0;
for k = 1:3
    for j = 3:n-2
        for i = 3:m-2
            sum = sum + (Temp(i,j,k)-J(i,j,k))^2;
        end
    end
end
MSE = sum/(m*n*o);
PSNR = 20*log10(double(max)/sqrt(MSE));

% Output filtered image
figure(3)
J = uint8(J);
image(J)
% Create a profile if desired
% improfile




