% Read in image
I2 = imread("ece6560_img4.jpg");
figure(1)
image(I2)
% Set dimensions
[m, n, o] = size(I2);
% Add noise if necessary, else uncomment line 9
J = imnoise(I2,'gaussian');
% J = I2;
figure(2)
image(J)
J = double(J);
% Set step sizes for x, y,and t
dx = 1/(m-1);
dy = 1/(n-1);
dt = 1e-7;
% Set number of iterations and value for lambda
numt = 20;
lambda = 0.5;
for l=1:numt
      for k =1:3
       for j = 2:m-1
           for i=2:n-1
              J(j,i,k) = J(j,i,k) + dt*lambda*(J(j,i,k)-I2(j,i,k)) + ... 
              dt/dx^2*(1-lambda)*(J(j+1,i,k)-2*J(j,i,k)+J(j-1,i,k))+dt/dy^2*(1-lambda)*(J(j,i+1,k)-2*J(j,i,k)+J(j,i-1,k)); 
           end
       end
      end 
end
% Calculate the PSNR
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
% Output the filtered image
figure(3)
J = uint8(J);
image(J)
% improfile


