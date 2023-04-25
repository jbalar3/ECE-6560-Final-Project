% Read in image
I2 = imread("ece6560_img4.jpg");
figure(1)
image(I2)
% Set dimensions
[m, n, o] = size(I2);
% Add noise or uncomment line 9
J = imnoise(I2,'gaussian');
% J = I2;
figure(2)
image(J)
J = double(J);
% Set step size for x, y, and t
dx = 1/(m-1);
dy = 1/(n-1);
dt = 1e-3;
% Set number of repetitions and K
numt = 5;
K = 1e5;
for l=1:numt
      for k =1:3
       for j = 2:m-1
           for i=2:n-1
                 grad_north = (J(j-1,i,k) - J(j,i,k))/(dx);
                 grad_south = (J(j+1,i,k) - J(j,i,k))/(dx);
                 grad_east = (J(j,i+1,k) - J(j,i,k))/(dy);
                 grad_west = (J(j,i-1,k) - J(j,i,k))/(dy);
                 alpha = 4;
                 c_n = exp(-(grad_north/K)^2);
                 c_s = exp(-(grad_south/K)^2);
                 c_e = exp(-(grad_east/K)^2);
                 c_w = exp(-(grad_west/K)^2);
              J(j,i,k) = J(j,i,k) + dt/alpha*(c_n*grad_north+c_s*grad_south+c_e*grad_east+c_w*grad_west); 
           end
       end
      end 
end

% Calculate PSNR
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

% Output filter image
figure(3)
J = uint8(J);
image(J)
% improfile


