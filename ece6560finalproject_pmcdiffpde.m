% Read in image
figure(1)
I2 = imread("ece6560_img4.jpg");
image(I2)
% Set dimensions
[m, n, o] = size(I2);
% Add noise if needed, else uncomment line 9
 J = imnoise(I2,'gaussian');
% J = I2;
% Output noisy image
figure(2)
image(J)

J = double(J);
% Set step size for x and y
dx = 1/(m-1);
dy = 1/(n-1);

% Set iteration number and time step and value for K
numt = 8;
dt = 2e-7;
K = 4e-6;
% Run the discrete model
for l=1:numt
       for k =1:3
       for i = 2:m-1
           for j=2:n-1
               grad = [(J(i+1,j,k) - J(i-1,j,k))/(2*dx),(J(i,j+1,k) - J(i,j-1,k))/(2*dy)];

               J(i,j,k) = J(i,j,k) + dt*((1/dx^2)*(J(i+1,j,k)-2*J(i,j,k)+J(i-1,j,k)) + (1/dy^2)*(J(i,j+1,k)-2*J(i,j,k)+J(i,j-1,k)) + ... 
                   K^2*((1/(dx^2))*((J(i+1,j,k)-J(i-1,j,k))/2)^2*(1/dy^2)*((J(i,j+1,k)-2*J(i,j,k)+J(i,j-1,k)))-2*(1/(dx))*((J(i+1,j,k)-J(i-1,j,k))/(2))*...
                   (1/(dy))*((J(i,j+1,k)-J(i,j-1,k))/(2))*(1/(dx*dy))*(0.25*(J(i-1,j-1,k)+J(i+1,j+1,k)-J(i+1,j-1,k)-J(i-1,j+1,k)))+(1/(dy^2))*((J(i,j+1,k)-J(i,j-1,k))/(2))^2*...
                   (1/dx^2)*((J(i+1,j,k)-2*J(i,j,k)+J(i-1,j,k)))))/(1+K^2*(norm(grad)))^2;
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

% Output filtered model
figure(3)
J = uint8(J);
image(J)
% improfile



