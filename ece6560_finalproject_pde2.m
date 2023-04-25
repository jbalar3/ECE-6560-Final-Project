% Read in image
I2 = imread("ece6560_img4.jpg");
figure(1)
image(I2)
% Set dimensions
[m, n, o] = size(I2);
% Add noise if necessary, else uncomment line 9
J = imnoise(I2,'gaussian');
% J = I2;
% Output noisy image
figure(2)
image(J)

J = double(J);
% Set step size for x and y direction and determine the time step
dx = 1/(m-1);
dy = 1/(n-1);
dt = 1e-7;

% Set number of iterations
numt = 20;
% Set constant K
K = 1e4;
for l=1:numt
%    t(k+1) = t(k) + dt;
      for k =1:3
       for j = 2:m-1
           for i=2:n-1
                 grad = [(J(j+1,i,k) - J(j-1,i,k))/(2*dx),(J(j,i+1,k) - J(j,i-1,k))/(2*dy)];
                 alpha = 1/(sqrt(1+(norm(grad)/K)^2));
              J(j,i,k) = J(j,i,k) + dt/dx^2*(alpha)*(J(j+1,i,k)-2*J(j,i,k)+J(j-1,i,k)) + ... 
              dt/dy^2*(alpha)*(J(j,i+1,k)-2*J(j,i,k)+J(j,i-1,k)); 
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

figure(3)
J = uint8(J);
image(J)
% Create improfile if desired
% improfile

