a = magic(5);
A = ones(5,5,2);
A(:,:,1) = a;
A(:,:,2) = a';
%%
x = padarray(A,[1,1]);
x(1,:,1) = x(2,:,1);
x(:,1,1) = x(:,2,1);
x(end,:,1) = x(end-1,:,1);
x(:,end,1)=x(:,end-1,1)
x(1,:,2) = x(2,:,2);
x(:,1,2) = x(:,2,2);
x(end,:,2) = x(end-1,:,2);
x(:,end,2)=x(:,end-1,2)