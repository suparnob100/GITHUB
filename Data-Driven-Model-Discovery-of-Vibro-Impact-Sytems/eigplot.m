function z=eigplot(A)
[~,lambda] = eig(A);
plot(lambda,'x','LineWidth',4,'MarkerSize',8)

hold on
theta = 0:0.01:2*pi;
plot(exp(1j*theta),'.-','LineWidth',2)
axis square
hold off

if sum(abs(diag(lambda))>1)
    disp('Unstable system');
    max(abs(diag(lambda)))
    z=1;
else
    z=0;
end
end