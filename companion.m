function[F] = companion(beta)
n=size(beta,2);
p=size(beta,1)/n;

F=zeros(n*p,n*p);

F(1:n,:)=beta';
F(n+1:end,1:n*p-n)=eye(size(F(n+1:end,1:n*p-n)));
end