function d_x = f(u, x)
x1=x(1:3); x2=x(4:6); x3=x(7:9); x4=x(10:12);
d_x=[x3; x4; g(x1,x3)+u; g(x2,x4); norm(x1-x2)/norm(x2)];
end