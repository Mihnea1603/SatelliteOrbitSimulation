function dd_r = g(r, d_r)
G=6.674e-11;
M=5.972e24;
R=6371000;
J2=1.08262668e-3;
w=7.2921e-5;


x=r(1); y=r(2); z=r(3);
norm_r=norm(r);
dd_r=-G*M*r/norm_r^3 ...
    -3/2*J2*G*M*R^2*[x-5*x*z^2; y-5*y*z^2; 3*z-5*z^3]/norm_r^7 ...
    +w^2*[x;y;0] ...
    +2*w*[d_r(2);-d_r(1);0];
end