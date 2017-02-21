function x = imspecfunc(x, func)
% Apply FUNC to the eigenvalues of the matrix field X

[~, ~, D, D] = size(x);

switch D
    case 1
        x = func(x);
    case 2
        a = squeeze(x(:, :, 1, 1));
        b = squeeze(x(:, :, 1, 2));
        c = squeeze(x(:, :, 2, 1));
        d = squeeze(x(:, :, 2, 2));

        delta = sqrt(4*b.*c+(a-d).^2);
        e1 = func((a+d+delta)/2);
        e2 = func((a+d-delta)/2);

        x(:, :, 1, 1) = 1./(2*delta).*((a-d+delta).*e1-(a-d-delta).*e2);
        x(:, :, 1, 2) = 1./(delta).*(b.*(e1-e2));
        x(:, :, 2, 1) = 1./(delta).*(c.*(e1-e2));
        x(:, :, 2, 2) = 1./(2*delta).*(-(a-d-delta).*e1+(a-d+delta).*e2);
    case 3
        a = squeeze(x(:, :, 1, 1));
        b = squeeze(x(:, :, 1, 2));
        c = squeeze(x(:, :, 1, 3));
        d = squeeze(x(:, :, 2, 1));
        e = squeeze(x(:, :, 2, 2));
        f = squeeze(x(:, :, 2, 3));
        g = squeeze(x(:, :, 3, 1));
        h = squeeze(x(:, :, 3, 2));
        i = squeeze(x(:, :, 3, 3));

        tr = a+e+i;
        x1 = a.^2+3*b.*d+e.^2+3*c.*g+3*f.*h-e.*i+i.^2-a.*(e+i);
        x2 = -2*a.^3 - 2*e.^3 + 18*c.*e.*g - 27*c.*d.*h - 9*e.*f.*h ...
             +3*(e.^2 - 3*c.*g - 3*f.*h).*i + 3*e.*i.^2 - 2*i.^3 + 3*a.^2.*(e + i)...
             -9*b.*(d.*e + 3*f.*g -2*d.*i) + 3*a.*(-3*b.*d + e.^2 -3*c.*g...
                                                   + 6*f.*h -4*e.*i + i.^2);
        x3 = (x2+sqrt(-4*x1.^3+x2.^2)).^(1/3);

        lambda1 = 1/6*(2*tr-2^(4/3)*x1./x3-2^(2/3)*x3);
        lambda2 = 1/12*(4*tr+2^(4/3)*(1+1i*sqrt(3))*x1./x3+2^(2/3)*(1-1i*sqrt(3))*x3);
        lambda3 = 1/12*(4*tr+2^(4/3)*(1-1i*sqrt(3))*x1./x3+2^(2/3)*(1+1i*sqrt(3))*x3);

        e1 = func(lambda1);
        e2 = func(lambda2);
        e3 = func(lambda3);

        v11 = (-i+lambda1+(h.*(f.*g-d.*(i-lambda1)))./(-d.*h+g.*(e-lambda1)))./g;
        v12 = -(f.*g-d.*(i-lambda1))./(-d.*h+g.*(e-lambda1));
        v13 = 1;
        v21 = (-i+lambda2+(h.*(f.*g-d.*(i-lambda2)))./(-d.*h+g.*(e-lambda2)))./g;
        v22 = -(f.*g-d.*(i-lambda2))./(-d.*h+g.*(e-lambda2));
        v23 = 1;
        v31 = (-i+lambda3+(h.*(f.*g-d.*(i-lambda3)))./(-d.*h+g.*(e-lambda3)))./g;
        v32 = -(f.*g-d.*(i-lambda3))./(-d.*h+g.*(e-lambda3));
        v33 = 1;
        n1 = sqrt(real(v11).^2+imag(v11).^2+real(v12).^2+imag(v12).^2+real(v13).^2+imag(v13).^2);
        n2 = sqrt(real(v21).^2+imag(v21).^2+real(v22).^2+imag(v22).^2+real(v23).^2+imag(v23).^2);
        n3 = sqrt(real(v31).^2+imag(v31).^2+real(v32).^2+imag(v32).^2+real(v33).^2+imag(v33).^2);
        v11 = v11./n1;
        v12 = v12./n1;
        v13 = v13./n1;
        v21 = v21./n2;
        v22 = v22./n2;
        v23 = v23./n2;
        v31 = v31./n3;
        v32 = v32./n3;
        v33 = v33./n3;

        x(:, :, 1, 1) = e1.*v11.*conj(v11) + e2.*v21.*conj(v21) + e3.*v31.*conj(v31);
        x(:, :, 1, 2) = e1.*v11.*conj(v12) + e2.*v21.*conj(v22) + e3.*v31.*conj(v32);
        x(:, :, 1, 3) = e1.*v11.*conj(v13) + e2.*v21.*conj(v23) + e3.*v31.*conj(v33);

        x(:, :, 2, 1) = e1.*v12.*conj(v11) + e2.*v22.*conj(v21) + e3.*v32.*conj(v31);
        x(:, :, 2, 2) = e1.*v12.*conj(v12) + e2.*v22.*conj(v22) + e3.*v32.*conj(v32);
        x(:, :, 2, 3) = e1.*v12.*conj(v13) + e2.*v22.*conj(v23) + e3.*v32.*conj(v33);

        x(:, :, 3, 1) = e1.*v13.*conj(v11) + e2.*v23.*conj(v21) + e3.*v33.*conj(v31);
        x(:, :, 3, 2) = e1.*v13.*conj(v12) + e2.*v23.*conj(v22) + e3.*v33.*conj(v32);
        x(:, :, 3, 3) = e1.*v13.*conj(v13) + e2.*v23.*conj(v23) + e3.*v33.*conj(v33);
    otherwise
        error(sprintf('Not implemented for %d x %d matrices', D));
end
