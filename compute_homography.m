function H = compute_homography(p1,p2)		
    % use SVD to solve for H as was done in the lecture
    % Note: p1 is a 4x2 matrix, the first column of which gives the x coordinates. 
    [r1,c1]=size(p1);
    A = zeros(8,9);
    for i = 1:r1
        ac = p1(i,:);%ac=[xi,yi]
        ac2 = [ac,1];%ac2=[xi,yi,1];
        for j = 2*i-1:2*i
            A(j,7:9) = ac2; 
            if mod(j,2)~=0%if odd row
                A(j,1:3) = ac2;
            else
                A(j,4:6) = ac2;
            end
        end
    end
    for k = 1:r1
        vx = -p2(k,1);
        vy = -p2(k,2);
        for m = 2*k-1:2*k
            if mod(m,2)~=0
                A(m,7:9) = A(m,7:9).*vx;
            else
                A(m,7:9) = A(m,7:9).*vy;
            end
        end
    end
    [u,d,v] = svd(A);
    X = v(:,end)/v(end,end);
    H=reshape(X,3,3)';
    
end
