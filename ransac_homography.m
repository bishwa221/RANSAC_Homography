function best_H = ransac_homography(p1,p2)
    thresh = sqrt(2); % threshold for inlier points
    p = 1-1e-4; % probability of RANSAC success
    w = 0.5; % fraction inliers
	
    % n: number of correspondences required to build the model (homography)
    n = 4;
    % number of iterations required
    % from the lecture given the probability of RANSAC success, and fraction of inliers
    k = ceil(log(1-p)/log(1-(w^n)));
	num_pts = size(p1,1);
    best_inliers = 4;
    best_H = eye(3);
    [np1,~] = size(p1);
    [np2,~] = size(p2);
    for iter = 1:k
        % randomly select n correspondences from p1 and p2
        % use these points to compute the homography
        temp = randperm(np1);
        p1_sample = p1(temp(1:n),:);
        p2_sample = p2(temp(1:n),:);
        H = compute_homography(p1_sample,p2_sample);
        % transform p2 to homogeneous coordinates
        p2_h = cart2hom(p2);
        % estimate the location of correspondences given the homography
        p1_hat = H*p2_h;
        % convert to image coordinates by dividing x and y by the third coordinate
        temp = p1_hat(3,:);
        p1_hat = p1_hat./temp;
        % compute the distance between the estimated correspondence location and the 
        % putative correspondence location
	dist = pdist2(p1_hat(1:2,:)', p1, 'euclidean');
        % inlying points have a distance less than the threshold thresh defined previously
	num_inliers = sum(sum(dist < thresh));
		
	if num_inliers > best_inliers
	    best_inliers = num_inliers;
	    best_H = H;
        end
    end
end