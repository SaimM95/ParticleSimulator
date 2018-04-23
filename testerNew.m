function [outputs] = testerNew(tMax)
	% P = [500,350;400,450];
	% V = [0,0;0,0];
	% M = [1.9*10^30;5.9*10^24];
	P = [500,500;800,500];
	V = [0.003,0.003;0,0];
	M = [15*10^6;50000];

	N = size(P)(1); % num rows in P
	t = 0;
	dt = 1;
	G = 6.673*10^-11;

	while (t < tMax)
		for i=1:N
			P(i,:) += V(i,:) * dt;
		end

		for i=1:N
			ax = 0;
			ay = 0;
			for j=1:N
				if (i != j)
					dx = P(i,1) - P(j,1);
					dy = P(i,2) - P(j,2);
					r = sqrt(dx*dx + dy*dy) / 100;
					a = (-1 * G * M(j)) / (r*r);
					ax += a/r * dx;
					ay += a/r * dy;
				end
			end
			V(i,1) += ax * dt;
			V(i,2) += ay * dt;
		end

		t += dt;
	end
	P
end
