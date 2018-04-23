function tester(iterations)
	% P = [500,300;500,500];
	% V = [0,0;0,0];
	% M = [1;10^17];
	H = 1;

	P = [1000,500; 
		 800,500; 
		 500,500];
	V = [0.1,0.1; 0.01,0.01; 0.003,0.003];
	M = [50000; 500000; 15*10^6];

	for t=1:iterations
		% fprintf('Time Step = %d\n', t);

		n = size(P)(1); % num rows in P
		forces = zeros(n,n,2);

		for i=1:n
			for j=i+1:n
				% fprintf('%d,%d\n', i, j);
				[f, H] = calcforce(M(i), M(j), P(i,1), P(i,2), P(j,1), P(j,2), H);
				forces(i,j,1) = f(1);
				forces(i,j,2) = f(2);
			end
		end

		for i=1:n
			for j=1:i-1
				forces(i,j,1) = -1 * forces(j,i,1);
				forces(i,j,2) = -1 * forces(j,i,2);
			end
		end

		% forces

		F = zeros(n,2);

		for i=1:n
			for j=1:n
				% F(i,1) = F(i,1) + forces(i,j,1);
				% F(i,2) = F(i,2) + forces(i,j,2);

				% F(j,1) = F(j,1) - forces(i,j,1);
				% F(j,2) = F(j,2) - forces(i,j,2);
				F(i,1) += forces(i,j,1);
				F(i,2) += forces(i,j,2);
			end
		end

		if (t == iterations)
			F
		end

		for i=1:n
			P(i,:) = P(i,:) + H*V(i,:);
			V(i,:) = V(i,:) + H*F(i,:)*(1/M(i));
		end

		% P
		% V
	end

	P
	V
end

function [force, Hout] = calcforce(m1, m2, p1x, p1y, p2x, p2y, Hin)
	pixelsPerMeter = 100;
	radius = 10;
	Hout = Hin;

	s1 = [p1x;p1y];
	s2 = [p2x;p2y];

	G = 6.673*10^-11;
	d = norm(s1-s2) / pixelsPerMeter;

	if (d*pixelsPerMeter < radius)
		d
		force = [0,0];
	else
		force = -G * ((m1*m2)/d^3) * (s1 - s2);
	end
end

function velocity = calcVelocity(velPrev, h, forcesPrev, m)
	velocity = velPrev + h*forcesPrev*(1/m);
end
