% The m-file loads a sparse matrix into Matlab, converts it to a full matrix, and then calculates
% a measure of the diagonal dominance.

tic();

% Load the sparse matrix and determine where the data ends and zeros begin
B=load('C/A.m');
for n=1:size(B,1)
	if B(n,1)==0
		n--;
		break
	end
end

% Create a full matrix and populate it with the sparse matrix information
A=0;
for s=1:n
	i=floor(B(s,1));
	j=floor(B(s,2));
	A(i,j) = B(s,3);
end

% Show a measure of symmetry. A symmetric matrix will give 0
disp("Symmetry Measure: sum(sum(abs(A-A')")
disp(sum(sum(abs(A-A'))))

n = size(A,1);
% Create a vector of the diagonals
diag = zeros(n,1);
for i=1:n, diag(i)=abs(A(i,i)); end

% Show a measure of diagonal dominance (ddom = diagonal dominance)
ddom = (sum(abs(A),2)-diag);

disp('Diagonal Dominance Measure: % of rows that are diagonal dominant')
disp(sum(diag > ddom)/n)

% Show a few random vector multiplications to determine Positive-Definiteness
PD = zeros(1,10);
for i=1:10
	vec = rand(n,1);
	PD(i)=vec'*A*vec;
end

disp('Positive-Definiteness of matrix - 10 examples')
disp(PD)

clear B i j s PD vec
toc();