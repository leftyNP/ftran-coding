% This script loads all the data from the tests run to benchmark the convergence for the CRASH review.
% There are two types of data, the spatial (4) and the temporal (1)
% Spatial results have a single, steady-state file.
% Temporal results have 20 output files for each run, and a single file that contains the 
% 3 error types, calculated in Fortran with the run.

close all
clc
tic();
%--------------------------------------------------------------------------------------------------%
loadfiles = 1;

%--------------------------------------------------------------------------------------------------%
% Load files needed
if (loadfiles == 1)

	clear all
	for i=1:5
		cd /Users/leftynm/Desktop/Research/ftran/Mixed-Cells/2011-11-1/Convergence
		if 		i==1
			cd 010
			size_str = '010';
		elseif 	i==2
			cd 020
			size_str = '020';
		elseif 	i==3
			cd 040
			size_str = '040';
		elseif 	i==4
			cd 080
			size_str = '080';
		elseif 	i==5
			cd 160
			size_str = '160';
		end

		cd 1
		for j=1:4
			if 		j==1
				cd ../1
			elseif 	j==2
				cd ../2
			elseif 	j==3
				cd ../3
			elseif 	j==4
				cd ../4
			end		

			c = load('C/c00001.m'); a = load('C/a00000.m');
			
			diffstr=sprintf('%2s%1i%1s%3s%5s','d_',j,'_',size_str,'=a-c;');
			eval(diffstr);
			clear a c
		end
	end
	cd /Users/leftynm/Desktop/Research/ftran/Mixed-Cells/2011-11-1/Convergence
	
	e5_0 = load('Time/0/Error_20x40.m');
	e5_1 = load('Time/1/Error_20x40.m');
	e5_2 = load('Time/2/Error_20x40.m');
	e5_3 = load('Time/3/Error_20x40.m');
	e5_4 = load('Time/4/Error_20x40.m');
end
%--------------------------------------------------------------------------------------------------%
% Calculate error
timeschecked = [3,5,12,17];
dx = 1E-1./2.^(0:4);
errtable = zeros(5,5);

for type=1:4			% simulation type: 	1-4 (represents 4 different problems)
	for i=1:5			% grid size (dx) : 	tsp=grid size n=0:4, dt=1E-1/2^n
		if 		i==1
			size_str = '010';
		elseif 	i==2
			size_str = '020';
		elseif 	i==3
			size_str = '040';
		elseif 	i==4
			size_str = '080';
		elseif 	i==5
			size_str = '160';
		end
		
		% create a string of the variable name
		varname=sprintf('%2s%1i%1s%3s','d_',type,'_',size_str);
		
		% create a string that calculates the error
		errcalc=sprintf('%14s%7s%8s%8.5f','sqrt(sum(sum((',varname,').^2)))*',dx(i));
		
		% populate a table with the result of evaluating the error-string 
		errtable(i,type) = eval(errcalc);
	end
end
%--------------------------------------------------------------------------------------------------%
% Include temporal error as the last column of the error table
errtable(1,5) = mean(e5_0(:,1));
errtable(2,5) = mean(e5_1(:,1));
errtable(3,5) = mean(e5_2(:,1));
errtable(4,5) = mean(e5_3(:,1));
errtable(5,5) = mean(e5_4(:,1));

%--------------------------------------------------------------------------------------------------%
% Calculate the convergence rate
conv = zeros(4,5);
for i=1:4		% time of simulation being evaluated
	for type=1:4	% type of simulation used
		conv(i,type) = log(errtable(i,type)/errtable(i+1,type))/log(dx(i)/dx(i+1));
	end
	type = 5;
	conv(i,type) = log(errtable(i,type)/errtable(i+1,type))/log(4);
end
%--------------------------------------------------------------------------------------------------%

%t=linspace(1e-2,5,100);
%dt=1.E-4./2.^(0:5);
%--------------------------------------------------------------------------------------------------%
% Plot each run's error WRT time
%figure(1,"visible","off")
%avg=zeros(6,3,2);
%for type = 1:2 		% 1=driven, 2=evolved
%	if type==1; ed='d'; end
%	if type==2; ed='e'; end
%	s=0;
%	for tstep=0:5		% time step size used
%		name=sprintf('%1s%1s%1i','e',ed,tstep);
%		plotcmd=sprintf('%11s%1s%3s%1s','semilogy(t,',name,')');
%		filename=sprintf('%10s%3s%4s','ErrorPlot-',name,'.eps');
%		
%		eval(plotcmd)
%		xlabel('time')
%		ylabel('Lp error')
%		legend('L_1','L_2','L_\infty')
%		print('-depsc2',filename)
%		
%		s++;
%		avg(s,1:3,type)=mean(eval(name));
%	end
%end
%clf
%figure(1,"visible","on")
%close
%%--------------------------------------------------------------------------------------------------%
%% calculate convergence factor - base for log is irrelvant due to ratio
%k=zeros(size(avg));
%for i=1:5				% time steps used (10, 20, 40, 80, 160, 320)
%	for j=1:3			% error type used (Emax, E2, Relative E2)
%		for type=1:2	% run type (driven or evolved)
%			k(i,j,type) = log(avg(i,j,type)/avg(i+1,j,type)) / log(dt(i)/dt(i+1));
%		end 
%	end
%end
%--------------------------------------------------------------------------------------------------%
toc()
