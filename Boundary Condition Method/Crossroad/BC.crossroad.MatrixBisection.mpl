MatrixBisection:=proc(Mat::Matrix,lsearchbnd::float,usearchbnd::float,widthparam::float:=1,Tol::float:=1e-10,Nmax::integer:=100)::float:
	description "Performs a bisection search for the roots of the determinant of the inputed matrix":
	option `Copyright (c) 2016, Benjamin Floyd`:
	uses LinearAlgebra:
	local Ma, Mb, Mc, dMa, dMb, dMc, a, b, c, N:

	# Initialize iteration index and search bounds
	N:=0:
	a:=lsearchbnd: b:=usearchbnd:


	# Bisection Search
	while (b-a)/2 > Tol or N<=Nmax do
		# Evaluates matrix at search bounds and width parameter
		Ma:=evalf(eval(Mat,[E=a,alpha=widthparam])):
		Mb:=evalf(eval(Mat,[E=b,alpha=widthparam])):

		# Compute determinants of matrices
		dMa:=Determinant(Ma):
		dMb:=Determinant(Mb):

		# Adaptive search bound algorithm. If search bounds do not cover root, attempt to correct
		if dMa*dMb > 0 then										
		 	if (dMb-dMa)/(b-a) > 0 then 							# Determinant has positive slope
				while dMa*dMb > 0 do
					if dMa > 0 and dMb > 0 then					# Both points to the right of root, shift lower bound to cover root
						a:=a-0.1:
						Ma:=evalf(eval(Mat,[E=a,alpha=widthparam])):
						dMa:=Determinant(Ma):
					elif dMa < 0 and dMb < 0 then					# Both points to the left of root, shift upper bound to cover root
						b:=b+0.1:
						Mb:=evalf(eval(Mat,[E=b,alpha=widthparam])):
						dMb:=Determinant(Mb):
					fi:
				od:
			elif (dMb-dMa)/(b-a) < 0 then							# Determinant has negative slope
				while dMa*dMb > 0 do
					if dMa < 0 and dMb < 0 then 					# Both points to the right of root, shift lower bound to cover root
						a:=a-0.1:
						Ma:=evalf(eval(Mat,[E=a,alpha=widthparam])):
						dMa:=Determinant(Ma):
					elif dMa > 0 and dMb > 0 then					# Both points to the left of root, shift upper bound to cover root
						b:=b+0.1:
						Ma:=evalf(eval(Mat,[E=b,alpha=widthparam])):
						dMb:=Determinant(Mb):
					fi:
				od:
			fi:
		fi:

		# Update index and compute midpoint
		N:=N+1:
		c:=(a+b)/2:

		# Compute determinant of matrix at midpoint
		Mc:=evalf(eval(Mat,[E=c,alpha=widthparam])):
		dMc:=Determinant(Mc):
	
		# Test new determinant value against tolerance threshold
		if abs(dMc) < Tol then
			return evalf(c)
		fi:

		# Test if root is covered from lower bound to midpoint and update bounds
		if dMa*dMc < 0 then
			b:=c:
		else a:=c:
		fi:
	od:

	return evalf((a+b)/2)
end proc:
