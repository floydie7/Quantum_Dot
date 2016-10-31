MatrixGenerator:=proc(equations::list,bounds::list,MatSize::integer)::Matrix:
	description "Generates a matrix from the linear system of equations describing the boundary conditions":
	option `Copyright (c) 2016, Benjamin Floyd`:
	uses LinearAlgebra:
	global sys, var:
	local Eq1, Eq2, Eq3, Eq4, mlbnd, nlbnd, M:


	# Extract equations and bounds from input
	Eq1:=equations[1]: Eq2:=equations[2]: Eq3:=equations[3]: Eq4:=equations[4]:
	mlbnd:=bounds[1]: nlbnd:=bounds[2]:

	# Sequence equations to build full system of equations
	sys:=[seq(eval(Eq1),n=nlbnd..MatSize),seq(eval(eval(Eq2,[mlb=mlbnd,ub=MatSize])),n=nlbnd..MatSize),seq(eval(Eq3),m=mlbnd..MatSize),seq(eval(eval(Eq4,[nlb=nlbnd,ub=MatSize])),m=mlbnd..MatSize)]:
	# Build list of coefficients 
	var:=[seq(a[n],n=nlbnd..MatSize),seq(b[m],m=mlbnd..MatSize),seq(c[n],n=nlbnd..MatSize),seq(d[m],m=mlbnd..MatSize)]:

	# Build coefficient matrix using GenerateMatrix procedure from LinearAlgebra package
	M:=GenerateMatrix(sys,var)[1]:
	
	return M
end proc:
