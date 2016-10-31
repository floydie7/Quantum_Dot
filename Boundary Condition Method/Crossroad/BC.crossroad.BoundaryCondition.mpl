BoundaryCondition:=proc(state::string,initlsearchbnd::float,initusearchbnd::float,widthparam::{integer,float}:=1,MinMatSize::integer:=0,MaxMatSize::integer:=0,step::float:=1e-17,Tol::float:=1e-10,Nmax::integer:=100)::list:
	description "Computes the eigenvalue corresponding to specified eigenstate of quantum dot":
	option `Copyright (c) 2016, Benjamin Floyd`:
	global eqns, bnds, M, lsearchbnd, usearchbnd, i:
	local ev0, ev, evs:

	evs:=[]:

	
	# Input specified state and generate the correct system of equations and the corresponding indexing lower bounds
	eqns:=EquationSelector(state)[1]:
	bnds:=EquationSelector(state)[2]:

	# Initial eigenvalue search
	
	# Build coefficient matrix based on system of equations
	M:=MatrixGenerator(eqns,bnds,MinMatSize):
	
	# Compute eigenvalue using a bisection search
	ev0:=MatrixBisection(M,initlsearchbnd,initusearchbnd,widthparam,Tol,Nmax):

	# Store value in output list
	evs:=[evs[],ev0]:
	

	# Test to see if multiple searches are requested
	if MaxMatSize > MinMatSize then
		for i from MinMatSize+1 to MaxMatSize do
			M:=MatrixGenerator(eqns,bnds,i):								# Build matrix
			lsearchbnd:=evs[-1]-step: usearchbnd:=evs[-1]+step:				# Update search bounds based on previous result
			ev:=MatrixBisection(M,lsearchbnd,usearchbnd,widthparam,Tol,Nmax):	# Compute new eigenvalue
			evs:=[evs[],ev]:											# Store new value in output list
			gc():													# Preform garbage collection
		od:
	fi:

	return evs
end proc:
