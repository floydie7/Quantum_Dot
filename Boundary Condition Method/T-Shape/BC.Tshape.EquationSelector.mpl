EquationSelector:=proc(state::string)::list,list:
	description "Selects system of equations describing boundary condition based off user input":
	option `Copyright (c) 2016, Benjamin Floyd`:
	local K, Kp, Eq1, Eq2, Eq3, Eq4, nlb, mlb:

	
	# Set system of equations
	if evalb(state = "ground") or evalb(state = "first") then 				# Even parity states
		# Wave numbers
		K:=j->sqrt((j*Pi)^2*(1/alpha)^2-E):
		Kp:=j->alpha*sqrt(((2*j+1)*Pi)^2-E):
		# System of equations satisfying boundary conditions
		Eq1:=a[n]*exp(-K(n)/2)-c[n]*cosh(K(n)/2)=0:
		Eq2:=a[n]*(-K(n))*exp(-K(n)/2)-c[n]*K(n)*sinh(K(n)/2)-sum(d[m]*(-1)^(m+n)*2*n*(2*m+1)*Pi^2*sinh(Kp(m))/(n^2*Pi^2+Kp(m)^2),m=mlb..ub)=0:
		Eq3:=b[m]*exp(-Kp(m))-d[m]*sinh(Kp(m))=0:
		Eq4:=b[m]*(-Kp(m))*exp(-Kp(m))-d[m]*Kp(m)*cosh(Kp(m))-sum(c[n]*(-1)^(m+n)*4*n*(2*m+1)*Pi^2*cosh(K(n)/2)/(((2*m+1)*Pi)^2+K(n)^2),n=nlb..ub)=0:
	else error "Parameter must be either 'ground' or 'first'."
	fi:

	# Set lower bounds on summations
	if evalb(state = "ground") then
		mlb:=1: nlb:=2:
	elif evalb(state = "first") then
		mlb:=2: nlb:=3:
	else error "Parameter must be either 'ground' or 'first'."
	fi:

	return [Eq1,Eq2,Eq3,Eq4],[mlb,nlb]
end proc:
