EquationSelector:=proc(state::string)::list,list:
	description "Selects system of equations describing boundary condition based off user input":
	option `Copyright (c) 2016, Benjamin Floyd`:
	local K, Kp, Eq1, Eq2, Eq3, Eq4, nlb, mlb:

	
	# Set system of equations
	if evalb(state = "ground") or evalb(state = "second")					# Even parity states
		then # Wave numbers
			K:=j->sqrt(((2*j+1)/2*Pi*1/alpha)^2-E):
			Kp:=j->alpha*sqrt(((2*j+1)/2*Pi)^2-E):
			# System of equations satisfying boundary conditions
			Eq1:=a[n]*exp(-K(n))-c[n]*cosh(K(n))=0:
			Eq2:=a[n]*(-K(n))*exp(-K(n))-c[n]*K(n)*sinh(K(n))-sum(d[m]*(-1)^(m+n+1)*2*(2*m+1)*(2*n+1)*Pi^2*cosh(Kp(m))/(4*Kp(m)^2+(2*n+1)^2*Pi^2),m=mlb..ub)=0:
			Eq3:=b[m]*exp(-Kp(m))-d[m]*cosh(Kp(m))=0:
			Eq4:=b[m]*(-Kp(m))*exp(-Kp(m))-d[m]*Kp(m)*sinh(Kp(m))-sum(c[n]*(-1)^(m+n+1)*2*(2*m+1)*(2*n+1)*Pi^2*cosh(K(n))/(4*K(n)^2+(2*m+1)^2*Pi^2),n=nlb..ub)=0:
	elif evalb(state = "first") or evalb(state = "third")					# Odd parity states
		then # Wave numbers
			K:=j->sqrt((j*Pi*1/alpha)^2-E):
			Kp:=j->alpha*sqrt((m*Pi)^2-E):
			# System of equations satisfying boundary conditions
			Eq1:=a[n]*exp(-K(n))-c[n]*sinh(K(n))=0:
			Eq2:=a[n]*(-K(n))*exp(-K(n))-c[n]*K(n)*cosh(K(n))-sum(d[m]*(-1)^(m+n+1)*2*m*n*Pi^2*sinh(Kp(m))/(Kp(m)^2+n^2*Pi^2),m=mlb..ub)=0:
			Eq3:=b[m]*exp(-Kp(m))-d[m]*sinh(Kp(m))=0:
			Eq4:=b[m]*(-Kp(m))*exp(-Kp(m))-d[m]*Kp(m)*cosh(Kp(m))-sum(c[n]*(-1)^(m+n+1)*2*m*n*Pi^2*sinh(K(n))/(K(n)^2+m^2*Pi^2),n=nlb..ub)=0:
	else error "Parameter must be either 'ground', 'first', 'second', or 'third'."
	fi:

	# Set lower bounds on summations
	if evalb(state = "ground")
		then nlb:=0: mlb:=0:
	elif evalb(state = "first")
		then nlb:=1: mlb:=1:
	elif evalb(state = "second")
		then nlb:=1: mlb:=1:
	elif evalb(state = "third")
		then nlb:=2: mlb:=2:
	else error "Parameter must be either 'ground', 'first', 'second', or 'third'."
	fi:

	return [Eq1,Eq2,Eq3,Eq4],[mlb,nlb]
end proc:
