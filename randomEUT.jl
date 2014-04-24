t0=time()
c = 299792458.0

function Efarfield(R,theta,phi,I,f)
	X=R*cos(phi)*sin(theta)
    	Y=R*sin(phi)*sin(theta)
    	Z=R*cos(theta)
	DX = X-I[1]
	DY = Y-I[2]
	DZ = Z-I[3]
	dist = sqrt(DX^2+DY^2+DZ^2)
	phase=2*pi*dist*f/c+I[7]
	ca    = cos(I[4])
	sa    = sin(I[4])
	cb    = cos(I[5])
	sb    = sin(I[5])
	distx = ((-sb)^2+(1-(-sb)^2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ
	disty = (-sb*cb*(1-ca))*DX+((cb)^2+(1-cb^2)*ca)*DY+(sb*sa)*DZ
	distz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
	distxy = sqrt(distx^2+disty^2)
	costheta = distz/dist
	sintheta = distxy/dist
	cosphi   = distx/distxy
	sinphi   = disty/distxy
	L =I[6]*1/dist*(f/c)^2*377
	Exx = exp(1im*phase)*L*((((-sb)^2+(1-(-sb)^2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta))))
	Eyy = exp(1im*phase)*L*(((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)^2+(1-(cb)^2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta))))
	Ezz = exp(1im*phase)*L*(((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta))))
	ETheta= Exx*cos(theta)*cos(phi)+Eyy*cos(theta)*sin(phi)-Ezz*sin(theta)
	EPhi= -Exx*sin(phi)+Eyy*cos(phi)
	#Er = Exx*sin(theta)*cos(phi)+Eyy*sin(theta)*sin(phi)+Ezz*cos(theta)
    return ETheta,EPhi
end


a=1. #EUT radius in m
nf=20;
freq= linspace(50e6,1e9,nf);#Frequency in Hz
ka=2*pi*freq./c*a; #electric size, ka=2*pi/lambda*a

#Radiation pattern
np=360;  #number of points along phi;
nt=180;  #number of points along theta;

phi=linspace(0,2*pi,np);
theta=linspace(0,pi,nt);
R = 100; #Measurement sphere radius in m


P=Array(Float64,np,nt,nf); # Power matrix

n=10; #number of radiating dipoles on the EUT
#generate the dipoles
theta_eut=acos(2*rand(n,1).-1); #uniformly random along theta
phi_eut=2*pi*rand(n,1); #uniformly random along phi
x=a*cos(phi_eut).*sin(theta_eut);
y=a*sin(phi_eut).*sin(theta_eut);
z=a*cos(theta_eut);
tilt=acos(2*rand(n,1).-1);
azimut=2*pi*rand(n,1);
ld=.1; #taille des dipôles en m
amplitude=rand(n,1)*ld; #random amplitude of the currents
phas=2*pi*rand(n,1); #random phase


for t=1:nt
	println(t)
	for p=1:np
		for f=1:nf
			Eth=0
			Eph=0
			for i=1:n
				Et,Ep=Efarfield(R,theta[t],phi[p],[x[i],y[i],z[i],tilt[i],azimut[i],amplitude[i],phas[i]],freq[f])
				Eth+=Et
				Eph+=Ep
			end
			P[p,t,f]=abs(Eth)^2+abs(Eph)^2
		end
	end

end


t1=time()
print(t1-t0)

using PyPlot
pygui(false)
#Radiation diagram   
figure(figsize=(10, 5))
for u=1:nf
    	Pnorm = P[:,:,u]./maximum(P[:,:,u]);
    	pcolor(phi,theta,Pnorm')
	xlim(0,2*pi)
	ylim(0,pi)
	grid()
	xlabel("\$\\phi\$")
	ylabel("\$\\theta\$")
	colorbar()
	clim(0,1)
	f=round(freq[u]/1e6)
	title("\$f=$f\$ MHz")
    	savefig("$u.png",bbox="tight")
    	clf()
end
