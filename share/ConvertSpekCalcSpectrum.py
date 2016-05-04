f=open("XRayCabinet_Spectrum.dat")
f2=open("XRayCabinet_Spectrum_MeV.mac","w")
lines = f.readlines()

for line in lines : 
	data = line.split()
	f2.write("/gps/hist/point %f %f\n"%(float(data[0])*0.001,float(data[1]))) 
