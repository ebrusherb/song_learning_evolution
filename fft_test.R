#set the sampling frequency 
 samplingFrequency = 1000;

#create the indexes to sample at
 timeInterval = 1/samplingFrequency;
 signalIndex = seq(0, 1, by=timeInterval);

#define num samples as the dimention of signalIndex. here N = 1000
N = length(signalIndex)

#amplitude of signal 1 
 a1 = 1;  

#amplitude of signal 2 
 a2 = 3;  

#frequency of signal 1 
 f1 = 10; 

#frequency of signal 2
 f2 = 20; 

#10 Hz Sine wave with peak amplitude 2
 signal1 = a1 * sin(2 * pi * f1 * signalIndex);

#20 Hz Sine wave with peak amplitude 3
 signal2 = a2 * sin(2 * pi * f2 * signalIndex);

#input signal is the sum of two signals in this case with frequencies 10 and 20 Hz
 inputSignal = signal1 + signal2 ;

#the input signal
# plot(signalIndex, inputSignal, type = 'l');

#perform the FFT. In this case the number of points (N) will be equal to 1000. Output will be the individual components of the FFT.
 fourierComponents = fft(inputSignal);

#get the absolute value of the coefficients  
 fourierCoefficients = abs(fourierComponents);
 

#Normalize coefficients fig 5 here N = 1000 samples so N/2 = 500
 normalizedFourierComponents = fourierCoefficients / (N/2);
 normalizedFourierComponents[1] = fourierCoefficients[1] / N

#get the first 50 coefficients fig 6
 mainCoeffs = normalizedFourierComponents[1:50];

plot(mainCoeffs[1:50])


##----

N = 100
delta = 1/(N)
x = seq(0,1,by=delta)

a1 = 2
f1 = 10
s1 = a1*sin(2*pi*f1*x)

a2 = 3
f2 = N - f1
s2 = a2*sin(2*pi*f2*x)

##----
s=2;
m=1;
layout(matrix(1:6,nrow=2,byrow=FALSE));
for(f in c(2,4,5)){

subset = 600:1200;
eq_pop = Pm_onepop[[s,f,m]][,Tsteps]
fourComps = fft(eq_pop[subset]);

#get the absolute value of the coefficients  
fourCoeffs = abs(fourComps);
 
#Normalize coefficients fig 5 here N = 1000 samples so N/2 = 500
coeffs = fourCoeffs / (Nm/2);
coeffs[1] = fourierCoefficients[1] / Nm
plot((coeffs[1:35]),type='o',ylim=c(0,0.1));plot(1:length(subset),eq_pop[subset],type='l')
}


##----
sd_vals = c(seq(0.01,1,by=0.005))
Nsd = length(sd_vals)
subset = 1:50
x = matrix(NA,Nsd,1)
logcoeffs = matrix(NA,Nsd,20)
mu = 0 
for(i in 1:Nsd){
	sd = sd_vals[i]
	fourComps = fft(dnorm(seq(-7,7,by=0.01),mu,sd));
	fourCoeffs= abs(fourComps);
	coeffs=fourCoeffs[1:floor(length(fourCoeffs)/2)];
	coeffs=coeffs/length(coeffs);
	l = lm(log(coeffs[2:21]) ~ poly(1:20,2,raw=TRUE));
	x[i] = q
	logcoeffs[i,] = log(coeffs[2:21])
}
l = lm(x[subset]~poly(sd_vals[subset],2,raw=TRUE));
a=as.vector(summary(l)$coefficients[,1]);
round(a,3)
plot(sd_vals[subset]^2,x[subset])

i=i+5;plot(1:20,logcoeffs[i,]);l = lm(logcoeffs[i,] ~poly(1:20,2,raw=TRUE));a=as.vector(summary(l)$coefficients[,1]);lines(a[1]+a[2]*(1:20)+a[3]*(1:20)^2)
	