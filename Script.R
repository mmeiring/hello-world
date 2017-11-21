
# Loading the necessary libraries

library(deSolve)
library(ggplot2)
library(Deriv)

# ===================================================================================

# Reading in the data. These data are from the MIDPOINT of the common carotid artery.
DCG_CC <- read.csv(file = "Boileau2015 Common Carotid/Data/DCG.txt", header = FALSE, sep = "	")
rxDataStep(inData = DCG_CC, outFile = "DCG_CC.xdf",
 	transforms = list(deltaRadiusMetres = V4*10^(-3)) ,overwrite = TRUE)

# Renaming the variables as per the accompanying text file.
varInfo <- list(V1 = list(newName = "time", description = "Time in seconds"),
	V2 = list(newName = "pressure", description = "Pressure measurement in kPa"),
	V3 = list(newName = "flowRate", description = "Flow rate in ml/s"),
	V4 = list(newName = "deltaRadius", description = "Change in radius from diastole (mm)"),
	deltaRadiusMetres = list(description = "Change in radius from diastole (m)"),
	V5 = list(newName = "deltaPressure", description = "Pressure gradient between the inlet and outlet in kPa"))

rxSetVarInfo(varInfo = varInfo, data = "DCG_CC.xdf")
rxGetInfo(data = "DCG_CC.xdf", numRows = 5)
rxGetVarInfo(data = "DCG_CC.xdf")

rxSummary(~., data = "DCG_CC.xdf")

# Converting the .xdf file into a dataframe so that we can create a ggplot2 object.
DCG_CC_data <- rxXdfToDataFrame(file = "DCG_CC.xdf")

ggCC <- ggplot(data = data.frame(DCG_CC_data))

ggCC + geom_line(aes(x = time, y = pressure), size = 0.5) + scale_x_continuous(limits = c(8.8, 9.9)) + scale_y_continuous(limits = c(9, 17), breaks = seq(9, 17, by = 1)) + theme_light()
ggCC + geom_line(aes(x = time, y = flowRate), size = 0.5) + scale_x_continuous(limits = c(8.8, 9.9), breaks = seq(8.8, 9.9, by = 0.1)) +theme_light()

ggCC + geom_line(aes(x = time, y = flowRate), size = 0.5) + theme_light()

# ============================================================================================

# Prescribed volumetric inflow in ml/s.
InFlow <- function(t, Ts){
	(6.5+3.294*sin(2*pi*t/Ts-0.023974)+1.9262*sin(4*pi*t/Ts-1.1801)-1.4219*sin(6*pi*t/Ts+0.92701)-0.66627*sin(8*pi*t/Ts-0.24118)-0.33933*sin(10*pi*t/Ts-0.27471)-0.37914*sin(12*pi*t/Ts-1.0557)+0.22396*sin(14*pi*t/Ts+1.22)+0.1507*sin(16*pi*t/Ts+1.0984)+0.18735*sin(18*pi*t/Ts+0.067483)+0.038625*sin(20*pi*t/Ts+0.22262)+0.012643*sin(22*pi*t/Ts-0.10093)-0.0042453*sin(24*pi*t/Ts-1.1044)-0.012781*sin(26*pi*t/Ts-1.3739)+0.014805*sin(28*pi*t/Ts+1.2797)+0.012249*sin(30*pi*t/Ts+0.80827)+0.0076502*sin(32*pi*t/Ts+0.40757)+0.0030692*sin(34*pi*t/Ts+0.195)-0.0012271*sin(36*pi*t/Ts-1.1371)-0.0042581*sin(38*pi*t/Ts-0.92102)-0.0069785*sin(40*pi*t/Ts-1.2364)+0.0085652*sin(42*pi*t/Ts+1.4539)+0.0081881*sin(44*pi*t/Ts+0.89599)+0.0056549*sin(46*pi*t/Ts+0.17623)+0.0026358*sin(48*pi*t/Ts-1.3003)-0.0050868*sin(50*pi*t/Ts-0.011056)-0.0085829*sin(52*pi*t/Ts-0.86463))
}

# Derivative of the prescribed inflow rate. This function is used in the Windkessel model.
DerivInFlow <- function(t, Ts){
	.e2 <- pi * t/Ts
	pi * (0.1043528 * cos(0.195 + 34 * .e2) + 0.1265184 * cos(48 * .e2 - 1.3003) + 0.2448064 * cos(0.40757 + 32 * .e2) + 0.2601254 * cos(0.17623 + 46 * .e2) + 0.278146 * cos(22 * .e2 - 0.10093) + 0.3597384 * cos(1.4539 + 42 * .e2) + 0.3602764 * cos(0.89599 + 44 * .e2) + 0.36747 * cos(0.80827 + 30 * .e2) + 0.41454 * cos(1.2797 + 28 * .e2) + 0.7725 * cos(0.22262 + 20 * .e2) + 2.4112 * cos(1.0984 + 16 * .e2) + 3.13544 * cos(1.22 + 14 * .e2) + 3.3723 * cos(0.067483 + 18 * .e2) + 6.588 * cos(2 * .e2 - 0.023974) + 7.7048 * cos(4 * .e2 - 1.1801) - (0.0441756 * cos(36 * .e2 - 1.1371) + 0.1018872 * cos(24 * .e2 - 1.1044) + 0.1618078 * cos(38 * .e2 - 0.92102) + 0.25434 * cos(50 * .e2 - 0.011056) + 0.27914 * cos(40 * .e2 - 1.2364) + 0.332306 * cos(26 * .e2 - 1.3739) + 0.4463108 * cos(52 * .e2 - 0.86463) + 3.3933 * cos(10 * .e2 - 0.27471) + 4.54968 * cos(12 * .e2 - 1.0557) + 5.33016 * cos(8 * .e2 - 0.24118) + 8.5314 * cos(0.92701 + 6 * .e2)))/Ts
}

# Time steps to be used in the simulations and plotting exercises.
times <- seq(0, 10, length.out = 1e5)

# Simply plotting the volumetric inflow boundary condition.
plot(times, InFlow(times, 1.155), type = "l")

# This is the "typical" 3-element Windkessel model: a resistor and capacitor in parallel connected to a second resistor in series and is used as a terminal boundary condition for all the models in the current paper.

Windkessel_3E <- function(t, x, params) { 
 with(as.list(c(params, x)), { 
 di <- DerivInFlow(t, Ts)
 dP <- R1*di - P/(R2*C) + (R2 + R1)*i/(R2*C)
 res <- c(di, dP) 
 list(res) 
 }) 
} 

# Parameters used in the above Windkessel model; these are obtained from table III on page 9 of the current paper. The original values given in the aforementioned table are
# R1 = 2.4875e8, R2 = 1.8697e9 and C = 1.7529e-10.
# These must, however, be adjusted in order to be converted to ml from m3 to be consistent with the inflow boundary condition.
params <- c(Ts = 1.1, R1 = 2.4875e2, R2 = 1.8697e3, C = 1.7529e-4)

# Start values for steady state 
y <- c(i = InFlow(0, 1.1), P = 0.12465e3)

# Solving 
out.3E <- lsoda(y, times, Windkessel_3E, params)

gg.3E <- ggplot(data = data.frame(out.3E))

gg.3E + geom_line(aes(x = time, y = P/1000), size = 0.5) +
geom_line(aes(x = time, y = pressure), colour = "red", size = 0.5, data = DCG_CC_data) + theme_light()

# =============================================================================

# This is a very simple 2-element Windkessel. This is included here to merely for comparison against the more accurate 3-element.
Windkessel_2E <- function(t, x, params) { 
 with(as.list(c(params, x)), { 
 dP <- (InFlow(t, Ts) - P/R)/C
 res <- c(dP) 
 list(res) 
 }) 
}

params <- c(Ts = 1.1, R = 1.8697e3, C = 1.7529e-4)

# Start values for steady state 
y <- c(P = 0.12465e3)

# Solving 
out.2E <- lsoda(y, times, Windkessel_2E, params)

gg.2E <- ggplot(data = data.frame(out.2E))

gg.2E + geom_line(aes(x = time, y = P/1000), colour = "blue", size = 0.5) +
geom_line(aes(x = time, y = P/1000), colour = "red", size = 0.5, data = data.frame(out.3E)) +
geom_line(aes(x = time, y = pressure), size = 0.5, data = DCG_CC_data) + theme_light()







