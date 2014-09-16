# Import the Mantid algorithms  we will be using
from mantid.simpleapi import LoadNexus, LoadSassena, Transpose, Rebin, Scale, SassenaFFT, ConvolveWorkspaces, Fit, mtd, DSFinterp

# Working directory. All files should be placed here
rootd=''  #Update this variable with the path to the workding directory
if not rootd:
	print 'Please update variable rootd with path of working directory'
	
# Input resolution file (the elastic line)
LoadNexus(Filename='{0}/elastic.nxs'.format(rootd), OutputWorkspace='elastic')
# Load experimental quasi-elastic spectra S(Q,E) at T=200K. It contains 9 spectra, each at different Q's (0.3, 0.5, 0.7,..,1.9)
# and defined in the energy domain [-0.15, 0.15]meV
LoadNexus(Filename='{0}/exp200K.nxs'.format(rootd), OutputWorkspace='exp200K')

# Load simulated spectra. We have 14 spectra at different values of the dihedral barrier K
parametervalues='0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15'
for K in parametervalues.split():
	# Load spectra into Mantid workspace. Each file contains the intermediate structure factor I(Q,t) spectra for
	# 9 different values of Q (0.3, 0.5, 0.7,..,1.9). The names of the resulting workspaces are incSM0.03, .., incSM0.15
	LoadSassena(Filename='{0}/K{1}/T_200/fqt_inc_run1_rms2first.h5'.format(rootd,K), TimeUnit=1.0, OutputWorkspace='incSMK{0}'.format(K))
	# What happens if one of the I(Q,t) files contains a different number of Q values? In the following special case, K=0.11, the I*Q,t) file contains spectra
	# for values of Q (0.02, 0.03, 0.04, 0.05,...,2.0). Thus, we have to rebin the spectra in the Q-coordinate
	if K=='0.11':
		Transpose(InputWorkspace='incSMK0.11_fqt.Re', OutputWorkspace='incSMK0.11_fqt.Re')   #from I(Q,t) to I(t,Q)
		Rebin(InputWorkspace='incSMK0.11_fqt.Re', Params=[0.2,0.2,2.0], OutputWorkspace='incSMK0.11_fqt.Re') #Rebin in Q to (0.3, 0.5,..,1.9)
		Transpose(InputWorkspace='incSMK0.11_fqt.Re', OutputWorkspace='incSMK0.11_fqt.Re')  # from I(t,Q) back to I(Q,t)
		# After rebin, the intensities I(Q,t=0) in the resulting spectra are different than for spectra at other K values. We rescale the intensities to match
		Scale(InputWorkspace='incSMK0.11_fqt.Re', factor=0.5, Operation='Multiply', OutputWorkspace='incSMK0.11_fqt.Re')
	# Fourier transform I(Q,t) to S(Q,E). Resulting workspaces are incSMK0.03_sqw .., incSMK0.15_sqw
	SassenaFFT(InputWorkspace='incSMK{0}'.format(K), FFTonlyRealpart=1, DetailedBalance=1, Temp=200)
	# Rebin S(Q,E) in E to the region [-0.2, 0.2] meV, of the same order as the experimental [-0.15, 0.15] but a bit broader
	Rebin(InputWorkspace='incSMK{0}_sqw'.format(K), Params=[-0.2,0.0004,0.2], OutputWorkspace='incSMK{0}_sqw'.format(K))
	# Important remark: the system simulated in one POSS molecule. thus, our simulated system does not take into account the broadening
	# due to motions of the molecule center of mass that take place in the crystalline phase. As a result, the elastic line in the simulation
	# is over-represented. To correct this shortcoming in the simulations we will remove the simulated elastic line from the simulated S(Q,E) 
	# and later include a term represented the elastic line in the fitting model.
	ws=mtd['incSMK{0}_sqw'.format(K)]  # simulated S(Q,E)
	# Remove the elastic line for the spectra at each Q. 
	for iQ in range(9):
		Q=0.03+0.02*iQ  #not neccessary but just to remind the Q-value
		#S(Q,E) for the previous Q. SQE is a list of 1000 intensities versus energy. The elastic line is in the middle of this list, at indexes 499 and 500
		#We remove these intensities just by setting these values to the adjacent ones, which are quasi-elastic in nature.
		SQE=ws.dataY(ih)
		SQE[499]=SQE[498]   
		SQE[500]=SQE[501]
#We match simulated and experimental resolution by convolving simulated S(Q,E) with the experimental resolution.
# The names of the convolved workspaces are simSMK0.03, .., simSMK0.15
for K in parametervalues.split():
	ConvolveWorkspaces(Workspace1='elastic', Workspace2='incSMK{0}_sqw'.format(K), OutputWorkspace='simSMK{0}'.format(K))

#The intensities of these convolved workspaces are order of magnitude higher than the intensities from the experimental S(Q,E). This
# by itself is not a problem, except that the fitting routine will not be able to compe with such a difference. Thus, we rescale the intensities
# of the convolved workspaces in order to bring both simulated and experimental intensities closer to each other.
for K in parametervalues.split():
	Scale(InputWorkspace='simSMK{0}'.format(K), Factor=1.0e-05,Operation='Multiply', OutputWorkspace='simSMK{0}'.format(K))

# Previous to the proper fitting, we will do a Chi2 plot versus dihedral barrier K. This will give us a feeling for how the goodness of the fit
# changes with the value of the dihedral value K. We use algorithm DSFinterp for this purpose.
# We will use the structure factors  simSMK0.03, .., simSMK0.15 from the simulations to construct S(Q,E,K) with
# K running from 0.03 to 0.13 every 0.001
parametervalues='0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14'.split()
inputworkspaces=['simSMK{0}'.format(K) for K in parametervalues]  # S(Q,E) from the simulations
targetvalues=[0.03+0.001*i for i in range(101)] #target K varies from 0.03 to 0.13 every 0.001
outputworkspaces=['out{0}'.format(K) for K in targetvalues] # interpolated S(Q,E)
#From the simulated structure factors, create interpolated structure factors for each target K
DSFinterp(ParameterValues=[float(K) for K in parametervalues], Workspaces=inputworkspaces, LoadErrors=0,\
	TargetParameters=targetvalues, OutputWorkspaces=outputworkspaces,\
	LocalRegression=1, RegressionWindow=4, RegressionType='quadratic')
# For each target K, fit the associated structure factor S(Q=1.9, E, K)to the experimenal S(Q=1.9, E), and extract the Chi2 value from the fitting.
# Recall that in S(Q,E) Q goes from 0.03 to 1.9. We select the structure factor with highest Q, Q=1.9, corresponding to L~3Angstroms, to do the fitting.
# Among the reported Q-values, Q=1.9 corresponds to a lenght scale most sensitive to methyl rotations
chi2values=[]
for K in targetvalues:
	# Our fitting model is a elastic line at Q=1.9 plus our simulated S(Q=1.9,E) plus a linear background
	fit_string    ='name=TabulatedFunction,Workspace=elastic,WorkspaceIndex=8,Scaling=1.0;'
	fit_string +='name=TabulatedFunction,Workspace=out{0},WorkspaceIndex=8,Scaling=1.0;'.format(K)
	fit_string +='name=LinearBackground,A0=0.0,A1=0.0'
	# Algorithm Fit carries out the fitting of the previous model to the experimental S(Q=1.9,E), which is index 8 of workspace exp200K.
	# Notice the energy range is [-0.13, 0.1] meV. Outside this energy range, the signal to noise ration of the experimental S(Q=1.9, E) 
	# is not very reliable.
	Fit(fit_string, InputWorkspace='exp200K', WorkspaceIndex=8, StartX=-0.13, EndX=0.1, CreateOutput = 1 )
	ws=mtd['exp200K_Parameters']  # Workspace resulting from the fitting containing the optimal value of the fitting parameters and the Chi2
	chi2=ws.row(4)['Value']
	print 'K=', K, ' Chi2=', chi2
	chi2values.append(chi2)
for outws in outputworkspaces:  DeleteWorkspace(outws)  #a bit of clean-up.
#Save chi2 values to file chi2vsK.dat for inspection. Plot this file with your favorite plotting tool
buf='#K Chi2\n'
for i in range(len(targetvalues)):  buf += '{0} {1}\n'.format(targetvalues[i], chi2values[i])
open('/tmp/chi2vsK.dat', 'w').write(buf)

# File chi2vsK.dat suggest that the dihedral value yielding a simulated S(Q=1.9,E) that compares optimally with the experimental
# S(Q=1.9, E) is Koptimal ~ 0.065. We will use this value as initial guess for fitting S(Q,E) at all Q values.
# Below are the initial guess of the fitting parameters. 'Scaling' is the intensity of the elastic line, 'Intensity' is the intensity of the
# simulated S(Q,E), and 'K' is the dihedral barrier.
guess={'Scaling':1.0, 'Intensity':1.0, 'K':0.065}
# We will use our simulated structure factors simSMK0.03, .., simSMK0.15 to search for the optimal dihedral barrier
parametervalues='0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14'
inputworkspaces=' '.join( ['simSMK{0}'.format(K) for K in parametervalues.split()])
# Cycle over S(Q,E) for different values of Q. We begin with the highest Q, Q=1.9, and work out way to the smallest Q, Q=0.3, doing
# a fit at S(Q,E) at each Q
for iw in [8,7,6,5,4,3,2,1,0]: 
	Q=0.3+iw*0.2
	# Our model is a elastic line plus a  structure factor interpolated from the simulated S(Q,E), plus a linear background
	fit_string = 'name=TabulatedFunction,Workspace=elastic,WorkspaceIndex={0},Scaling={1},constraints=(0.0001<Scaling);'.format(iw,guess['Scaling'])+\
		'name=DSFinterp1DFit,InputWorkspaces="{0}",ParameterValues="{1}",'.format(inputworkspaces,parametervalues) +\
		'LoadErrors=0,LocalRegression=1,RegressionType=quadratic,RegressionWindow=4,' +\
		'WorkspaceIndex={0},Intensity={1},TargetParameter={2},'.format(iw,guess['Intensity'],guess['K']) +\
		'constraints=(0.0001<Intensity);' +\
		'name=LinearBackground,A0=0.0,A1=0.0'
	Fit(fit_string, InputWorkspace='exp200K', WorkspaceIndex=iw, StartX=-0.13, EndX=0.1, CreateOutput = 1 )
	ws=mtd['exp200K_Parameters'] # Workspace resulting from the fitting containing the optimal value of the fitting parameters and the Chi2
	Koptimal=ws.row(2)['Value']
	print 'Q=', Q, 'L={0:5.1f}'.format( 6.2832/Q), 'K={0:6.4f}'.format(Koptimal), 'Chi2={0:3.1f}'.format(ws.row(5)['Value'])
	# Use the optimal intensities of the elastic line and simulated quasi-elastic line as guesses for the fitting of the S(Q,E) for the next Q
	guess['Scaling']=ws.row(0)['Value']
	guess['Intensity']=ws.row(1)['Value']

#The output of the previous print statement:
# Q= 1.9 L=  3.3 K=0.0649 Chi2=4.1
# Q= 1.7 L=  3.7 K=0.0652 Chi2=3.6
# Q= 1.5 L=  4.2 K=0.0651 Chi2=2.9
# Q= 1.3 L=  4.8 K=0.0650 Chi2=2.5
# Q= 1.1 L=  5.7 K=0.0647 Chi2=2.4
# Q= 0.9 L=  7.0 K=0.0640 Chi2=2.7
# Q= 0.7 L=  9.0 K=0.0631 Chi2=2.2
# Q= 0.5 L= 12.6 K=0.0635 Chi2=2.4
# Q= 0.3 L= 20.9 K=0.0635 Chi2=1.8
# Simulations thus predict optimal K in the range [0.0635, 0.0649], a narrow range if we take into account the wide range of lenght scales.


