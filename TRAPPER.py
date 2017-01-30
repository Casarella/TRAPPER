#RedTrProb.py
#TRAPPER - TRAnsition Probability Processing/computER
#c. Jan 27, 2017 - Clark Casarella

# Updated to output to a LaTeX friendly table to the output file
# Does not take uncertainty in mixing into account

import math as m
import scipy.constants as sc
import numpy as np


csvpath=input("Enter path to csv file (including file extension): ")
print('Using the input parameters from',str(csvpath)+':')
outpath=input("Enter output path/filename (will be a text file): ")
print('Output placed at:',outpath)
output_file=open(outpath,'w')

#see_weisskopf_units=input("Do you want to see the Weisskopf unit conversion? [Y/N]")
#csvpath='162Dy_GRID.csv'
#csvpath='../162Dy_GRID/162Dy_GRID.csv'
#output_file=open('out.TEST','w')

dtype_full=[('E_g','f8'),('E_g_error','f8'),
            ('I_g','f8'),('I_g_error','f8'),('I_g_total','f8'),
            ('delta_mixing','f8'),('delta_upper','f8'),('delta_lower','f8'),
            ('tau','f8'),('tau_up','f8'),('tau_down','f8'),
            ('alpha_conv','f8'),('alpha_conv_error','f8'),
            ('A','int'),('multipolarity','S6'),
            ('E_level','f8')]
ndtype=str
npnames=['E_g','E_g_error','I_g','I_g_error','I_g_total',
	 'delta_mixing','delta_upper','delta_lower',
	 'tau','tau_up','tau_down','alpha_conv','alpha_conv_error',
	 'A','multipolarity','E_level']

csvfile = np.genfromtxt(csvpath,delimiter=",",skip_header=1,names=npnames,dtype=dtype_full)
#print('array:',csvfile)

#Test single input section
#E_g=0.888157
#I_g=174.8
#I_tot=369.3
#delta=0
#tau=2830*10**-15
#alpha_conv=0.0032
#multipolarity='E2'
#A=162

#E_g=1.31303
#I_g=0.428
#I_tot=1
#delta=0.28
#delta_up=0.34
#tau=320*10**-15
#alpha_conv=0
#multipolarity='E2(M1)'
#A=160

def set_multipolarity():
  """
  Returns and extracts the multipolarity of a transition 
  Decoded from UTF-8 encoding on the string keyed 'multipolarity'
  """
  #multipolarity=csvfile[0][14].decode("utf-8")[-1]
  if multipolarity[-1]=='1':
    return 1
  elif multipolarity[-1]=='2':
    return 2
  else:
      #return 2
      return 'E2(M1)'

def BwE(A):
  """
  Weisskopf estimate for an electric type, multipolarity 
  l transition, in units of e^2*fm^l 
  (1 W.u. = XXX e^2fm^l)
  """
  l=set_multipolarity()
  if l=='E2(M1)':
      l=2
  return 0.12**(2*l)/(4*m.pi)*(3/(l+3))**2*A**(2*l/3)

def BwM(A):
  """
  Weisskopf estimate for an magnetic type, multipolarity l transition,
  in units of mu_N^2
  """
  l=set_multipolarity()
  if l=='E2(M1)':
      l=1
  return 0.12**(2*(l-1))*10/m.pi*(3/(l+3))**2*A**(2*(l-1))

def doublefactorial(n):
  """
  Double factorial (every other n factorialed)
  """
  if n <=0:
    return 1
  else:
    return n*doublefactorial(n-2)
  
  
def mult_coefficient():
  """
  This coefficient removes angular momentum mixing from the transition
  probabilities.
  """
  l=set_multipolarity()
  if l=='E2(M1)':
      l=2
  return l*(doublefactorial(2*l+1))**2/(l+1)
  #return l*(2*l+1)**2/(l+1)

#print('Coefficient for L:',mult_coefficient())
#print('BwE:',BwE(162))

def mixing_fraction(delta):
  """
  Multipole mixing fraction for any mixed-multipolarity transitions
  Unitless, and calculates relative E2 strength to M1 B(E2) strength
  """
  #delta=csvfile[1][14][5]
  l=set_multipolarity()
  if l=='E2':
      l=2
  if delta==0 or l==1:
    return 1
  elif delta!=0 and l=='E2(M1)':
    return delta**2/(1+delta**2)


#print(mixing_fraction(0.64))
def BR():
    """
    Returns branching ratio (ratio of measured intensity to total intensity leaving the state)
    """
    return I_g/I_tot

#print('Mixing Fraction Delta:',mixing_fraction(delta))

#units from scipy - generally helps with precision
m_p=sc.value('proton mass energy equivalent in MeV')
hc=sc.value('Planck constant over 2 pi times c in MeV fm')
hbar=sc.value('Planck constant over 2 pi in eV s')/10**6
barn=10**2

def B_coefficients():
    """
    Calculates coefficients for the final B(pl) calculation.
    Makes an exception for E1 transitions, traditionally reported in mW.u.
    """
    l=set_multipolarity()
    if l=='E2(M1)':
      l=2
    if multipolarity=='E1':
        return hbar/(8*m.pi)*mult_coefficient()*hc**(1+2*l)*1000
    else:
        return hbar/(8*m.pi)*mult_coefficient()*hc**(1+2*l)

def units():
    """
    Corrects the units from e^2b^l to W.u.
    """
    l=set_multipolarity()
    if l=='E2(M1)':
      l=2
    if multipolarity[0]=='E':
        return barn**l*sc.alpha*sc.hbar/10**-9*sc.c/sc.e*BwE(A) # check here
    elif multipolarity[0]=='M':
        return hc*sc.alpha*BwM(A)*(hc/(2*m_p))**2 #check here again
    
#print('Units from MeVfm to W.u.:',units())
def latex_friendly_units():
    """
    Returns LaTeX-friendly units for copying-pasting into LaTeX documents
    """
    l=multipolarity
    if l=='E1':
        return 'mW.u.'
    elif l=='E2':
        return 'W.u.'
    elif l=='M1':
        return '$\mu_N^2$'
    else:
        return 'W.u. (mixed)'

def B(tau):
    """
    Calculation of transition probability B(pl) from all inputs necessary
    """
    l=set_multipolarity()
    if l=='E2(M1)':
      l=2
    if l==1:
        return round(mixing_fraction(delta)*BR()/(tau*10**-15*(1+alpha_conv)*E_g**(2*l+1))*B_coefficients()/units(),3)
    else:
        return round(mixing_fraction(delta)*BR()/(tau*10**-15*(1+alpha_conv)*E_g**(2*l+1))*B_coefficients()/units(),2)


#determine delta_upper bounds on error

def mixing_upper_bounds():
    """
    Determines which bound should be used for a particular mixing fraction
    - Used in error propagation -
    If delta <0, then the most E2 mixing will occur at the most negative number
    (delta-delta_lower)
    if delta >0, then the most E2 mixing will occur at the most positive number
    (delta+delta_upper)
    """
    if delta<0:
        return delta-delta_lower
    elif delta>0:
        return delta+delta_upper
    else:
        return 0
def mixing_lower_bounds():
    """
    Performs a similar function to finding the upper bounds on mixing, 
    Only on the lower bounds
    """
    if delta<0:
        return delta+delta_upper
    elif delta>0:
        return delta-delta_lower
    else:
        return 0

#Error propagation for symmetric quantities:
def dBdE():
  #"""
  #Uncertainty in B with respect to gamma
  #ray energy
  #"""
    l=set_multipolarity()
    if l=='E2(M1)':
      l=2
    return round((-B(tau)/E_g*(2*l+1)*E_g_error)**2,3)

def dBdI():
    """
    Uncertainty in B with respect to gamma
    ray intensity
    """
    l=set_multipolarity()
    if l=='E2(M1)':
      l=2
    return round((B(tau)/I_g*I_g_error)**2,3)**2
def dBdalpha():
    """
    Uncertainty in B with respect to internal
    conversion coefficient
    """
    l=set_multipolarity()
    if l=='E2(M1)':
      l=2
    return round((-B(tau)/(1+alpha_conv)*alpha_conv_error)**2,3)


"""
Asymmetric error is calculated via a 'consistent addition
technique' where B is calculated from the 'highest' value
and then subtracting the nominal value, etc
"""

def dBdtau_up():
   """
   Calculation of B for the longest lifetime, 
   for use in error propagation
   """
   l=set_multipolarity()
   if l=='E2(M1)':
      l=2
   return round(B(tau_down),3)
def dBdtau_down():
    """
    Calculation of B for the shortest lifetime, 
    for use in error propagation
    """
    l=set_multipolarity()
    if l=='E2(M1)':
      l=2
    return round(B(tau_up),3)
def uncertainty_tau_upper():
    return round((-B(tau)+dBdtau_up())**2,3)
def uncertainty_tau_lower():
    return round((B(tau)-dBdtau_down())**2,3)
  
#def calc_unc_delta_upper():
  #"""
  #This is an odd section, I need to calculate B under two 
  #delta conditions, upper and nominal, 
  #and subtract the two like I did for tau
  #"""
  #l=set_multipolarity()
  #if l=='E2(M1)':
    #tempB=B(tau)
    #delta=mixing_upper_bounds()
    #return -tempB+B(tau)
  #else:
    #return 0


#Aggregate uncertainty (upper bound)
def upper_uncertainty():
    """
    Returns the upper bound for final, added in quadrature
    uncertainty in B from any sources of uncertainty
    in measured quantities.
    """
    return round((dBdE()+dBdI()+dBdalpha()+uncertainty_tau_upper())**0.5,3)
#Aggregate uncertainty (lower bound)
def lower_uncertainty():
    """
    Returns the lower bound for final, added in quadrature
    uncertainty in B from any sources of uncertainty
    in measured quantities.
    """
    return round((dBdE()+dBdI()+dBdalpha()+uncertainty_tau_lower())**0.5,3)

#LaTeX Table header
output_file.write('\\begin{table}[ht]\n')
output_file.write('\\begin{tabular}{l|l|l|l|l|l|l}\n')
header1='E$_{lev}$ (keV) & E$_{\gamma}$ (keV) & I$_{\gamma}$ & ${\\tau}$ (fs)'
header2='               & $\pi\\ell$ & $\delta$            & B($\pi\\ell$) (W.u.)        '

#Terminal Outputs - Not LaTeX friendly
output_file.write(header1+header2+'\\\\\hline\hline\n')
for row in list(range(len(csvfile))):
    E_g=csvfile[row]['E_g']
    E_g_error=csvfile[row]['E_g_error']
    I_g=csvfile[row]['I_g']
    I_g_error=csvfile[row]['I_g_error']
    I_tot=csvfile[row]['I_g_total']
    delta=csvfile[row]['delta_mixing']
    delta_upper=csvfile[row]['delta_upper']
    delta_lower=csvfile[row]['delta_lower']
    tau=csvfile[row]['tau']
    tau_up=csvfile[row]['tau_up']-tau
    tau_down=tau-csvfile[row]['tau_down']
    alpha_conv=csvfile[row]['alpha_conv']
    alpha_conv_error=csvfile[row]['alpha_conv_error']
    A=csvfile[row]['A']
    multipolarity=csvfile[row]['multipolarity'].decode("utf-8")
    E_lev=csvfile[row]['E_level']
    #print('mixing',calc_unc_delta_upper(),tempB)
    lineEnergy=str(round(E_lev,2)).ljust(16,' ')+'& '+(str(round(E_g*1000,2))+' ('+str(int(E_g_error*1000))+')').ljust(19,' ')+'& '
    lineIntensity=(str(round(I_g,1))+' ('+str(int(I_g_error*10))+')').ljust(13,' ')+'& '+(str(int(tau))+'$^{+'+str(tau_up+tau)+'}_{'+str(tau_down-tau)+'}$').ljust(28,' ')+'& '
    lineLifetime=str(multipolarity).ljust(10,' ')+'& '
    lineDelta=(str(delta)+' $^{+'+str(delta_upper)+'}_{-'+str(delta_lower)+'}$').ljust(20,' ')+'& '
    lineMult=(str(round(B(tau),2))+' $^{+'+str(round(-upper_uncertainty()+B(tau),2))+'}_{'+str(round(B(tau)-lower_uncertainty(),2))+'}$ '+latex_friendly_units()).ljust(30,' ')+'\\\\ \n'
    output_file.write(lineEnergy+lineIntensity+lineLifetime+lineDelta+lineMult)
    print('B('+multipolarity+')=',B(tau),'p\m',upper_uncertainty(),latex_friendly_units(),'for the',E_g*1000,'keV transition leaving the',E_lev,'keV state')
    
    
output_file.write('\\end{tabular}\n')
output_file.write('\caption{REMEMBER TO CHANGE TABLE CAPTION AND REFERENCE TAG HERE! \label{tab:BE2}}\n')
output_file.write('\\end{table}')
output_file.close()
