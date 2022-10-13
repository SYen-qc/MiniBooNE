import math
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib.pyplot as plt

#variable definition
ml = 0.105658
Br_e = 1.6 * 1e-5
Br_mu = 0.636
POT = 18.75 * 1e20
me = 0.510999 * 1e-3
mmu = 105.658 * 1e-3
mK = 493.677 * 1e-3
mN = 35 * 1e-3
ma = 10 * 1e-3
x_l = mmu/mK
x_N = mN/mK
U_l4_sq = (1.59*1e-4)**2
rho = U_l4_sq * ((x_N**2 + x_l**2 - (x_N**2 - x_l**2)**2) * np.sqrt((1 - (x_N + x_l)**2) * (1 - (x_N - x_l)**2)))/(x_l**2 * (1-x_l**2)**2)
AMB = np.pi * (500)**2
L1 = 500/(1.9733 * 1e-16)
L2 = 510/(1.9733 * 1e-16)
Delta_L = 10/(1.9733 * 1e-16)
cN = 0.2
ce = 1
fa = 50
gagg = 2.32 * 1e-3 * (ce/fa)
ΓN = ((cN**2 * U_l4_sq * mN**3)/(128*math.pi*fa**2)) * (1 - (ma**2/mN**2))**2
pa_min = np.sqrt((1 + 0.97)/(1 - 0.97)) * np.sqrt(ma**2 - 4 * me**2)
Pnu0 = (0.493677**2 - 0.105658**2)/(2 * 0.493677)
PN0 = math.sqrt((mK**2 - (mN - ml)**2)*(mK**2 - (mN + ml)**2)) / (2*mK)
EN0 = math.sqrt(PN0**2 + mN**2)

tv = 500/(299792458)
Δtmax = 8*1e-9

def pNmin(x):
    return x*(tv/Δtmax)/(np.sqrt(1+2*(tv/Δtmax)))
print('pNmin')
print(pNmin(mN))

#eff
eff = np.loadtxt('/Users/shihyentseng/Documents/ＭiniBooNE/eff.txt')
effnew = np.empty([5011,2])

feff = interp1d(eff[:,0], eff[:,1], kind='linear')

xeff = np.linspace(0, 5.01, num=5010, endpoint=True)

for i in range(5011):
   effnew[i,0] = i * 1/1000
   effnew[i,1] = feff(i * 1/1000)
feffnew = interp1d(effnew[:,0], effnew[:,1], kind='linear')

#p_a distribution
p_a = np.loadtxt('/Users/shihyentseng/Documents/ＭiniBooNE/pa_MN35_Ma10.txt')
p_anew = np.empty([5501,2])
count_p_a = len(p_a[:,0])

for i in range(count_p_a):
    p_a[i,0] = p_a[i,0]/1000

fp_a = interp1d(p_a[:,0], p_a[:,1], kind='linear')

for i in range(5501):
   p_anew[i,0] = i * 1/1000
   p_anew[i,1] = fp_a(i * 1/1000)
   #p_anew[i,1] = fp_a(i * 1/1000)/pa_total
fp_a_fin = InterpolatedUnivariateSpline(p_anew[:,0], p_anew[:,1], k=5)

#heavy neutrino flux

X = np.loadtxt('/Users/shihyentseng/Documents/ＭiniBooNE/flux_fig29_v6.txt')
X2 = np.loadtxt('/Users/shihyentseng/Documents/ＭiniBooNE/fig_29_2.txt')
Y = np.loadtxt('/Users/shihyentseng/Documents/ＭiniBooNE/flux_fig31_antinu.txt')
Y2 = np.loadtxt('/Users/shihyentseng/Documents/ＭiniBooNE/flux_fig31_antinu_2.txt')

count_nu_mode = len(X[:,0])
count_nu_mode_2 = len(X2[:,0])
count_antinu_mode = len(Y[:,0])
count_antinu_mode_2 = len(Y2[:,0])

#Kaon
Xnew = np.empty([10001,2])
X2new = np.empty([10001,2])
Xtotal = np.empty([10001,2])
Ynew = np.empty([10001,2])
Y2new = np.empty([10001,2])
Ytotal = np.empty([10001,2])
#neutrino mode
NfromKfwd = np.empty([10001,2])
NfromKbwd = np.empty([10001,2])
fNfromKfwdnew = np.empty([10001,2])
fNfromKbwdnew = np.empty([10001,2])
pKfwd = np.empty([10001])
pKbwd = np.empty([10001])
ffwd = np.empty([10001])
fbwd = np.empty([10001])
Phi_N = np.empty([10001,2])
ffwdnewfin = np.empty([10001,2])
fbwdnewfin = np.empty([10001,2])
#antineutrino mode
NfromKfwdanti = np.empty([10001,2])
NfromKbwdanti = np.empty([10001,2])
fNfromKfwdantinew = np.empty([10001,2])
fNfromKbwdantinew = np.empty([10001,2])
pKfwdanti = np.empty([10001])
pKbwdanti = np.empty([10001])
ffwdanti = np.empty([10001])
fbwdanti = np.empty([10001])
ffwdnewantifin = np.empty([10001,2])
fbwdnewantifin = np.empty([10001,2])
'''
Kaon flux
'''
for i in range(count_nu_mode):
    X[i,1] = X[i,1] * (1/mK) * (2 * Pnu0) * (1 + (Pnu0**2/X[i,0]**2))**-1
for i in range(count_nu_mode):
    X[i,0] = (mK/2) * (X[i,0]/Pnu0 - Pnu0/X[i,0])

for i in range(count_nu_mode_2):
    X2[i,1] = X2[i,1] * (1/mK) * (2 * Pnu0) * (1 + (Pnu0**2/X2[i,0]**2))**-1
for i in range(count_nu_mode_2):
    X2[i,0] = (mK/2) * (X2[i,0]/Pnu0 - Pnu0/X2[i,0])

for i in range(count_antinu_mode):
    Y[i,1] = Y[i,1] * (1/mK) * (2 * Pnu0) * (1 + (Pnu0**2/Y[i,0]**2))**-1
for i in range(count_antinu_mode):
    Y[i,0] = (0.493677/2) * (Y[i,0]/((0.493677**2 - 0.105658**2)/(2 * 0.493677)) - ((0.493677**2 - 0.105658**2)/(2 * 0.493677))/Y[i,0])

for i in range(count_antinu_mode_2):
    Y2[i,1] = Y2[i,1] * (1/mK) * (2 * Pnu0) * (1 + (Pnu0**2/Y2[i,0]**2))**-1
for i in range(count_antinu_mode_2):
    Y2[i,0] = (mK/2) * (Y2[i,0]/Pnu0 - Pnu0/Y2[i,0])

f1 = interp1d(X[:,0], X[:,1], kind='linear')
f1_2 = interp1d(X2[:,0], X2[:,1], kind='linear')
f2 = interp1d(Y[:,0], Y[:,1], kind='linear')
f2_2 = interp1d(Y2[:,0], Y2[:,1], kind='linear')

xnew = np.linspace(0, 10.0, num=10000, endpoint=True)

for i in range(10001):
    Xnew[i,0] = i/1000
    X2new[i,0] = i/1000
    Ynew[i,0] = i/1000
    Y2new[i,0] = i/1000
    if Xnew[i,0] < np.min(X[:,0]) or Xnew[i,0] > np.max(X[:,0]):
        Xnew[i,1] = 0
    else: Xnew[i,1] = f1(i/1000)
    if X2new[i,0] < np.min(X2[:,0]) or X2new[i,0] > np.max(X2[:,0]):
        X2new[i,1] = 0
    else: X2new[i,1] = f1_2(i/1000)
    if Ynew[i,0] < np.min(Y[:,0]) or Ynew[i,0] > np.max(Y[:,0]):
        Ynew[i,1] = 0
    else: Ynew[i,1] = f2(i/1000)
    if Y2new[i,0] < np.min(Y2[:,0]) or Y2new[i,0] > np.max(Y2[:,0]):
        Y2new[i,1] = 0
    else: Y2new[i,1] = f2_2(i/1000)

for i in range(10001):
    Xtotal[i,0] = i/1000
    Ytotal[i,0] = i/1000
    Xtotal[i,1] = Xnew[i,1] + X2new[i,1]
    Ytotal[i,1] = Ynew[i,1] + Y2new[i,1]

f1fin = InterpolatedUnivariateSpline(Xtotal[:,0], Xtotal[:,1], k=1)
f2fin = InterpolatedUnivariateSpline(Ytotal[:,0], Ytotal[:,1], k=1)
#print("Kaon flux")
#print(f1fin(5.3))
for i in range(10001):
    if 0 < i/1000 and i/1000 < 0.1:
        if Xtotal[i-1,1] > Xtotal[i,1] and Xtotal[i,1] < Xtotal[i+1,1]:
            index_nu_mode = i
        if Xtotal[i-1,1] < Xtotal[i,1] and Xtotal[i,1] > Xtotal[i+1,1]:
            index_nu_mode_max = i


for i in range(10001):
    if 0 < Ynew[i,0] and Ynew[i,0] < 0.1:
        if Ytotal[i-1,1] > Ytotal[i,1] and Ytotal[i,1] < Ytotal[i+1,1]:
            index_antinu_mode = i
        if Ytotal[i-1,1] < Ytotal[i,1] and Ytotal[i,1] > Ytotal[i+1,1]:
            index_antinu_mode_max = i

'''
Heavy neutrino flux
'''
#neutrino mode
for i in range(10001):
    NfromKfwd[i,1] = Xtotal[i,1] * ((1/mK) * (EN0 + (Xnew[i,0]/math.sqrt(Xnew[i,0]**2 + mK**2) * PN0)))**-1
for i in range(10001):
    NfromKfwd[i,0] = (Xnew[i,0]/mK) * EN0 + (math.sqrt(Xnew[i,0]**2 + mK**2)/mK) * PN0
for i in range(10001):
    if NfromKfwd[i,0] > PN0:
        NfromKfwd[i,1] = NfromKfwd[i,1]
    else: NfromKfwd[i,1] = 0.0

for i in range(10001):
    NfromKbwd[i,1] = Xtotal[i,1] * ((1/mK) * (EN0 - (Xnew[i,0]/math.sqrt(Xnew[i,0]**2 + mK**2) * PN0)))**-1
for i in range(10001):
    NfromKbwd[i,0] = (Xnew[i,0]/mK) * EN0 - (math.sqrt(Xnew[i,0]**2 + mK**2)/mK) * PN0
for i in range(10001):
    if NfromKbwd[i,0] > 0.0:
        NfromKbwd[i,1] = NfromKbwd[i,1]
    else: NfromKbwd[i,1] = 0.0

for i in range(10001):
    if NfromKfwd[i,0] > 0.0:
        pKfwd[i] = (mK/mN**2) * (-math.sqrt(mN**2 + NfromKfwd[i,0]**2) * PN0 + EN0 * NfromKfwd[i,0])
        ffwd[i] = (pKfwd[i] * EN0 + math.sqrt(mK**2 + pKfwd[i]**2) * PN0)/(PN0 * (pKfwd[i] + math.sqrt(mK**2 + pKfwd[i]**2)))
    if NfromKbwd[i,0] > 0.0:
        pKbwd[i] = (mK/mN**2) * (math.sqrt(mN**2 + NfromKbwd[i,0]**2) * PN0 + EN0 * NfromKbwd[i,0])
        fbwd[i] = (pKbwd[i] * EN0 - math.sqrt(mK**2 + pKbwd[i]**2) * PN0)/(PN0 * (pKbwd[i] + math.sqrt(mK**2 + pKbwd[i]**2)))

fNfromKfwd = InterpolatedUnivariateSpline(NfromKfwd[:,0], NfromKfwd[:,1], k=3)
fNfromKbwd = InterpolatedUnivariateSpline(NfromKbwd[:,0], NfromKbwd[:,1], k=3)
ffwdnew = interp1d(NfromKfwd[:,0], ffwd[:], kind='cubic', axis = 0, fill_value="extrapolate")
fbwdnew = interp1d(NfromKbwd[:,0], fbwd[:], kind='linear', axis = 0, fill_value="extrapolate")

for i in range(10001): 
    ffwdnewfin[i,0] = i/1000
    if ffwdnewfin[i,0] < PN0:
        ffwdnewfin[i,1] = 0
    else: ffwdnewfin[i,1] = ffwdnew(i/1000)
    fbwdnewfin[i,1] = fbwdnew(i/1000)
    if ffwdnewfin[i,0] < np.min(NfromKfwd[:,0]) or ffwdnewfin[i,0] > np.max(NfromKfwd[:,0]):
        fNfromKfwdnew[i,1] = 0
    else: fNfromKfwdnew[i,1] = fNfromKfwd(i/1000)
    if ffwdnewfin[i,0] < np.min(NfromKbwd[:,0]) or ffwdnewfin[i,0] > np.max(NfromKbwd[:,0]):
        fNfromKbwdnew[i,1] = 0
    else: fNfromKbwdnew[i,1] = fNfromKbwd(i/1000)

#antineutrino mode
for i in range(10001):
    NfromKfwdanti[i,1] = Ytotal[i,1] * ((1/mK) * (EN0 + (Ynew[i,0]/math.sqrt(Ynew[i,0]**2 + mK**2) * PN0)))**-1
for i in range(10001):
    NfromKfwdanti[i,0] = (Ynew[i,0]/mK) * EN0 + (math.sqrt(Ynew[i,0]**2 + mK**2)/mK) * PN0
for i in range(10001):
    if NfromKfwdanti[i,0] > PN0:
        NfromKfwdanti[i,1] = NfromKfwdanti[i,1]
    else: NfromKfwdanti[i,1] = 0.0

for i in range(10001):
    NfromKbwdanti[i,1] = Ytotal[i,1] * ((1/mK) * (EN0 - (Ynew[i,0]/math.sqrt(Ynew[i,0]**2 + mK**2) * PN0)))**-1
for i in range(10001):
    NfromKbwdanti[i,0] = (Ynew[i,0]/mK) * EN0 - (math.sqrt(Ynew[i,0]**2 + mK**2)/mK) * PN0
for i in range(10001):
    if NfromKbwdanti[i,0] > 0.0:
        NfromKbwdanti[i,1] = NfromKbwdanti[i,1]
    else: NfromKbwdanti[i,1] = 0.0

for i in range(10001):
    if NfromKfwdanti[i,0] > 0.0:
        pKfwdanti[i] = (mK/mN**2) * (-math.sqrt(mN**2 + NfromKfwdanti[i,0]**2) * PN0 + EN0 * NfromKfwdanti[i,0])
        ffwdanti[i] = (pKfwdanti[i] * EN0 + math.sqrt(mK**2 + pKfwdanti[i]**2) * PN0)/(PN0 * (pKfwdanti[i] + math.sqrt(mK**2 + pKfwdanti[i]**2)))
    if NfromKbwdanti[i,0] > 0.0:
        pKbwdanti[i] = (mK/mN**2) * (math.sqrt(mN**2 + NfromKbwdanti[i,0]**2) * PN0 + EN0 * NfromKbwdanti[i,0])
        fbwdanti[i] = (pKbwdanti[i] * EN0 - math.sqrt(mK**2 + pKbwdanti[i]**2) * PN0)/(PN0 * (pKbwdanti[i] + math.sqrt(mK**2 + pKbwdanti[i]**2)))

fNfromKfwdanti = InterpolatedUnivariateSpline(NfromKfwdanti[:,0], NfromKfwdanti[:,1], k=3)
fNfromKbwdanti = InterpolatedUnivariateSpline(NfromKbwdanti[:,0], NfromKbwdanti[:,1], k=3)
ffwdantinew = interp1d(NfromKfwdanti[:,0], ffwdanti[:], kind='cubic', axis = 0, fill_value="extrapolate")
fbwdantinew = interp1d(NfromKbwdanti[:,0], fbwdanti[:], kind='linear', axis = 0, fill_value="extrapolate")

for i in range(10001): 
    ffwdnewantifin[i,0] = i/1000
    if ffwdnewantifin[i,0] < PN0:
        ffwdnewantifin[i,1] = 0
    else: ffwdnewantifin[i,1] = ffwdantinew(i/1000)
    fbwdnewantifin[i,1] = fbwdantinew(i/1000)
    if ffwdnewantifin[i,0] < np.min(NfromKfwdanti[:,0]) or ffwdnewantifin[i,0] > np.max(NfromKfwdanti[:,0]):
        fNfromKfwdantinew[i,1] = 0
    else: fNfromKfwdantinew[i,1] = fNfromKfwdanti(i/1000)
    if ffwdnewantifin[i,0] < np.min(NfromKbwdanti[:,0]) or ffwdnewantifin[i,0] > np.max(NfromKbwdanti[:,0]):
        fNfromKbwdantinew[i,1] = 0
    else: fNfromKbwdantinew[i,1] = fNfromKbwdanti(i/1000)

xNfromK = np.linspace(0, 10.0, num=10000, endpoint=True)

flux_neu = InterpolatedUnivariateSpline(ffwdnewfin[:,0], (ffwdnewfin[:,1]*fNfromKfwdnew[:,1]+fbwdnewfin[:,1]*fNfromKbwdnew[:,1])/(ffwdnewfin[:,1]+fbwdnewfin[:,1]), k=1)
flux_antineu = InterpolatedUnivariateSpline(ffwdnewfin[:,0], (ffwdnewantifin[:,1]*fNfromKfwdantinew[:,1]+fbwdnewantifin[:,1]*fNfromKbwdantinew[:,1])/(ffwdnewantifin[:,1]+fbwdnewantifin[:,1]), k=1)
#print("phi_N")
#print(flux_neu(1))
phi_neu = np.empty([10001,2])
phi_antineu = np.empty([10001,2])

for i in range(10001):
   phi_neu[i,0] = i * 1/1000
   phi_neu[i,1] = flux_neu(i * 1/1000)
   phi_antineu[i,0] = i * 1/1000
   phi_antineu[i,1] = flux_antineu(i * 1/1000)

#angular distribution
#functoins
def EN(x):
    return np.sqrt(mN**2+x**2)

def phi_neu_fun(x):
    return flux_neu(x)

def phi_antineu_fun(x):
    return flux_antineu(x)

def P_N_dec(x):
    return np.exp(-L1*ΓN*mN/x) - np.exp(-L2*ΓN*mN/x)

def w_time(x):
    t0 = (500/299792458)
    dt = 1.6 * 1e-6
    tN = EN(x) * t0 / x
    if tN < dt + t0:
       return (t0 + dt - tN)/dt
    if tN >= dt + t0:
       return 0

def P_a_dec(x):
    Γ_aee = ((ce**2*me**2*ma)/(8*np.pi*fa**2)) * np.sqrt(1 - (4*me**2)/(ma**2))
    Γ_agg = (gagg**2 * ma**3)/(64*np.pi)
    Γ_a = Γ_aee + Γ_agg
    return 1 - np.exp(-Delta_L * Γ_a * (ma/x))

def ΓN_lab(x):
    return (mN/EN(x)) * ((cN**2 * U_l4_sq * mN**3)/(128*np.pi*fa**2)) * (1 - (ma**2/mN**2))**2

def R(x):
    return (mN**2 - ma**2)/(2 * ma * x)

def pa_plus(x,z):
    return ((mN**2 + ma**2) * x * z + EN(x) * np.sqrt((mN**2 - ma**2)**2 - 4*ma**2*x**2*(1 - z**2)))/(2 * (np.square(EN(x)) - x**2 * z**2))

def pa_minus(x,z):
    return ((mN**2 + ma**2) * x * z - EN(x) * np.sqrt((mN**2 - ma**2)**2 - 4*ma**2*x**2*(1 - z**2)))/(2 * (np.square(EN(x)) - x**2 * z**2))

def eff_fun(x):
    return feffnew(x)

def p_a_fun(x):
    return fp_a_fin(x)

def dΓ_dz(x,z):
    return ((cN**2 * U_l4_sq)/(128*np.pi*fa**2)) * ((mN**2 * (mN**2 - ma**2) * np.square(pa_plus(x,z)))/(EN(x) * np.abs(pa_plus(x,z) * EN(x) - x * np.sqrt(ma**2+np.square(pa_plus(x,z))) * z)))

def integrand_neu(x,z):
    if ((mN**2 - ma**2)**2 - 4*ma**2*x**2*(1 - z**2)) >= 0:
        if R(x) >= 1:
            if pa_plus(x,z) >= pa_min and pa_minus(x,z) < pa_min:
            #if 1.25 >= pa_plus(x,z) >= pa_min and pa_minus(x,z) < pa_min:
                if 1 >= z >= -1:
                    return phi_neu_fun(x)*P_N_dec(x)*w_time(x)*P_a_dec(pa_plus(x,z)) * (1/ΓN_lab(x)) * ((cN**2 * U_l4_sq)/(128*np.pi*fa**2)) * ((mN**2 * (mN**2 - ma**2) * np.square(pa_plus(x,z)))/(EN(x) * np.abs(pa_plus(x,z) * EN(x) - x * np.sqrt(ma**2+np.square(pa_plus(x,z))) * z))) * eff_fun(pa_plus(x,z)) * p_a_fun(pa_plus(x,z))
                    #return phi_antineu_fun(x)*P_N_dec(x)*w_time(x)*P_a_dec(pa_plus(x,z)) * (1/ΓN_lab(x)) * ((cN**2 * U_l4_sq)/(128*np.pi*fa**2)) * ((mN**2 * (mN**2 - ma**2) * np.square(pa_plus(x,z)))/(EN(x) * np.abs(pa_plus(x,z) * EN(x) - x * np.sqrt(ma**2+np.square(pa_plus(x,z))) * z))) * eff_fun(pa_plus(x,z)) * p_a_fun(pa_plus(x,z))
                else:
                    return 0
            elif pa_plus(x,z) < pa_min:
                return 0
            else:
                return 0
        elif R(x) < 1:
            if pa_plus(x,z) >= pa_min and pa_minus(x,z) < pa_min:
            #if 1.25 >= pa_plus(x,z) >= pa_min and pa_minus(x,z) < pa_min:
                if 1 >= z >= np.sqrt(1 - R(x)**2) or -np.sqrt(1 - R(x)**2) >= z >= -1:
                    return phi_neu_fun(x)*P_N_dec(x)*w_time(x)*P_a_dec(pa_plus(x,z)) * (1/ΓN_lab(x)) * ((cN**2 * U_l4_sq)/(128*np.pi*fa**2)) * ((mN**2 * (mN**2 - ma**2) * np.square(pa_plus(x,z)))/(EN(x) * np.abs(pa_plus(x,z) * EN(x) - x * np.sqrt(ma**2+np.square(pa_plus(x,z))) * z))) * eff_fun(pa_plus(x,z)) * p_a_fun(pa_plus(x,z))
                    #return phi_antineu_fun(x)*P_N_dec(x)*w_time(x)*P_a_dec(pa_plus(x,z)) * (1/ΓN_lab(x)) * ((cN**2 * U_l4_sq)/(128*np.pi*fa**2)) * ((mN**2 * (mN**2 - ma**2) * np.square(pa_plus(x,z)))/(EN(x) * np.abs(pa_plus(x,z) * EN(x) - x * np.sqrt(ma**2+np.square(pa_plus(x,z))) * z))) * eff_fun(pa_plus(x,z)) * p_a_fun(pa_plus(x,z))
                else:
                    return 0
            elif pa_plus(x,z) < pa_min:
                return 0
            else:
                return 0
        elif R(x) < 1:
            if pa_plus(x,z) >= pa_min and pa_minus(x,z) >= pa_min:
            #if 1.25 >= pa_plus(x,z) >= pa_min and 1.25 >= pa_minus(x,z) >= 0.15:
                if 1 >= z >= np.sqrt(1 - R(x)**2) or -np.sqrt(1 - R(x)**2) >= z >= -1:
                    return phi_neu_fun(x)*P_N_dec(x)*w_time(x)*P_a_dec(pa_plus(x,z)) * (1/ΓN_lab(x)) * ((cN**2 * U_l4_sq)/(128*np.pi*fa**2)) * ((mN**2 * (mN**2 - ma**2) * np.square(pa_plus(x,z)))/(EN(x) * np.abs(pa_plus(x,z) * EN(x) - x * np.sqrt(ma**2+np.square(pa_plus(x,z))) * z))) * eff_fun(pa_plus(x,z)) * p_a_fun(pa_plus(x,z)) + phi_neu_fun(x)*P_N_dec(x)*w_time(x)*P_a_dec(pa_minus(x,z)) * (1/ΓN_lab(x)) * ((cN**2 * U_l4_sq)/(128*np.pi*fa**2)) * ((mN**2 * (mN**2 - ma**2) * np.square(pa_minus(x,z)))/(EN(x) * np.abs(pa_minus(x,z) * EN(x) - x * np.sqrt(ma**2+np.square(pa_minus(x,z))) * z))) * eff_fun(pa_minus(x,z)) * p_a_fun(pa_minus(x,z))
                    #return phi_antineu_fun(x)*P_N_dec(x)*w_time(x)*P_a_dec(pa_plus(x,z)) * (1/ΓN_lab(x)) * ((cN**2 * U_l4_sq)/(128*np.pi*fa**2)) * ((mN**2 * (mN**2 - ma**2) * np.square(pa_plus(x,z)))/(EN(x) * np.abs(pa_plus(x,z) * EN(x) - x * np.sqrt(ma**2+np.square(pa_plus(x,z))) * z))) * eff_fun(pa_plus(x,z)) * p_a_fun(pa_plus(x,z))
                else:
                    return 0
            elif pa_plus(x,z) < pa_min:
                return 0
            else:
                return 0
        else:
            return 0
    elif ((mN**2 - ma**2)**2 - 4*ma**2*x**2*(1 - z**2)) < 0:
        return 0
    else:
        return 0

def N(x,z):
    if (mN**2 - ma**2)**2 - 4*ma**2*x**2*(1 - z**2) > 0:
        if R(x) >= 1:
            if pa_plus(x,z) >= pa_min and pa_minus(x,z) < pa_min:
            #if 1.25 >= pa_plus(x,z) >= pa_min and pa_minus(x,z) < pa_min:
                if 1 >= z >= -1:
                    return p_a_fun(pa_plus(x,z))
                else:
                    return 0
            elif pa_plus(x,z) < pa_min:
                return 0
            else:
                return 0
        elif R(x) < 1:
            if pa_plus(x,z) >= pa_min and pa_minus(x,z) < pa_min:
            #if 1.25 >= pa_plus(x,z) >= pa_min and pa_minus(x,z) < pa_min:
                if 1 >= z >= np.sqrt(1 - R(x)**2) or -np.sqrt(1 - R(x)**2) >= z >= -1:
                    return p_a_fun(pa_plus(x,z))
                else:
                    return 0
            elif pa_plus(x,z) < pa_min:
                return 0
            else:
                return 0
        elif R(x) < 1:
            if pa_plus(x,z) >= pa_min and pa_minus(x,z) >= pa_min:
            #if 1.25 >= pa_plus(x,z) >= pa_min and 1.25 >= pa_minus(x,z) >= pa_min:
                if 1 >= z >= np.sqrt(1 - R(x)**2) or -np.sqrt(1 - R(x)**2) >= z >= -1:
                    return p_a_fun(pa_plus(x,z)) + p_a_fun(pa_minus(x,z))
                else:
                    return 0
            elif pa_plus(x,z) < pa_min:
                return 0
            else:
                return 0
        else:
            return 0
    elif (mN**2 - ma**2)**2 - 4*ma**2*x**2*(1 - z**2) <= 0:
        return 0
    else:
        return 0

xnew = np.arange(-0.9, 0.9, 0.2)
partial_int_neu = lambda z: integrate.quad(integrand_neu, PN0, 5, args=(z,), limit=60, full_output=1)[0]
partial_N = lambda z: integrate.quad(N, PN0, 5, args=(z,), limit=60, full_output=1)[0]

def result(z):
    if  partial_N(z) > 0:
        return POT*rho*AMB*(partial_int_neu(z))/partial_N(z)
    else:
        return 0

integration_bin = np.empty([20,2])
for i in range(20):
    integration_bin[i,0] = -0.95 + i/10
    integration_bin[i,1] = integrate.quad(lambda z:result(z), -1 + i/10, -1 + (i+1)/10, limit=100, full_output=1)[0]

print("result")
print(integration_bin[:,0])
print(integration_bin[:,1])


sum = 0
for i in range(20):
    sum = sum + integration_bin[i,1]
print(sum)

angular_data = np.loadtxt('/Users/shihyentseng/Documents/ＭiniBooNE/angular_data.txt')
angular_data_array = np.empty([20,2])
for i in range(20):
    angular_data_array[i,0] = -0.95 + i/10
    angular_data_array[i,1] = angular_data[i,1]

angular_sim = np.loadtxt('/Users/shihyentseng/Documents/ＭiniBooNE/angular_sim.txt')
angular_sim_array = np.empty([22,2])
angular_sim_array[0,0] = -1.0
angular_sim_array[0,1] = angular_sim[0,1]
angular_sim_array[21,0] = 1.0
angular_sim_array[21,1] = angular_sim[19,1]
for i in range(20):
    angular_sim_array[i+1,0] = -0.95 + i/10
    angular_sim_array[i+1,1] = angular_sim[i,1]

angular_excess = np.empty([22,2])
angular_excess[0,0] = -1.0
angular_excess[0,1] = integration_bin[0,1] + angular_sim[0,1]
angular_excess[21,0] = 1.0
angular_excess[21,1] = integration_bin[19,1] + angular_sim[19,1]
for i in range(20):
    angular_excess[i+1,0] = -0.95 + i/10
    angular_excess[i+1,1] = integration_bin[i,1] + angular_sim[i,1]

plt.figure(figsize=(10.0, 9.0))
#main results
plt.scatter(angular_data_array[:,0], angular_data_array[:,1],c='k')
plt.fill_between(angular_sim_array[:,0], angular_sim_array[:,1] + 1, step="mid", alpha=0.5)
plt.step(angular_excess[:,0], angular_excess[:,1] + 1, where='mid', c='k', linestyle='dashed')

plt.xticks(fontsize=20, rotation=0)
plt.yticks(fontsize=20, rotation=0)
plt.xlabel('cos$\\theta$', fontsize=25)
plt.ylabel('Events', fontsize=25)
plt.xlim(-1, 1)
plt.ylim(0, 1000)
plt.grid(linestyle='dotted')
#plt.savefig('/Users/shihyentseng/Documents/ＭiniBooNE/angular_dis_v6.pdf', dpi=300)
plt.show()

'''
data1 = np.array([ffwdnewfin[:,0], ffwdnewfin[:,1]*fNfromKfwdnew[:,1]+fbwdnewfin[:,1]*fNfromKbwdnew[:,1]])
data1 = data1.T #here you transpose your data, so to have it in two columns
#datafile_path = "/Users/shihyentseng/Documents/ＭiniBooNE/N_flux_250_nu_mode_result.txt"
datafile_path = "/Users/shihyentseng/Documents/ＭiniBooNE/N_flux_350_nu_mode_result.txt"
with open(datafile_path, 'w+') as datafile_id: #here you open the ascii file
    np.savetxt(datafile_id, data1, fmt=['%.8f','%.8e']) #here the ascii file is written.

data2 = np.array([ffwdnewfin[:,0], ffwdnewantifin[:,1]*fNfromKfwdantinew[:,1]+fbwdnewantifin[:,1]*fNfromKbwdantinew[:,1]])
data2 = data2.T 
#datafile_path = "/Users/shihyentseng/Documents/ＭiniBooNE/N_flux_250_antinu_mode_result.txt"
datafile_path = "/Users/shihyentseng/Documents/ＭiniBooNE/N_flux_350_antinu_mode_result.txt"
with open(datafile_path, 'w+') as datafile_id:
    np.savetxt(datafile_id, data2, fmt=['%.8f','%.8e'])
'''