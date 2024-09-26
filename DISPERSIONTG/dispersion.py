from DISPERSIONTG.decay import markersize, capsize, DecaySubPlot
import numpy as np
import matplotlib.pyplot as plt
import lmfit
import scipy.stats

gammaH = 267.52218744*10**6
gammaN = -27.116*10**6
gammaC = 67.2828*10**6
pi = 22/7
hbar = 6.626069*(10**(-34))/(2*pi)
mu0 = 4*pi*10**(-7)

# ====== Dispersion plots without fitting (for all nuclii) ==============
class DispPlotwithoutFit:
    def __init__(self,spin,decayfitout,dispfitflag,fittype,wR,MC,pcutoff,prefilename,ymax1):
        # If you don not need dispersion fit
        s1 = ToDisp(decayfitout)
        filename1 = prefilename + '_' + 'dispersion_withoutfit.pdf'
        ncol  = 5
        nrows = (s1.nres//ncol)+1
        nrf = s1.nrf
        self.fig1, self.axes1 = plt.subplots(nrows,ncol,figsize=(ncol*2.3,1.8*nrows))
        xmax = 1.15*np.amax(s1.effrf)
        xmin = -0.1*xmax
        if (ymax1 == -1 ):
            ymax = 1.15*np.amax(s1.r1rhocor,0)
        elif (ymax1 == 0):
            ymax = np.amax(s1.r1rhocor)
        else:
            ymax = ymax1
        ymin = -0.1*ymax
        k=0
        for i in range(nrows):
            for j in range(ncol):
                if (k < s1.nres):
                    #self.axes[i,j].errorbar(self.effrf1[:,k],s.r1rho[:,k], yerr=s.r1rhoerr[:,k],
                    #                        marker='o', color ='red',linestyle = '',markersize=markersize,
                    #                        ecolor = 'black',capsize =capsize,lw=capsize)
                    self.axes1[i,j].errorbar(np.array(s1.effrf)[:,k],np.array(s1.r1rhocor)[:,k],
                                    yerr=np.array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                    markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                    markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes1[i,j].set_xlim([-5.0,xmax])
                    self.axes1[i,j].set_ylim([ymin,ymax])
                    self.axes1[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    if (spin == '1H'):
                        self.axes1[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    DecaySubPlot(self.axes1[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                              s1.seqname[k])
                    k=k+1
                else:
                    self.axes1[i,j].set_axis_off()
                    k=k+1
        self.fig1.subplots_adjust(wspace=0.4, hspace=0.9)
        # save plots to pdf
        self.fig1.savefig(filename1,dpi=1200)





# ===== Decay fit results to dispersion fitting module  =======================
class ToDisp:
    def __init__(self, rffitlist):
        """
        Initializes the ToDisp class by processing the list of RF fit results.

        Parameters:
        -----------
        rffitlist : list
            List of RF fit results from which to extract the data.
        """
        self.nrf = len(rffitlist)  # Number of RF fit results
        self.nres = rffitlist[0].nres  # Number of residues (peak sequences)
        self.seqname = rffitlist[0].seqname  # Sequence names

        # Initialize empty lists to hold the extracted data
        effrf, angles, times, normvolumes, normnoises = [], [], [], [], []
        ntimes, r1rho, r1rhoerr = [], [], []

        # Loop over each RF fit result and extract relevant data
        for i in range(self.nrf):
            ntimes.append(len(rffitlist[i].timelist))  # Number of time points for each RF fit
            angles.append(rffitlist[i].angles)  # Extract angles
            effrf.append(rffitlist[i].effrfkhz)  # Extract effective RF frequencies
            times.append(rffitlist[i].timelist)  # Extract time lists
            normvolumes.append(rffitlist[i].normvolumes)  # Extract normalized volumes
            normnoises.append(rffitlist[i].normnoises)  # Extract normalized noise values
            r1rho.append(rffitlist[i].mcfitvalues[:, 1])  # Extract R1rho values
            r1rhoerr.append(rffitlist[i].mcfiterr[:, 1])  # Extract R1rho errors

        # Create empty arrays to store combined data across all RF fits
        normvolumes_combined = np.zeros((self.nres, sum(ntimes)))  # Combined normalized volumes
        normnoises_combined = np.zeros((self.nres, sum(ntimes)))  # Combined normalized noise values
        rfangletime_combined = np.zeros((self.nres, sum(ntimes), 3))  # Combined RF, angle, and time data

        # Loop over each residue and combine data across all RF fits
        for k in range(self.nres):
            combined_rf_angle_time = []  # Temporary list to store combined RF, angle, and time data
            combined_volumes = []  # Temporary list to store combined normalized volumes
            combined_noises = []  # Temporary list to store combined normalized noise values

            # Loop over each RF fit and extract data for each time point
            for i in range(self.nrf):
                for j in range(ntimes[i]):
                    # Combine RF, angle, and time data for each time point
                    combined_rf_angle_time.append([effrf[i][k], angles[i][j], times[i][j] * 0.001])
                    combined_volumes.append(normvolumes[i][k][j])  # Store normalized volume
                    combined_noises.append(normnoises[i][k][j])  # Store normalized noise

            # Store the combined data for the current residue
            rfangletime_combined[k, :, :] = np.array(combined_rf_angle_time)
            normvolumes_combined[k, :] = np.array(combined_volumes)
            normnoises_combined[k, :] = np.array(combined_noises)

        # Set class attributes
        self.times = [[x * 0.001 for x in y] for y in times]  # Convert times to seconds
        self.effrf = effrf  # Effective RF frequencies
        self.angles = angles  # Angles
        self.normvolumes = normvolumes  # Normalized volumes
        self.normnoises = normnoises  # Normalized noise values
        self.ntimes = ntimes  # Number of time points
        self.r1rho = r1rho  # R1rho values
        self.r1rhoerr = r1rhoerr  # Errors in R1rho values

        # Correct R1rho values and errors using angles
        self.r1rhocor = self.r1rho / (np.sin(self.angles) ** 2)
        self.r1rhoerrcor = self.r1rhoerr / (np.sin(self.angles) ** 2)

        # Store combined RF, angle, and time data, and normalized volumes/noises
        self.rfangletime = rfangletime_combined
        self.normvolumes1 = normvolumes_combined
        self.normnoises1 = normnoises_combined








# ====== All residues Dispersion fitting and right model estimation for Proton  ==============
class HDispFit:
    def __init__(self,rfouttodisp,fittype,wR,B0,rHHlist,dwHcsa,r20values,MC,pcutoff):
        s=rfouttodisp
        fitvalues1,fiterr1,mcfitvalues1,mcfiterr1 = [],[],[],[]
        reffitvalues1,reffiterr1,refmcfitvalues1,refmcfiterr1 = [],[],[],[]
        ftestresults,ysimref,ysim  = [],[],[]
        r1rhosimref,r1rhosim = [],[]
        x4 = np.arange(0.2,150.0,0.2)  # effrf singlet
        x6 = np.array([90.0 for i in x4]) #angle singlet
        if np.array(rHHlist).size == 1:
            rHH = [rHHlist for i in range(s.nres)]
        else:
            rHH = np.array(rHHlist)
        try:
            for k in range(s.nres):
                print(f"Residue {k} fit is running ...")
                s2=HDispFitsr(rfouttodisp,k,fittype,wR,B0,rHH[k],dwHcsa,r20values,MC,pcutoff)
                fitvalues1.append(s2.fitvalues1)
                fiterr1.append(s2.fiterr1)
                mcfitvalues1.append(s2.mcfitvalues1)
                mcfiterr1.append(s2.mcfiterr1)
                r1rhosim.append(s2.r1rhosim)
                reffitvalues1.append(s2.reffitvalues1)
                reffiterr1.append(s2.reffiterr1)
                refmcfitvalues1.append(s2.refmcfitvalues1)
                refmcfiterr1.append(s2.refmcfiterr1)
                ysimref.append(s2.ysimref)
                r1rhosimref.append(s2.r1rhosimref)
                ftestresults.append(s2.ftestresults)
            self.fitvalues1,self.fiterr1 = fitvalues1,fiterr1
            self.mcfitvalues1,self.mcfiterr1 = mcfitvalues1,mcfiterr1
            self.refmcfitvalues1,self.refmcfiterr1 = refmcfitvalues1,refmcfiterr1
            self.effrf = x4
            self.r1rhosim = np.array(r1rhosim)
            self.r1rhosimref = np.array(r1rhosimref)
            self.ftestresults = np.array(ftestresults)
            self.rHH = rHH
        except IndexError:
            print("rHH list is not same length as no. of peakks")



# ====== Single residueDispersion fitting and right model estimation for Proton ==============
class HDispFitsr:
    def __init__(self,rfouttodisp,k,fittype,wR,B0,rHH,dwHcsa,r20values,MC,pcutoff):
        s=rfouttodisp
        # s = experimental r1rho values
        #s2 = no dispersion
        #s3 = dispersion
        # x1, x2, x3 are triplets with rf, angle and time
        x1 = s.rfangletime[k,:,0] * np.sin(s.rfangletime[k,:,1])
        x2 = s.rfangletime[k,:,1]*180/np.pi # angle triplet
        x3 = s.rfangletime[k,:,2] # time triplet
        y  = s.normvolumes1[k,:]
        yerr = s.normnoises1[k,:]
        logy = np.log(abs(y))
        logyerr = abs(np.divide(yerr,y))
        tau0,s20,r200,a0,dwHex0 = 0.01,0.98,r20values[0],1,0.5
        taumin,s2min,r20min,a0min,dwHexmin = 0.00001,0.1,r20values[1],0.1,0.001
        taumax,s2max,r20max,a0max,dwHexmax = 0.5,0.999,r20values[2],10,20
        s2= NoDispDecayFit(x1,x2,x3,logy,logyerr,r200,r20min,r20values[3],a0,a0min,a0max,MC)
        # NERRD fitting
        if (fittype == "NERRD"):
            s3=HNERRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,rHH,dwHcsa,[tau0,taumin,taumax],
                                  [s20,s2min,s2max],[r200,r20min,r20max],[a0,a0min,a0max],MC)
        else:
            print("Fit type not currently supported")
        self.fitvalues1=s3.popt
        self.fiterr1=s3.perr
        self.mcfitvalues1=s3.mcfitvalues2
        self.mcfiterr1=s3.mcfiterr2
        self.r1rhosim=s3.tmpr1rhosim
        self.reffitvalues1=s2.poptref
        self.reffiterr1=s2.perrref
        self.refmcfitvalues1=s2.refmcfitvalues2
        self.refmcfiterr1=s2.refmcfiterr2
        self.ysimref=s2.tmpysimref
        self.r1rhosimref=s2.tmpr1rhosimref
        # Run F-test for the fit
        tmpftestresults = FRatioToPvalue(logy,logyerr,s3.tmpysim,s2.tmpysimref,len(s3.popt),
                                                 len(s2.poptref),pcutoff)
        self.ftestresults=tmpftestresults



# ===== Class for No Dispersion Decay Fit ============================
class NoDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,r200,r20min,r20max,a0,a0min,a0max,MC):
        x4 = np.arange(0.2,150.0,0.2)  # effrf singlet
        x6 = np.array([90.0 for i in x4]) #angle singlet
        fitmodel = lmfit.Model(NoDispDecay,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('R20', value = r200, min=r20min, max=r20max)
        fitmodel.set_param_hint('a', value = a0, min=a0min, max=a0max)
        params = fitmodel.make_params()
        #print(params)
        fitresult=fitmodel.fit(logy, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
        self.poptref =[fitresult.params['R20'].value,fitresult.params['a'].value]
        self.perrref = [fitresult.params['R20'].stderr,fitresult.params['a'].stderr]
        # Monte-corlo fitting of decays and rf fields together
        tmpfit = []
        for j in range(MC):
            tmpysimref = (NoDispDecay(x1,x2,x3,self.poptref[0],self.poptref[1]) +
                          np.random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(logy, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsimref =[fitresult.params['R20'].value,fitresult.params['a'].value]
                tmpfit.append(poptsimref)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.refmcfitvalues2=np.mean(tmpfit,0)
        self.refmcfiterr2=np.std(tmpfit,0)
        self.tmpysimref = NoDispDecay(x1,x2,x3,np.mean(tmpfit,0)[0],np.mean(tmpfit,0)[1])
        self.tmpr1rhosimref = NoDisp(x4,x6,np.mean(tmpfit,0)[0])


# ===== Proton r1rho decay (with No dispersion) calculations ===================
def NoDispDecay(w1,tilt,t,R20,a):
    tmp1 = R20*np.sin(tilt*np.pi/180)*np.sin(tilt*np.pi/180)
    return np.log(abs(a))-(t*tmp1)
# ===== Proton r1rho (with No dispersion) calculations =========================
def NoDisp(w1,tilt,R20):
    return R20*np.sin(tilt*np.pi/180)*np.sin(tilt*np.pi/180)


# ===== Fratio to P value estimatiion ==========================================
def FRatioToPvalue(y,yerr,ymodel1,ymodel2,dof1,dof2,pcutoff):
    n =len(y)
    dof1 = n-dof1
    dof2 = n-dof2
    meany = np.mean(y)
    variance1 = sum([((y[i]-meany)/yerr[i])**2 for i in range(n)])
    chi2model1 = sum([((y[i]-ymodel1[i])/yerr[i])**2 for i in range(n)])
    chi2model2 = sum([((y[i]-ymodel2[i])/yerr[i])**2 for i in range(n)])
    redchi2model1 = chi2model1/dof1
    redchi2model2 = chi2model2/dof2
    # r2 value close to one is good fit
    r2model1 = 1 - (chi2model1/variance1)
    r2model2 = 1 - (chi2model1/variance1)
    if (dof1 == dof2):
        fratio = chi2model2/chi2model1
    else:
        fratio = (((chi2model2-chi2model1)/(dof2-dof1))/(chi2model1/dof1))
    single_tailed_pval = 1-scipy.stats.f.cdf(fratio,dof2,dof1)
    double_tailed_pval = single_tailed_pval * 2
    if ((chi2model1 < chi2model2) & (single_tailed_pval < pcutoff)):
        modelnum = 1
    else :
        modelnum = 2
    return [modelnum,single_tailed_pval,double_tailed_pval,chi2model1,chi2model2,
            redchi2model1,redchi2model2,dof1,dof2,fratio]



# ===== Class for Proton NERRD Dispersion Decay Fit ============================
class HNERRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,rHH,dwcsa,tauvalues,s2values,r2values,a0values,MC):
        x4 = np.arange(0.2,150.0,0.2)  # effrf singlet
        x6 = np.array([90.0 for i in x4]) #angle singlet
        tau0,taumin,taumax = tauvalues[0],tauvalues[1],tauvalues[2]
        s20,s2min,s2max = s2values[0],s2values[1],s2values[2]
        r200,r20min,r20max=r2values[0],r2values[1],r2values[2]
        a0,a0min,a0max = a0values[0],a0values[1],a0values[2]
        fitmodel = lmfit.Model(ProtonR1rhoNERRDDecay1,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('tau', value = tau0, min=taumin, max=taumax)
        fitmodel.set_param_hint('S2', value = s20, min=s2min, max=s2max)
        fitmodel.set_param_hint('R20', value = r200, min=r20min, max=r20max)
        fitmodel.set_param_hint('a', value = a0, min=a0min, max=a0max)
        fitmodel.set_param_hint('B0',value=B0,vary=False)
        fitmodel.set_param_hint('wR',value=wR,vary=False)
        fitmodel.set_param_hint('rHH',value=rHH,vary=False)
        fitmodel.set_param_hint('dwcsa',value=dwcsa,vary=False)
        params = fitmodel.make_params()
        #print(params)
        fitresult=fitmodel.fit(logy, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
        self.popt =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                  fitresult.params['R20'].value,fitresult.params['a'].value]
        self.perr = [fitresult.params['tau'].stderr,fitresult.params['S2'].stderr,
                  fitresult.params['R20'].stderr,fitresult.params['a'].stderr]
        # Monte-corlo fitting of decays and rf fields together
        tmpfit = []
        for j in range(MC):
            tmpysimref = (ProtonR1rhoNERRDDecay1(x1,x2,x3,B0,wR,rHH,dwcsa,self.popt[0],
                                self.popt[1],self.popt[2],self.popt[3]) + np.random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                             fitresult.params['R20'].value,fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=np.mean(tmpfit,0)
        self.mcfiterr2=np.std(tmpfit,0)
        self.tmpysim =ProtonR1rhoNERRDDecay1(x1,x2,x3,B0,wR,rHH,dwcsa,
                            np.mean(tmpfit,0)[0],np.mean(tmpfit,0)[1],np.mean(tmpfit,0)[2],np.mean(tmpfit,0)[3])
        self.tmpr1rhosim = ProtonR1rhoNERRD1(x4,x6,B0,wR,rHH,dwcsa,np.mean(tmpfit,0)[0],
                            np.mean(tmpfit,0)[1],np.mean(tmpfit,0)[2])
        




# ===== Proton R1rho (all contributions) calculation ===========================
# for AA spin systems
def ProtonR1rho(wH,wN,rHH,rHN,dw,tau,S2,dwhex,pE,R1,R20,wR,w1,tilt):
    """ All frequencies are in kHz and times are in milliseconds
    wH = larmor frequency proton, wR= spinning frequency, w1 = spin lock field
    tilt = spinlock tilt angle, rHH  = proton-proton distances (in Angstroms)
    """
    pi=np.pi
    tau = tau*10**(-3)
    wH = -2*pi*wH*10**3
    wN = 2*pi*wN*10**3
    wR = 2*pi*wR*10**3
    tilt = tilt*pi/180
    dw = wH*dw*10**(-6)
    wE = (2*pi*w1*10**3)/(np.sin(tilt)+0.00001)
    # print(f"tau {tau} wH {wH} wN {wN} wR {wR} tilt {tilt} dw {dw} wE {wE}")
    
    # Constants for HH contributions
    tmp1 = (4*pi*pi)*(10**6)*Dcnst(rHH,"1H","1H")*Dcnst(rHH,"1H","1H")/4
    tmp2 = (1+(3*np.cos(2*tilt)*np.cos(2*tilt)))/48
    tmp3 = (3/4)*(np.sin(tilt)**4)
    tmp4 = (3/16)*(np.sin(2*tilt)**2)
    tmp5 = (3/4)*(7-(3*np.cos(2*tilt)))
    tmp6 = (3/2)*(5+(3*np.cos(2*tilt)))
    # Constants for HN contributions
    tmp7 = (4*pi*pi)*(10**6)*Dcnst(rHN,"1H","15N")*Dcnst(rHN,"1H","15N")/4
    tmp8 = np.sin(tilt)*np.sin(tilt)/3
    tmp9 = (3+np.cos(2*tilt))/4
    # Constants for CSA contributions
    tmp10 = (dw*dw)*(3/4)
    tmp11 = np.sin(tilt)*np.sin(tilt)/9
    tmp12 = (3+np.cos(2*tilt))/4
    
    
    
    
    # Constants for exchange contributions
    pG = 1-pE
    dwhex = dwhex*2*pi*10**3
    # Spectral density functions
    jw2wh = SDF(tau,S2,2*wH)
    jwwh = SDF(tau,S2,wH)
    jwwn = SDF(tau,S2,wN)
    jw2wr = SDF(tau,S2,2*wR)
    jwwr = SDF(tau,S2,wR)
    jw2wep2wr = SDF(tau,S2,2*(wE+wR))
    jwwep2wr = SDF(tau,S2,wE+(2*wR))
    jw2wepwr = SDF(tau,S2,(2*wE)+wR)
    jwwepwr = SDF(tau,S2,(wE+wR))
    jw2wem2wr = SDF(tau,S2,2*(wE-wR))
    jwwem2wr = SDF(tau,S2,wE-(2*wR))
    jw2wemwr = SDF(tau,S2,(2*wE)-wR)
    jwwemwr = SDF(tau,S2,(wE-wR))
    jwwhpwn = SDF(tau,S2,wH+wN)
    jwwhmwn = SDF(tau,S2,wH-wN)
    #calculating different contributions
    hhcontri =tmp1*( (tmp3*((2*jw2wepwr)+jw2wep2wr + jw2wem2wr + (2*jw2wemwr)))
                + (tmp5*jwwh) + (tmp6*jw2wh) +  (tmp4*((2*jwwepwr) + jwwep2wr + jwwem2wr + (2*jwwemwr))))
    hncontri = ((tmp7*tmp8*(jwwem2wr + (2*jwwemwr) + (2*jwwepwr) + jwwep2wr + (9*jwwn)))
          + (tmp7*tmp9*((3*jwwh) + (6*jwwhpwn) + jwwhmwn)))
    csacontri = ((tmp10*tmp11*(jwwem2wr + (2*jwwemwr) + (2*jwwepwr) + jwwep2wr))
          + (tmp10*tmp12*jwwh))
    r2contri = R20*np.sin(tilt)*np.sin(tilt)
    r1contri = R1*np.cos(tilt)*np.cos(tilt)
    excontri = (np.sin(tilt)*np.sin(tilt)*pG*pE*dwhex*dwhex*tau)/(1+(wE*wE*tau*tau))
    # print(f"r1contri {r1contri} r2contri {r2contri} hhcontri {hhcontri} hncontri {hncontri} csacontri {csacontri} excontri {excontri}")

    return r1contri+r2contri+hhcontri+hncontri+csacontri+excontri



# ===== Proton r1rho (NERRD) calculations =======================================
def ProtonR1rhoNERRD1(w1,tilt,B0,wR,rHH,dwcsa,tau,S2,R20):
    wH,wN = larmor_freq_calc(B0,"1H"),larmor_freq_calc(B0,"15N")
    R1,rHN,dwhex,pE = 0.0,1.02,0.0,0.0
    return ProtonR1rho(wH,wN,rHH,rHN,dwcsa,tau,S2,dwhex,pE,R1,R20,wR,w1,tilt)

# ===== Proton r1rho decay (NERRD) calculations =================================
def ProtonR1rhoNERRDDecay1(w1,tilt,t,B0,wR,rHH,dwcsa,tau,S2,R20,a):
    tmp1 = ProtonR1rhoNERRD1(w1,tilt,B0,wR,rHH,dwcsa,tau,S2,R20)
    return np.log(abs(a))-(t*tmp1)


# ===== Spectral density function =============================================
def SDF(tau,S2,w):
    return (2*tau*(1-S2))/(5*(1+(tau*tau*w*w)))


# ===== Dipolar constant calculation ==========================================
def Dcnst(r,nuc1,nuc2):
    """ r is in Angstroms and output coupling is given in Radians per second """
    if (nuc1 == "1H"):
        g1 = gammaH
    elif (nuc1 == "15N"):
        g1 = gammaN
    elif(nuc1 == "13C"):
        g1 = gammaC
    else:
        print("Invalid Nucleus")
        g1=None
    if (nuc2 == "1H"):
        g2 = gammaH
    elif (nuc2 == "15N"):
        g2 = gammaN
    elif(nuc2 == "13C"):
        g2 = gammaC
    else:
        print("Invalid Nucleus")
        g2=None
    return (g1*g2*mu0*hbar*10**(27))/(8*pi*pi*r**(3))

# ===== Larmor frequceny calculation ===========================================
def larmor_freq_calc(spectrometer_MHz,spin):
    """Returns Larmor frequency in kHz"""
    # larmor_freq_H=2*np.pi*spectrometer_MHz*10**3
    larmor_freq_H=spectrometer_MHz*10**3
    if (spin == "1H"):
        return larmor_freq_H
    if (spin == "15N"):
        larmor_freq_N = larmor_freq_H * (gammaN/gammaH)
        return larmor_freq_N
    if (spin == "13C"):
        larmor_freq_C= larmor_freq_H * (gammaC/gammaH)
        return larmor_freq_C
    else:
        print("Invalid Nucleus")
        return None
    

# ====== Dispersion plots for proton ==============
class HDispPlotwithFit:
    def __init__(self,decayfitout,dispfitout,fittype,wR,MC,pcutoff,prefilename,ymax1):
        # If you don not need dispersion fit
        s1,s2 = ToDisp(decayfitout),dispfitout
        ncol  = 5
        nrows = (s1.nres//ncol)+1
        nrf = s1.nrf
        xmax = 1.15*np.amax(s1.effrf)
        xmin = -0.1*xmax
        if (ymax1 == -1 ):
            ymax = 1.15*np.amax(s1.r1rhocor,0)
        elif (ymax1 == 0):
            ymax = np.amax(s1.r1rhocor)
        else:
            ymax = ymax1
        ymin = -0.1*ymax
        k=0
        filename2 = prefilename + '_' + fittype + 'nodispersion_withfit.pdf'
        filename3 = prefilename + '_' + fittype + 'dispersion_withfit.pdf'
        filename4 = prefilename + '_' + fittype + 'dispersion_withfit_ftest.pdf'
        #  fig2,fig3,fig4
        self.fig2, self.axes2 = plt.subplots(nrows,ncol,figsize=(ncol*2.3,1.8*nrows))
        self.fig3, self.axes3 = plt.subplots(nrows,ncol,figsize=(ncol*2.3,1.8*nrows))
        self.fig4, self.axes4 = plt.subplots(nrows,ncol,figsize=(ncol*2.3,1.8*nrows))
        for i in range(nrows):
            for j in range(ncol):
                if (k < s1.nres):
                    self.axes2[i,j].errorbar(np.array(s1.effrf)[:,k],np.array(s1.r1rhocor)[:,k],
                                yerr=np.array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes3[i,j].errorbar(np.array(s1.effrf)[:,k],np.array(s1.r1rhocor)[:,k],
                                yerr=np.array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes4[i,j].errorbar(np.array(s1.effrf)[:,k],np.array(s1.r1rhocor)[:,k],
                                yerr=np.array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    # Plot fitted curves
                    self.axes2[i,j].plot(s2.effrf[:],s2.r1rhosimref[k,:],color='red',linewidth=1.5)
                    self.axes3[i,j].plot(s2.effrf[:],s2.r1rhosim[k,:],color='green',linewidth=1.5)
                    #Write fitted parameters on the plots
                    tmptxt1 = r'R$_2^0$=' + str(round(s2.refmcfitvalues1[k][0],2)) \
                              + r' $\pm$ ' + str(round(s2.refmcfiterr1[k][0],2)) + r' s$^{-1}$'
                    tmptxt2 = r'R$_2^0$=' + str(round(s2.mcfitvalues1[k][2],2)) \
                              + r' $\pm$ ' + str(round(s2.mcfiterr1[k][2],2)) + r' s$^{-1}$'
                    tmptxt3 = r'$\tau_c$ =' + str(round((s2.mcfitvalues1[k][0])*1000,2)) \
                              + r' $\pm$ ' + str(round((s2.mcfiterr1[k][0])*1000,2)) + r' $\mu$s'
                    tmptxt5 = r'r$_{HH}^{eff}$=' + str(round(s2.rHH[k],2)) + r' $\AA$'
                    if ((fittype == 'NERRD') or (fittype == 'BMRDNERRD') ) :
                        xpos = 0.8*xmin
                        tmptxt4 = r'S$^2$=' + str(round(s2.mcfitvalues1[k][1],2)) \
                              + r' $\pm$ ' + str(round(s2.mcfiterr1[k][1],2))
                    if (fittype == 'BMRD'):
                        xpos = -xmin
                        tmptxt4 = r'$\omega_{ex}$=' + str(round(s2.mcfitvalues1[k][1],2)) \
                              + r' $\pm$ ' + str(round(s2.mcfiterr1[k][1],2))
                    if (fittype == 'BMRDNERRD'):
                        xpos = -xmin
                        tmptxt2 = r'R$_2^0$=' + str(round(s2.mcfitvalues1[k][3],2)) \
                                + r' $\pm$ ' + str(round(s2.mcfiterr1[k][3],2)) + r' s$^{-1}$'
                    self.axes2[i,j].text(xpos,0.85*ymax,tmptxt1,fontsize = 6)
                    self.axes3[i,j].text(xpos,0.85*ymax,tmptxt3,fontsize = 6)
                    self.axes3[i,j].text(xpos,0.70*ymax,tmptxt4,fontsize = 6)
                    self.axes3[i,j].text(xpos,0.55*ymax,tmptxt2,fontsize = 6)
                    if not (fittype == 'BMRD'):
                        self.axes3[i,j].text(xpos,0.40*ymax,tmptxt5,fontsize = 6)
                    if (s2.ftestresults[k,0]==1):
                        self.axes4[i,j].plot(s2.effrf[:],s2.r1rhosim[k,:],color='green',linewidth=1.5)
                        self.axes4[i,j].text(xpos,0.85*ymax,tmptxt3,fontsize = 6)
                        self.axes4[i,j].text(xpos,0.70*ymax,tmptxt4,fontsize = 6)
                        self.axes4[i,j].text(xpos,0.55*ymax,tmptxt2,fontsize = 6)
                        if not (fittype == 'BMRD'):
                            self.axes4[i,j].text(xpos,0.40*ymax,tmptxt5,fontsize = 6)
                    else :
                        self.axes4[i,j].plot(s2.effrf[:],s2.r1rhosimref[k,:],color='red',linewidth=1.5)
                        self.axes4[i,j].text(xpos,0.85*ymax,tmptxt1,fontsize = 6)
                    self.axes2[i,j].set_xlim([xmin,xmax])
                    self.axes3[i,j].set_xlim([xmin,xmax])
                    self.axes4[i,j].set_xlim([xmin,xmax])
                    self.axes2[i,j].set_ylim([ymin,ymax])
                    self.axes3[i,j].set_ylim([ymin,ymax])
                    self.axes4[i,j].set_ylim([ymin,ymax])
                    self.axes2[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes3[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes4[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes2[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    self.axes3[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    self.axes4[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    DecaySubPlot(self.axes2[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    DecaySubPlot(self.axes3[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    DecaySubPlot(self.axes4[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    k=k+1
                else:
                    self.axes2[i,j].set_axis_off()
                    self.axes3[i,j].set_axis_off()
                    self.axes4[i,j].set_axis_off()
                    k=k+1
        # adjust space between the subplots
        self.fig2.subplots_adjust(wspace=0.4, hspace=0.9)
        self.fig3.subplots_adjust(wspace=0.4, hspace=0.9)
        self.fig4.subplots_adjust(wspace=0.4, hspace=0.9)
        # save plots to pdf
        self.fig2.savefig(filename2,dpi=1200)
        self.fig3.savefig(filename3,dpi=1200)
        self.fig4.savefig(filename4,dpi=1200)