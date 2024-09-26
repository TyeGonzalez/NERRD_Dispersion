import time
from numpy import array, amax, log, divide, random, mean, std, arange
import lmfit
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt

#Fitting Class for all residues
class DecayFit:
    def __init__(self, rfout, MC, ncpucore):
        start_time = time.time()
        
        self.normvolumes = rfout.normvolumes
        self.normnoises = rfout.normnoises
        self.effrfkhz = rfout.effrfkhz
        self.timelist = rfout.timelist
        self.seqname = rfout.seqname
        self.seqnum = rfout.seqnum
        self.nres = rfout.nres
        self.offsets = rfout.offsets
        self.angles = rfout.angles

        fitvalues1 = []
        fiterr1 = []
        mcfitvalues1 = []
        mcfiterr1 = []

        # Initialize multiprocessing pool
        pool_start_time = time.time()
        with multiprocessing.Pool(ncpucore) as pool:
            pool_initialization_time = time.time() - pool_start_time
            
            processes_start_time = time.time()
            processes = [pool.apply_async(DecayFitsr, args=(rfout, i, MC)) for i in range(rfout.nres)]
            processes_creation_time = time.time() - processes_start_time

            for p in processes:
                s2 = p.get()
                fitvalues1.append(s2.popt)
                
                # Handle fit errors
                if not (s2.perr[0] or s2.perr[1]):
                    fiterr1.append([0.0, 0.0])
                else:
                    fiterr1.append(s2.perr)
                
                mcfitvalues1.append(s2.mcfitvalues1)
                mcfiterr1.append(s2.mcfiterr1)

        data_conversion_start_time = time.time()
        # Convert to numpy arrays
        self.fitvalues = array(fitvalues1)
        self.fitvalues[:, 1] *= 1000
        self.fiterr = array(fiterr1)
        self.fiterr[:, 1] *= 1000
        self.mcfitvalues = array(mcfitvalues1)
        self.mcfitvalues[:, 1] *= 1000
        self.mcfiterr = array(mcfiterr1)
        self.mcfiterr[:, 1] *= 1000
        data_conversion_time = time.time() - data_conversion_start_time

        total_time = time.time() - start_time
        
        # Print timing information
        print(f"Total execution time: {total_time:.2f} seconds")
        print(f"Pool initialization time: {pool_initialization_time:.2f} seconds")
        print(f"Processes creation time: {processes_creation_time:.2f} seconds")
        print(f"Data conversion time: {data_conversion_time:.2f} seconds")
    def __str__(self):
        return f"DecayFit object with eff. field: {self.effrfkhz[0]:.2f}, timelist: {self.timelist} \n "
    def __repr__(self):
        return self.__str__()
    

# Fit the decay curves with exponential function with monte-carlo error analysis
#Fitting function for single residue
class DecayFitsr:
    def __init__(self,rfout,i,MC):
        s=rfout
        x = s.timelist
        xmax = 1.2*amax(s.timelist)
        xt = arange(0,xmax)
        y = abs(s.normvolumes[i])
        logy = log(y)
        yerr = abs(s.normnoises[i])
        logyerr = abs(divide(yerr,y))
        #print(logyerr)
        #:::::::::::::::::::::::::::::
        # Fitting with "curve_fit function". This changed to lmfit. in the following section
        #self.popt,self.pcov = curve_fit(MonoExpLin, x, logy, p0 = [1.0,10.0],sigma=logyerr,
        #                       absolute_sigma=True,bounds=([0,0],[1000,1000]))
        #self.perr = sqrt(diag(self.pcov))
        #:::::::::::::::::::::::::::::
        fitmodel =lmfit.Model(MonoExpLin,independent_vars=['x'])
        fitmodel.set_param_hint('a0',value=1.0,min=0,max=1000)
        fitmodel.set_param_hint('R0',value=1.0,min=0,max=1000)
        params=fitmodel.make_params()
        fitresult = fitmodel.fit(logy,params,x=x,weights=1./logyerr)
        self.popt = [fitresult.params['a0'].value,fitresult.params['R0'].value]
        self.perr = [fitresult.params['a0'].stderr,fitresult.params['R0'].stderr]
        #print(self.popt[0],self.popt[1])
        tmpfit = []
        for j in range(MC):
            #print(popt[0],popt[1])
            ysim = MonoExpLin(x,self.popt[0],self.popt[1])+random.normal(0,logyerr)
            #print(ysim)
            try:
                #::::::::::::::::::::::::::::::::
                # Fitting with "curve_fit function". This changed to lmfit.
                #poptsim, pcovsim = curve_fit(MonoExpLin,x,ysim, p0=[1.0,10.0])
                #::::::::::::::::::::::::::::::::
                fitresultsim = fitmodel.fit(ysim,params,x=x,weights=1./logyerr)
                poptsim = [fitresultsim.params['a0'].value,fitresultsim.params['R0'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues1=mean(tmpfit,0)
        self.mcfiterr1=std(tmpfit,0)
    
    
# ===== Mono-exponential function (in logscale)definition =====================
def MonoExpLin(x, a0, R0,):
    return  log(a0)-(x*R0)

# Mono-exponential function definition
def MonoExp(x, a0, R0,):
    return  a0*np.exp(-x*R0)

def DecaySubPlot(ax,xlabel,ylabel,title):
    ax.set_xlabel(xlabel,fontsize=8)
    ax.set_ylabel(ylabel,fontsize=8)
    ax.set_title(title,fontsize=10)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.minorticks_on()
    return

# Default plot parameters
textsize = 8
markersize = 4
capsize = markersize/2
linesize = 2


class DecayPlots:
    def __init__(self, rfout, MC, fitplotflag, ncpucore):
        self.s = rfout
        self.ncol = 5
        self.nrows = (self.s.nres // self.ncol) + 1

        self.x = self.s.timelist
        self.xmax = 1.15 * np.amax(self.s.timelist)
        self.xt = np.arange(0, self.xmax)
        self.rfout = rfout
        self.MC = MC
        self.fitplotflag = fitplotflag
        self.ncpucore = ncpucore
        
    def plot_decays(self):

        self.fig, self.axes = plt.subplots(self.nrows, self.ncol, figsize=(self.ncol * 2.3, 1.8 * self.nrows))

        if self.fitplotflag == 'yes':
            fitout = DecayFit(self.rfout, self.MC, self.ncpucore)
            self.fitvalues = fitout.fitvalues
            self.fiterr = fitout.fiterr
            self.mcfitvalues = fitout.mcfitvalues
            self.mcfiterr = fitout.mcfiterr

        max_normvolume = np.amax(self.s.normvolumes[0])
        min_normvolume = np.amin(self.s.normvolumes[0])

        for k in range(self.s.nres):
            i, j = divmod(k,self. ncol)
            self.axes[i, j].errorbar(self.s.timelist, self.s.normvolumes[k], yerr=abs(self.s.normnoises[k]), marker='o',
                                     color='green', linestyle='', markersize=markersize, ecolor='black',
                                     capsize=capsize, lw=capsize)
            if self.fitplotflag == 'yes':
                a0 = fitout.mcfitvalues[k, 0]
                R0 = fitout.mcfitvalues[k, 1]
                a0err = fitout.mcfiterr[k, 0]
                R0err = fitout.mcfiterr[k, 1]
                self.axes[i, j].plot(self.xt, MonoExp(self.xt, a0, R0 / 1000), 'r-')
                tmptxt = f"R1\u03C1={round(R0, 2)}\u00B1{round(R0err, 2)}"
                self.axes[i, j].text(0.2 * self.xmax, 1.3, tmptxt, fontsize=7)

            self.axes[i, j].set_xlim([-10.0, 1.2 * self.xmax])
            self.axes[i, j].set_ylim([0, 1.5 * max_normvolume])
            DecaySubPlot(self.axes[i, j], 'Delay(in ms.)', 'rel.int.', self.s.seqname[k])

        # Turn off any remaining empty subplots
        for k in range(self.s.nres, self.nrows * self.ncol):
            i, j = divmod(k, self.ncol)
            self.axes[i, j].set_axis_off()

        self.fig.subplots_adjust(wspace=0.4, hspace=0.9)
        self.fig.savefig(self.s.prefilename + 'decay.pdf')
        plt.close()