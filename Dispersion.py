# load all necessary pacakages
from defaultpackages import *
from Globalparam import *


def MyPlot(axs,x,y,xlim = None, ylim = None,
           xlabel="",ylabel="",title="",legend=None,xscale = "linear",yscale="linear",
           marker='o', color='darkmagenta',
           markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = ''):
    x = array(x)
    y = array(y)
    #fig,axs = plt.subplots()
    if (len(y.shape) == 2):
        nplots=y.shape[0]
        for i in range(nplots):
            axs.plot(x,y[i])
    else:
         axs.plot(x,y,marker=marker, color=color,
         markeredgecolor = markeredgecolor,markeredgewidth=0.5,
         linestyle = linestyle)
    axs.set_xscale(xscale)
    axs.set_yscale(yscale)
    if not (legend == None):
        axs.legend(legend)
    axs.set_xlabel(xlabel,fontsize=16)
    axs.set_ylabel(ylabel,fontsize=16)
    axs.set_title(title,fontsize=16)
    if not (xlim == None):
        axs.set_xlim(xlim)
    if not (ylim == None):
        axs.set_ylim(ylim)
    axs.tick_params(axis='both', which='major', labelsize=14)
    return axs

def MyContourPlot(axs,x,y,z,nlevels,xlabel,ylabel,zlabel,xscale='linear',yscale='linear'):
    fig = plt.gcf()
    zmax = round(1.1*z.max())
    zmin = round(1.1*z.min())
    incr = round((zmax-zmin)/nlevels)
    if (incr==0):
        zmax = 1.1*z.max()
        zmin = 1.1*z.min()
        incr = (zmax-zmin)/nlevels
    levels = arange(zmin, zmax, incr)
    norm = cm.colors.Normalize(vmax=zmax, vmin=zmin)
    #cmap = cm.PRGn
    cmap = cm.turbo
    #fig, axs = plt.subplots()
    cset1 = axs.contourf(x, y, z, levels, norm=norm,cmap=cm.get_cmap(cmap, len(levels) - 1))
    cset2 = axs.contour(x, y, z, cset1.levels, colors='k',linewidths=1)
    #cset3 = axs.contour(x, y, z, colors='g', linewidths=2)
    axs.set_xscale(xscale)
    axs.set_yscale(yscale)
    axs.set_title(zlabel,fontsize=16)
    axs.set_xlabel(xlabel,fontsize=16)
    axs.set_ylabel(ylabel,fontsize=16)
    axs.tick_params(axis='both', which='major', labelsize=14)
    fig.colorbar(cset1,ax=axs,ticks=levels)
    return axs


# Import peaklist file into a data frame
# '\s+' for white space and '\t' for tabs
def ImportFile(filename, filetype):
    #X = pd.read_csv(filename, sep="\t", header=0)
    if (filetype == "ccpnmr"):
        skiprows = 1
        sep1 = '\s+ '
    elif ((filetype == "xeasy") or ((filetype == "sparky"))):
        skiprows = 3
        sep1 = '\s+ '
    else :
        skiprows = 0
        sep1 = '\s+ '
    #tmp = read_table(filename, sep=sep1, header=header1,index_col=False,engine='python')
    tmp = read_table(filename,delim_whitespace=True,skiprows=skiprows,header=None)
    X = DataFrame(tmp)
    return X

#===== Begin: Definition of functions for CCPNMR loading =======================
#ccpnmrprojdir = '/Users/Suresh/Documents/home_workspace/SMILE/SH3_NUSrelax.ccpn'
#filename = os.listdir(ccpnmrprojdir+'/ccpnv3/ccp/nmr/Nmr/')[0]
class loadccpnmrproject:
    def __init__(self,ccpnmrprojdir):
        peaklistnames,spectranames,heights,volumes,seqname,residues = [],[],[],[],[],[]
        si1,dw1,ref1,refpos1,sf1 = {},{},{},{},{}
        si2,dw2,ref2,refpos2,sf2 = {},{},{},{},{}
        si3,dw3,ref3,refpos3,sf3 = {},{},{},{},{}
        atomcodes1,atomcodes2,ndim,scale = {}, {}, {},{}
        assign1,assign2,assign3,pos1,pos2,pos3 = [],[],[],[],[],[]
        filename = os.listdir(ccpnmrprojdir+'/ccpnv3/ccp/nmr/Nmr/')[0]
        tree = ET.parse(ccpnmrprojdir+'/ccpnv3/ccp/nmr/Nmr/'+filename)
        root = tree.getroot()
        rootname1 = "NMR.NmrProject/NMR.NmrProject.experiments/NMR.Experiment"
        rootname2 = rootname1 + "/NMR.Experiment.dataSources/NMR.DataSource"
        # find sf1 sf2 and sf3 for all spectra
        for elem1 in root.findall(rootname1):
            try:
                spectraname = (elem1.find
                               ("NMR.Experiment.dataSources/NMR.DataSource/NMR.DataSource.name/IMPL.Line")).text
                ndim[spectraname] = int((elem1.find
                               ("NMR.Experiment.dataSources/NMR.DataSource")).attrib["numDim"])
                #print(spectraname,ndim)
                for elem2 in elem1.findall('NMR.Experiment.expDims/NMR.ExpDim'):
                    if (elem2.attrib["dim"] == '1'):
                       sf1[spectraname]=float((elem2.find("NMR.ExpDim.expDimRefs/NMR.ExpDimRef")).attrib["sf"])
                    if (elem2.attrib["dim"] == '2'):
                       sf2[spectraname]=float((elem2.find("NMR.ExpDim.expDimRefs/NMR.ExpDimRef")).attrib["sf"])
                    if (elem2.attrib["dim"] == '3'):
                       sf3[spectraname]=float((elem2.find("NMR.ExpDim.expDimRefs/NMR.ExpDimRef")).attrib["sf"])
            except AttributeError:
                pass
        # Find residue names and residue atoms
        tmp12 = 999
        print(filename, tree, root)
        for elem1 in root[0].find('NMR.NmrProject.resonanceGroups'):
            try:
                tmp11 = elem1.attrib['seqCode']+elem1.attrib['residueType']
            except KeyError:
                tmp11 = str(tmp12)+'NAN'
                tmp12 = tmp12 -1
            residues.append(tmp11)
            try:
                tmp13 = elem1.find('NMR.ResonanceGroup.resonances').text.split()
            except AttributeError:
                tmp13 = ['-999']
            tmp14 = dict.fromkeys(tmp13,tmp11)
            atomcodes1.update(tmp14)
        # find atom types
        try:
            for elem1 in root[0].find('NMR.NmrProject.resonances'):
                tmp11 = elem1.attrib['_ID']
                try:
                    tmp12=atomcodes1[tmp11]+'-'+elem1.attrib['implName']
                    tmp13 = elem1.find('NMR.Resonance.peakDimContribs').text.split()
                except KeyError:
                    tmp12=atomcodes1[tmp11]+'-'+elem1.attrib['isotopeCode'][-1]+'?'
                except AttributeError:
                    tmp13 = ['-999']
                tmp14 = dict.fromkeys(tmp13,tmp12)
                atomcodes2.update(tmp14)
        except TypeError:
            atomcodes2.update({'-999':'NAN'})
        #print(atomcodes2)
        #Find spectra and peaks
        for elem1 in root.findall(rootname2):
            # Here each elem1 represents each spectrum in ccpnmr project
            spectraname = (elem1.find("NMR.DataSource.name/IMPL.Line")).text
            spectranames.append(spectraname)
            try:
                scale[spectraname] = float(elem1.attrib["scale"])
            except KeyError:
                scale[spectraname] = 1.0
            #print(spectraname)
            for elem2 in elem1.findall("NMR.DataSource.dataDims/NMR.FreqDataDim"):
                tmp31 = elem2.find("NMR.FreqDataDim.dataDimRefs/NMR.DataDimRef").attrib["refValue"]
                tmp32 = elem2.find("NMR.FreqDataDim.dataDimRefs/NMR.DataDimRef").attrib["refPoint"]
                if (elem2.attrib["dim"]=="1"):
                    dw1[spectraname] = float(elem2.attrib["valuePerPoint"])
                    si1[spectraname] = float(elem2.attrib["numPoints"])
                    ref1[spectraname] = float(tmp31)
                    refpos1[spectraname] = float(tmp32)
                if (elem2.attrib["dim"]=="2"):
                    dw2[spectraname] = float(elem2.attrib["valuePerPoint"])
                    si2[spectraname] = float(elem2.attrib["numPoints"])
                    ref2[spectraname] = float(tmp31)
                    refpos2[spectraname] = float(tmp32)
                if (elem2.attrib["dim"]=="3"):
                    dw3[spectraname] = float(elem2.attrib["valuePerPoint"])
                    si3[spectraname] = float(elem2.attrib["numPoints"])
                    ref3[spectraname] = float(tmp31)
                    refpos3[spectraname] = float(tmp32)
            for elem2 in elem1.findall("NMR.DataSource.peakLists/NMR.PeakList"):
                peaklistnum = elem2.attrib["serial"]
                tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9 = [],[],[],[],[],[],[],[],[]
                for elem3 in elem2.iter("NMR.Peak"):
                    # Load peak heights
                    tmp1.append(float(elem3.attrib["height"]))
                    # Load peak volumes
                    try:
                        tmp2.append(float(elem3.attrib["volume"]))
                    except KeyError:
                        tmp2.append(-999.999)
                    #Load peak shifts and assignements
                    for elem4 in elem3.findall("NMR.Peak.peakDims/NMR.PeakDim"):
                        if (elem4.attrib["dim"] == '1'):
                            tmp4.append(float(elem4.attrib["position"]))
                            # for non-assigning peaks put -999
                            try:
                                tmp7.append(atomcodes2[elem4.find(
                                    'NMR.PeakDim.peakDimContribs/NMR.PeakDimContrib').attrib['_ID']])
                            except AttributeError:
                                # tmp7.append(atomcodes2['-999'])
                                tmp7.append(atomcodes2.setdefault('-999', '-999'))
                        if (elem4.attrib["dim"] == '2'):
                            tmp5.append(float(elem4.attrib["position"]))
                            # for non-assigning peaks put -999
                            try:
                                tmp8.append(atomcodes2[elem4.find(
                                    'NMR.PeakDim.peakDimContribs/NMR.PeakDimContrib').attrib['_ID']])
                            except AttributeError:
                                # tmp8.append(atomcodes2['-999'])
                                tmp8.append(atomcodes2.setdefault('-999', '-999'))
                        if (elem4.attrib["dim"] == '3'):
                            tmp6.append(float(elem4.attrib["position"]))
                            # for non-assigning peaks put -999
                            try:
                                tmp19 = elem4.findall('NMR.PeakDim.peakDimContribs/NMR.PeakDimContrib')
                                tmp29 = ' '
                                for i in range(len(tmp19)):
                                    tmp29 = tmp29 + atomcodes2[tmp19[i].attrib['_ID']]+ ' '
                                tmp9.append(tmp29)
                                #tmp9.append(atomcodes2[elem4.find(
                                #    'NMR.PeakDim.peakDimContribs/NMR.PeakDimContrib').attrib['_ID']])
                            except AttributeError:
                                tmp9.append(atomcodes2['-999'])
                peaklistnames.append(spectraname+"."+peaklistnum)
                #print(peaklistnames)
                heights.append(array(tmp1))
                # print(len(volumes))
                # print(len(tmp2))
                volumes.append(array(tmp2))
                assign1.append(tmp7)
                assign2.append(tmp8)
                assign3.append(tmp9)
                pos1.append(array(tmp4))
                pos2.append(array(tmp5))
                pos3.append(array(tmp6))
        self.peaklistnames,self.spectranames = peaklistnames,spectranames
        self.heights = heights
        self.volumes = volumes
        self.seqname = seqname
        self.ndim,self.scale = ndim,scale
        self.ref1,self.dw1,self.si1,self.pos1,self.sf1,self.assign1=ref1,dw1,si1,pos1,sf1,assign1
        self.ref2,self.dw2,self.si2,self.pos2,self.sf2,self.assign2=ref2,dw2,si2,pos2,sf2,assign2
        self.ref3,self.dw3,self.si3,self.pos3,self.sf3,self.assign3=ref3,dw3,si3,pos3,sf3,assign3
        self.refpos1,self.refpos2,self.refpos3 = refpos1,refpos2,refpos3
        self.nspectra = len(peaklistnames)
#ccpnmrdata = loadccpnmrproject(ccpnmrprojdir)

class loadccpnmrpeaks:
    def __init__(self,ccpnmrdata,filename,timevalue,noisevalue,nscans,centerppm,spin):
        self.seqname= []
        self.seqnum = []
        self.noises = []
        self.heights= []
        self.offsets= []
        self.timevalue=timevalue
        try:
            tmpindex = [ccpnmrdata.peaklistnames.index(i)
                        for i in ccpnmrdata.peaklistnames if filename in i][0]
            print(filename)
            self.spectraname = filename[:-2]
            self.ndim = ccpnmrdata.ndim[self.spectraname]
            if (spin=='1H'):
                tmpseqname = array(ccpnmrdata.assign1[tmpindex])
            else:
                tmpseqname = array(ccpnmrdata.assign2[tmpindex])
            try:
                tmpseqnum =array([int(findall(r'\d+',i)[0]) for i in tmpseqname])
            except IndexError:
                tmpseqnum =array([i for i in range(len(tmpseqname))])
            tmpsort = tmpseqnum.argsort()
            self.seqname = tmpseqname[tmpsort]
            self.nres = len(self.seqname)
            #self.seqnum=[int(re.findall(r'\d+',self.seqname[i])[0]) for i in range(self.nres)]
            self.seqnum = tmpseqnum[tmpsort]
            self.volumes =(ccpnmrdata.scale[self.spectraname])*(1/nscans)*(ccpnmrdata.volumes[tmpindex][tmpsort])
            self.heights =(ccpnmrdata.scale[self.spectraname])*(1/nscans)*(ccpnmrdata.heights[tmpindex][tmpsort])
            self.noises = [noisevalue/nscans for i in range(self.nres)]
            dw1 =ccpnmrdata.dw1[self.spectraname]/ccpnmrdata.sf1[self.spectraname]
            dw2 =ccpnmrdata.dw2[self.spectraname]/ccpnmrdata.sf2[self.spectraname]
            if (ccpnmrdata.refpos1[self.spectraname]==1):
                highref1 = ccpnmrdata.ref1[self.spectraname]
            else:
                highref1 = (ccpnmrdata.si1[self.spectraname]*dw1/2) + ccpnmrdata.ref1[self.spectraname]
            if (ccpnmrdata.refpos2[self.spectraname]==1):
                highref2 = ccpnmrdata.ref2[self.spectraname]
            else:
                highref2 = (ccpnmrdata.si2[self.spectraname]*dw2/2) + ccpnmrdata.ref2[self.spectraname]
            self.freq1 = highref1 - (ccpnmrdata.pos1[tmpindex][tmpsort]*dw1)
            self.freq2 = highref2 - (ccpnmrdata.pos2[tmpindex][tmpsort]*dw2)
            self.freq1 = highref1 - (ccpnmrdata.pos1[tmpindex][tmpsort]*dw1)
            self.freq2 = highref2 - (ccpnmrdata.pos2[tmpindex][tmpsort]*dw2)
            self.assign1 = array(ccpnmrdata.assign1[tmpindex])[tmpsort]
            self.assign2 = array(ccpnmrdata.assign2[tmpindex])[tmpsort]
            if (self.ndim == 3):
                dw3 =ccpnmrdata.dw3[self.spectraname]/ccpnmrdata.sf3[self.spectraname]
                if (ccpnmrdata.refpos3[self.spectraname]==1):
                    highref3 = ccpnmrdata.ref3[self.spectraname]
                else:
                    highref3 = (ccpnmrdata.si3[self.spectraname]*dw3/2) + ccpnmrdata.ref3[self.spectraname]
                self.freq3 = highref3 - (ccpnmrdata.pos3[tmpindex][tmpsort]*dw3)
                self.assign3 = array(ccpnmrdata.assign3[tmpindex])[tmpsort]
            if (spin == '1H'):
                self.offsets = 0.001*ccpnmrdata.sf1[self.spectraname]*(self.freq1 - centerppm)
            else:
                self.offsets = 0.001*ccpnmrdata.sf2[self.spectraname]*(self.freq2 - centerppm)
        except IndexError:
            print("Peaklist = '' %s ''does not exist in the ccpnmr project: Please Check"%filename)

#tmp = loadccpnmrpeaks(ccpnmrdata,filenamelist1[0],timelist1[0],noisevolumelist1[0],nscanlist1[0],centerppm,spin)

# ===== Class definition for multiple spectra peaks into a single object =======
class loadccpnmrpeakslist:
    def __init__(self,ccpnmrdata,filenamelist,rfkHz,offset,timelist,noisevaluelist,
                 nscanlist,centerppm,prefilename,spin):
        nfiles = len(filenamelist)
        ntimes = len(timelist)
        self.rfkHz = rfkHz
        self.offset = offset
        self.prefilename = prefilename + '_' +str(rfkHz) + 'kHz_' + str(offset) + 'kHz_'
        self.volumes = []
        self.noises  = []
        try:
            for i in range(nfiles):
                # print(i)
                # print(nfiles)
                tmp1 = loadccpnmrpeaks(ccpnmrdata,filenamelist[i],timelist[i],
                                       noisevaluelist[i],nscanlist[i],centerppm,spin)
                self.volumes.append(tmp1.heights)
                self.noises.append(tmp1.noises)
            for sub_array in self.volumes:
                print(len(sub_array))
            self.assign1 = tmp1.assign1
            self.assign2 = tmp1.assign2
            self.freq1 = tmp1.freq1
            self.freq2 = tmp1.freq2
            self.ndim = tmp1.ndim
            if (self.ndim == 3):
                self.assign3 = tmp1.assign3
                self.freq3 = tmp1.freq3
            self.seqname= tmp1.seqname
            self.seqnum = tmp1.seqnum
            self.offsets = tmp1.offsets + offset
            self.timelist = array(timelist)
            self.volumes = array(self.volumes)
            self.noises  = array(self.noises)
            self.nres    = len(self.seqnum)
            self.ntimes  = ntimes
            self.angles  = arctan(divide(rfkHz, self.offsets))
            self.normvolumes = transpose(divide(self.volumes,self.volumes[0]))
            self.normnoises  = transpose(divide(self.noises,self.volumes[0]))
            self.effrfkhz = sqrt(((rfkHz)**2) + ((self.offsets)**2))
        except IndexError:
            if not (nfiles == ntimes):
                print("Total no.of peaklist doesnot match with total no. of decay times: Please Check")
            if not (nfiles == len(noisevaluelist)):
                print("Total no.of peaklist doesnot match with total no. of noise values: Please Check")
            if not (nfiles == len(nscanlist)):
                print("Total no.of peaklist doesnot match with total no. of nscan values: Please Check")

#tmp1=loadccpnmrpeakslist(ccpnmrdata,filenamelist1,spinlockrf1, offset1,timelist1,
#                          noisevolumelist1,nscanlist1,centerppm,prefilename,spin)

#===== End: Definition of functions for CCPNMR loading =========================

# Class definition for single spectrum peaks into a single object
class loadpeaks:
    def __init__(self,filename,timevalue,noisevalue,nscans,centerppm,filetype,
                 namecolumn,volumecolumn,freqcolumn,spin,B0):
        self.seqname= []
        self.seqnum = []
        self.noises = []
        self.offsets= []
        self.timevalue=timevalue
        tmp2 = ImportFile(filename,filetype)
        self.nres = len(tmp2)
        tmp3=tmp2.iloc[0:self.nres,namecolumn-1].to_numpy()
        tmp4 = []
        for i in range(self.nres):
            tmp = findall(r'\d+',str(tmp3[i]))[0]
            tmp4.append(int(tmp))
        tmp5 = concat([tmp2,DataFrame({'seqnum':tmp4})],axis=1)
        tmp6 = tmp5.sort_values(by=[tmp5.columns[len(tmp2.columns)]])
        tmp3 = tmp6.iloc[0:self.nres,namecolumn-1].to_numpy()
        for i in range(self.nres):
            tmp = str(tmp3[i]).strip()
            self.seqname.append(tmp)
            self.seqnum.append(tmp6.iloc[i])
            self.noises.append(noisevalue)
        self.volumes = (1/nscans)*tmp6.iloc[0:self.nres,volumecolumn-1].to_numpy()
        if (spin == '1H'):
            self.offsets = 0.001*B0*((tmp6.iloc[0:self.nres,freqcolumn-1].to_numpy()) - centerppm)
        if (spin == '15N'):
            self.offsets = 0.001*(gammaN/gammaH)*B0*((tmp6.iloc[0:self.nres,freqcolumn-1].to_numpy()) - centerppm)
        if (spin == '13C'):
            self.offsets = 0.001*(gammaC/gammaH)*B0*((tmp6.iloc[0:self.nres,freqcolumn-1].to_numpy()) - centerppm)
        self.noises  = divide(self.noises,nscans)

# Class definition for multiple spectra peaks into a single object
class loadpeakslist:
    def __init__(self,filenamelist,rfkHz,offset,timelist,noisevaluelist,nscanlist,centerppm,filetype,namecolumn,
                 volumecolumn,freqcolumn,prefilename,spin,B0):
        nfiles = len(filenamelist)
        ntimes = len(timelist)
        self.rfkHz = rfkHz
        self.offset = offset
        self.prefilename = prefilename + '_' +str(rfkHz) + 'kHz_' + str(offset) + 'kHz_'
        self.volumes = []
        self.noises  = []
        if not (nfiles == ntimes):
            print('files list and time list are not equal in length')
            exit()
        for i in range(nfiles):
            #print(i)
            tmp1 = loadpeaks(filenamelist[i],timelist[i],noisevaluelist[i],nscanlist[i],centerppm,filetype,
                             namecolumn,volumecolumn,freqcolumn,spin,B0)
            self.volumes.append(tmp1.volumes)
            self.noises.append(tmp1.noises)
        self.seqname= tmp1.seqname
        self.seqnum = tmp1.seqnum
        self.offsets = tmp1.offsets + offset
        self.timelist = array(timelist)
        self.volumes = array(self.volumes)
        self.noises  = array(self.noises)
        self.nres    = len(self.seqnum)
        self.ntimes  = ntimes
        self.angles  = arctan(divide(rfkHz, self.offsets))
        try:
            self.normvolumes = transpose(divide(self.volumes,self.volumes[0]))
            self.normnoises  = transpose(divide(self.noises,self.volumes[0]))
            self.effrfkhz = sqrt(((rfkHz)**2) + ((self.offsets)**2))
        except ValueError:
            print("Number of peaks in the %s file does not match with other files.Please fix this"
                  %filenamelist[0])
            raise

# Fit the decay curves with exponential function with monte-carlo error analysis
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


import time

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

# ===== Decay fit results to dispersion fitting module  =======================
class ToDisp:
    def __init__(self, rffitlist):
        self.nrf = len(rffitlist)
        self.nres = rffitlist[0].nres
        self.seqname = rffitlist[0].seqname
        effrf,angles,times,normvolumes,normnoises = [],[],[],[],[]
        ntimes,r1rho,r1rhoerr = [],[],[]
        for i in range(self.nrf):
            ntimes.append(len(rffitlist[i].timelist))
            angles.append(rffitlist[i].angles)
            effrf.append(rffitlist[i].effrfkhz)
            times.append(rffitlist[i].timelist)
            normvolumes.append(rffitlist[i].normvolumes)
            normnoises.append(rffitlist[i].normnoises)
            r1rho.append(rffitlist[i].mcfitvalues[:,1])
            r1rhoerr.append(rffitlist[i].mcfiterr[:,1])
        normvolumes1 = zeros((self.nres,sum(ntimes)))
        normnoises1 = zeros((self.nres,sum(ntimes)))
        rfangletime = zeros((self.nres,sum(ntimes),3))
        for k in range(self.nres):
            tmp1,tmp2,tmp3 = [],[],[]
            for i in range(self.nrf):
                for j in range(ntimes[i]):
                    tmp1.append([effrf[i][k],angles[i][j],times[i][j]*0.001])
                    tmp2.append(normvolumes[i][k][j])
                    tmp3.append(normnoises[i][k][j])
            rfangletime[k,:,:] = array(tmp1)
            normvolumes1[k,:] = array(tmp2)
            normnoises1[k,:] = array(tmp3)
        self.times = [[x*0.001 for x in y] for y in times]
        self.effrf = effrf
        self.angles = angles
        self.normvolumes = normvolumes
        self.normnoises = normnoises
        self.ntimes = ntimes
        self.r1rho = r1rho
        self.r1rhoerr = r1rhoerr
        self.r1rhocor = self.r1rho/(sin(self.angles)*sin(self.angles))
        self.r1rhoerrcor = self.r1rhoerr/(sin(self.angles)*sin(self.angles))
        self.rfangletime = rfangletime
        self.normvolumes1 = normvolumes1
        self.normnoises1 = normnoises1

# ===== Class for No Dispersion Decay Fit ============================
class NoDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,r200,r20min,r20max,a0,a0min,a0max,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        fitmodel = lmfit.Model(NoDispDecay1,independent_vars=['w1','tilt','t'])
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
            tmpysimref = (NoDispDecay1(x1,x2,x3,self.poptref[0],self.poptref[1]) +
                          random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(logy, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsimref =[fitresult.params['R20'].value,fitresult.params['a'].value]
                tmpfit.append(poptsimref)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.refmcfitvalues2=mean(tmpfit,0)
        self.refmcfiterr2=std(tmpfit,0)
        self.tmpysimref = NoDispDecay1(x1,x2,x3,mean(tmpfit,0)[0],mean(tmpfit,0)[1])
        self.tmpr1rhosimref = NoDisp1(x4,x6,mean(tmpfit,0)[0])

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
        xmax = 1.15*amax(s1.effrf)
        xmin = -0.1*xmax
        if (ymax1 == -1 ):
            ymax = 1.15*amax(s1.r1rhocor,0)
        elif (ymax1 == 0):
            ymax = amax(s1.r1rhocor)
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
                    self.axes1[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                    yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                    markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                    markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes1[i,j].set_xlim([-5.0,xmax])
                    self.axes1[i,j].set_ylim([ymin,ymax])
                    self.axes1[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    if (spin == '1H'):
                        self.axes1[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    MySubPlot(self.axes1[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                              s1.seqname[k])
                    k=k+1
                else:
                    self.axes1[i,j].set_axis_off()
                    k=k+1
        self.fig1.subplots_adjust(wspace=0.4, hspace=0.9)
        # save plots to pdf
        self.fig1.savefig(filename1,dpi=1200)
#DispPlotwithoutFit(spin,decayfitout,dispfitflag,fittype,wR,MC,pcutoff,prefilename,220)

#===============================================================================
# ===== Proton fitting functions ===============================================
#===============================================================================

# ===== Class for Proton NERRD Dispersion Decay Fit ============================
class HNERRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,rHH,dwcsa,tauvalues,s2values,r2values,a0values,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
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
                                self.popt[1],self.popt[2],self.popt[3]) + random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                             fitresult.params['R20'].value,fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=mean(tmpfit,0)
        self.mcfiterr2=std(tmpfit,0)
        self.tmpysim =ProtonR1rhoNERRDDecay1(x1,x2,x3,B0,wR,rHH,dwcsa,
                            mean(tmpfit,0)[0],mean(tmpfit,0)[1],mean(tmpfit,0)[2],mean(tmpfit,0)[3])
        self.tmpr1rhosim = ProtonR1rhoNERRD1(x4,x6,B0,wR,rHH,dwcsa,mean(tmpfit,0)[0],
                            mean(tmpfit,0)[1],mean(tmpfit,0)[2])

# ===== Class for Proton BMRD Dispersion Decay Fit ============================
class HBMRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,tauvalues,dwHexvalues,r2values,a0values,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        tau0,taumin,taumax = tauvalues[0],tauvalues[1],tauvalues[2]
        dwHex0,dwHexmin,dwHexmax = dwHexvalues[0],dwHexvalues[1],dwHexvalues[2]
        r200,r20min,r20max=r2values[0],r2values[1],r2values[2]
        a0,a0min,a0max = a0values[0],a0values[1],a0values[2]
        fitmodel = lmfit.Model(ProtonR1rhoBMRDDecay1,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('tau', value = tau0, min=taumin, max=taumax)
        fitmodel.set_param_hint('dwHex', value = dwHex0, min=dwHexmin, max=dwHexmax)
        fitmodel.set_param_hint('R20', value = r200, min=r20min, max=r20max)
        fitmodel.set_param_hint('a', value = a0, min=a0min, max=a0max)
        fitmodel.set_param_hint('B0',value=B0,vary=False)
        fitmodel.set_param_hint('wR',value=wR,vary=False)
        params = fitmodel.make_params()
        #print(params)
        fitresult=fitmodel.fit(logy, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
        self.popt =[fitresult.params['tau'].value,fitresult.params['dwHex'].value,
                  fitresult.params['R20'].value,fitresult.params['a'].value]
        self.perr = [fitresult.params['tau'].stderr,fitresult.params['dwHex'].stderr,
                  fitresult.params['R20'].stderr,fitresult.params['a'].stderr]
        # Monte-corlo fitting of decays and rf fields together
        tmpfit = []
        for j in range(MC):
            tmpysimref = (ProtonR1rhoBMRDDecay1(x1,x2,x3,B0,wR,self.popt[0],
                                self.popt[1],self.popt[2],self.popt[3]) + random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['dwHex'].value,
                             fitresult.params['R20'].value,fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=mean(tmpfit,0)
        self.mcfiterr2=std(tmpfit,0)
        self.tmpysim =ProtonR1rhoBMRDDecay1(x1,x2,x3,B0,wR,
                            mean(tmpfit,0)[0],mean(tmpfit,0)[1],mean(tmpfit,0)[2],mean(tmpfit,0)[3])
        self.tmpr1rhosim = ProtonR1rhoBMRD1(x4,x6,B0,wR,mean(tmpfit,0)[0],
                            mean(tmpfit,0)[1],mean(tmpfit,0)[2])

# ===== Class for Proton BMRD+NERRD Dispersion Decay Fit ============================
class HBMRDNERRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,rHH,dwcsa,tauvalues,s2values,dwHexvalues,r2values,
                 a0values,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        tau0,taumin,taumax = tauvalues[0],tauvalues[1],tauvalues[2]
        s20,s2min,s2max = s2values[0],s2values[1],s2values[2]
        dwHex0,dwHexmin,dwHexmax = dwHexvalues[0],dwHexvalues[1],dwHexvalues[2]
        r200,r20min,r20max=r2values[0],r2values[1],r2values[2]
        a0,a0min,a0max = a0values[0],a0values[1],a0values[2]
        fitmodel = lmfit.Model(ProtonR1rhoBMRDNERRDDecay1,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('tau', value = tau0, min=taumin, max=taumax)
        fitmodel.set_param_hint('S2', value = s20, min=s2min, max=s2max)
        fitmodel.set_param_hint('dwHex', value = dwHex0, min=dwHexmin, max=dwHexmax)
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
                    fitresult.params['dwHex'].value,fitresult.params['R20'].value,
                    fitresult.params['a'].value]
        self.perr = [fitresult.params['tau'].stderr,fitresult.params['S2'].stderr,
                     fitresult.params['dwHex'].stderr,fitresult.params['R20'].stderr,
                     fitresult.params['a'].stderr]
        # Monte-corlo fitting of decays and rf fields together
        tmpfit = []
        for j in range(MC):
            tmpysimref = (ProtonR1rhoBMRDNERRDDecay1(x1,x2,x3,B0,wR,rHH,dwcsa,self.popt[0],
                       self.popt[1],self.popt[2],self.popt[3],self.popt[4]) + random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                          fitresult.params['dwHex'].value,fitresult.params['R20'].value,
                          fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=mean(tmpfit,0)
        self.mcfiterr2=std(tmpfit,0)
        self.tmpysim =ProtonR1rhoBMRDNERRDDecay1(x1,x2,x3,B0,wR,rHH,dwcsa,
                                    mean(tmpfit,0)[0],mean(tmpfit,0)[1],mean(tmpfit,0)[2],
                                    mean(tmpfit,0)[3],mean(tmpfit,0)[4])
        self.tmpr1rhosim = ProtonR1rhoBMRDNERRD1(x4,x6,B0,wR,rHH,dwcsa,mean(tmpfit,0)[0],
                            mean(tmpfit,0)[1],mean(tmpfit,0)[2],mean(tmpfit,0)[3])


# ====== Single residueDispersion fitting and right model estimation for Proton ==============
class HDispFitsr:
    def __init__(self,rfouttodisp,k,fittype,wR,B0,rHH,dwHcsa,r20values,MC,pcutoff):
        s=rfouttodisp
        # x1, x2, x3 are triplets with rf, angle and time
        x1 = s.rfangletime[k,:,0] * sin(s.rfangletime[k,:,1])
        x2 = s.rfangletime[k,:,1]*180/pi # angle triplet
        x3 = s.rfangletime[k,:,2] # time triplet
        y  = s.normvolumes1[k,:]
        yerr = s.normnoises1[k,:]
        logy = log(abs(y))
        logyerr = abs(divide(yerr,y))
        tau0,s20,r200,a0,dwHex0 = 0.01,0.98,r20values[0],1,0.5
        taumin,s2min,r20min,a0min,dwHexmin = 0.00001,0.1,r20values[1],0.1,0.001
        taumax,s2max,r20max,a0max,dwHexmax = 0.5,0.999,r20values[2],10,20
        s2= NoDispDecayFit(x1,x2,x3,logy,logyerr,r200,r20min,r20values[3],a0,a0min,a0max,MC)
        # NERRD fitting
        if (fittype == "NERRD"):
            s3=HNERRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,rHH,dwHcsa,[tau0,taumin,taumax],
                                  [s20,s2min,s2max],[r200,r20min,r20max],[a0,a0min,a0max],MC)
        # BMRD fitting
        if (fittype == "BMRD"):
            s3=HBMRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,[tau0,taumin,taumax],
                                 [dwHex0,dwHexmin,dwHexmax],[r200,r20min,r20max],[a0,a0min,a0max],MC)
        # NERRD+BMRD fitting
        if (fittype == "BMRDNERRD"):
            s3=HBMRDNERRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,rHH,dwHcsa,[tau0,taumin,taumax],
                                      [s20,s2min,s2max],[dwHex0,dwHexmin,dwHexmax],[r200,r20min,r20max],
                                      [a0,a0min,a0max],MC)
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

# ====== All residues Dispersion fitting and right model estimation for Proton  ==============
class HDispFit:
    def __init__(self,rfouttodisp,fittype,wR,B0,rHHlist,dwHcsa,r20values,MC,pcutoff):
        s=rfouttodisp
        fitvalues1,fiterr1,mcfitvalues1,mcfiterr1 = [],[],[],[]
        reffitvalues1,reffiterr1,refmcfitvalues1,refmcfiterr1 = [],[],[],[]
        ftestresults,ysimref,ysim  = [],[],[]
        r1rhosimref,r1rhosim = [],[]
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        if array(rHHlist).size == 1:
            rHH = [rHHlist for i in range(s.nres)]
        else:
            rHH = array(rHHlist)
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
            self.r1rhosim = array(r1rhosim)
            self.r1rhosimref = array(r1rhosimref)
            self.ftestresults = array(ftestresults)
            self.rHH = rHH
        except IndexError:
            print("rHH list is not same length as no. of peakks")
#dispfitout = HDispFit(rfouttodisp,fittype,wR,B0,rHH,dwHcsa,MC,pcutoff)

# ====== Dispersion plots for proton ==============
class HDispPlotwithFit:
    def __init__(self,decayfitout,dispfitout,fittype,wR,MC,pcutoff,prefilename,ymax1):
        # If you don not need dispersion fit
        s1,s2 = ToDisp(decayfitout),dispfitout
        ncol  = 5
        nrows = (s1.nres//ncol)+1
        nrf = s1.nrf
        xmax = 1.15*amax(s1.effrf)
        xmin = -0.1*xmax
        if (ymax1 == -1 ):
            ymax = 1.15*amax(s1.r1rhocor,0)
        elif (ymax1 == 0):
            ymax = amax(s1.r1rhocor)
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
                    self.axes2[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes3[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes4[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
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
                    MySubPlot(self.axes2[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    MySubPlot(self.axes3[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    MySubPlot(self.axes4[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
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
#HDispPlotwithFit(decayfitout,dispfitout,fittype,wR,MC,pcutoff,prefilename,220)

class HDispPlotwithFit1:
    def __init__(self,decayfitout,dispfitout,fittype,wR,MC,pcutoff,prefilename,ymax1):
        # If you don not need dispersion fit
        s1,s2 = ToDisp(decayfitout),dispfitout
        ncol  = 5
        nrows = (s1.nres//ncol)+1
        nrf = s1.nrf
        xmax = 1.15*amax(s1.effrf)
        xmin = -0.1*xmax
        if (ymax1 == -1 ):
            ymax = 1.15*amax(s1.r1rhocor,0)
        elif (ymax1 == 0):
            ymax = amax(s1.r1rhocor)
        else:
            ymax = ymax1
        ymin = -0.1*ymax
        k=0
        filename2 = prefilename + '_' + fittype + 'nodispersion_withfit.pdf'
        filename3 = prefilename + '_' + fittype + 'dispersion_withfit.pdf'
        filename4 = prefilename + '_' + fittype + 'dispersion_withfit_ftest.pdf'
        filename5 = prefilename + '_' + fittype + 'dispersion_withfit_ftest_noMC.pdf'
        #  fig2,fig3,fig4
        self.fig2, self.axes2 = plt.subplots(nrows,ncol,figsize=(ncol*2.3,1.8*nrows))
        self.fig3, self.axes3 = plt.subplots(nrows,ncol,figsize=(ncol*2.3,1.8*nrows))
        self.fig4, self.axes4 = plt.subplots(nrows,ncol,figsize=(ncol*2.3,1.8*nrows))
        self.fig5, self.axes5 = plt.subplots(nrows,ncol,figsize=(ncol*2.3,1.8*nrows))
        for i in range(nrows):
            for j in range(ncol):
                if (k < s1.nres):
                    self.axes2[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes3[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes4[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes5[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
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
                    tmptxt4 = r'S$^2$=' + str(round(s2.mcfitvalues1[k][1],2)) \
                              + r' $\pm$ ' + str(round(s2.mcfiterr1[k][1],2))
                    tmptxt21 = r'R$_2^0$=' + str(round(s2.fitvalues1[k][2],2)) \
                              + r' $\pm$ ' + str(round(s2.fiterr1[k][2],2)) + r' s$^{-1}$'
                    tmptxt31 = r'$\tau_c$ =' + str(round((s2.fitvalues1[k][0])*1000,2)) \
                              + r' $\pm$ ' + str(round((s2.fiterr1[k][0])*1000,2)) + r' $\mu$s'
                    tmptxt41 = r'S$^2$=' + str(round(s2.fitvalues1[k][1],2)) \
                              + r' $\pm$ ' + str(round(s2.fiterr1[k][1],2))
                    tmptxt5 = r'r$_{HH}^{eff}$=' + str(round(s2.rHH[k],2)) + r' $\AA$'
                    self.axes2[i,j].text(-1,0.85*ymax,tmptxt1,fontsize = 6)
                    self.axes3[i,j].text(-1,0.85*ymax,tmptxt3,fontsize = 6)
                    self.axes3[i,j].text(-1,0.70*ymax,tmptxt4,fontsize = 6)
                    self.axes3[i,j].text(-1,0.55*ymax,tmptxt2,fontsize = 6)
                    self.axes3[i,j].text(-1,0.40*ymax,tmptxt5,fontsize = 6)
                    if (s2.ftestresults[k,0]==1):
                        self.axes4[i,j].plot(s2.effrf[:],s2.r1rhosim[k,:],color='green',linewidth=1.5)
                        self.axes4[i,j].text(-1,0.85*ymax,tmptxt3,fontsize = 6)
                        self.axes4[i,j].text(-1,0.70*ymax,tmptxt4,fontsize = 6)
                        self.axes4[i,j].text(-1,0.55*ymax,tmptxt2,fontsize = 6)
                        self.axes4[i,j].text(-1,0.40*ymax,tmptxt5,fontsize = 6)
                        self.axes5[i,j].plot(s2.effrf[:],s2.r1rhosim[k,:],color='green',linewidth=1.5)
                        self.axes5[i,j].text(-1,0.85*ymax,tmptxt31,fontsize = 6)
                        self.axes5[i,j].text(-1,0.70*ymax,tmptxt41,fontsize = 6)
                        self.axes5[i,j].text(-1,0.55*ymax,tmptxt21,fontsize = 6)
                        self.axes5[i,j].text(-1,0.40*ymax,tmptxt5,fontsize = 6)
                    else :
                        self.axes4[i,j].plot(s2.effrf[:],s2.r1rhosimref[k,:],color='red',linewidth=1.5)
                        self.axes4[i,j].text(-1,0.85*ymax,tmptxt1,fontsize = 6)
                        self.axes5[i,j].plot(s2.effrf[:],s2.r1rhosimref[k,:],color='red',linewidth=1.5)
                        self.axes5[i,j].text(-1,0.85*ymax,tmptxt1,fontsize = 6)
                    self.axes2[i,j].set_xlim([-5.0,xmax])
                    self.axes3[i,j].set_xlim([-5.0,xmax])
                    self.axes4[i,j].set_xlim([-5.0,xmax])
                    self.axes5[i,j].set_xlim([-5.0,xmax])
                    self.axes2[i,j].set_ylim([ymin,ymax])
                    self.axes3[i,j].set_ylim([ymin,ymax])
                    self.axes4[i,j].set_ylim([ymin,ymax])
                    self.axes5[i,j].set_ylim([ymin,ymax])
                    self.axes2[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes3[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes4[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes5[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes2[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    self.axes3[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    self.axes4[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    self.axes5[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    MySubPlot(self.axes2[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    MySubPlot(self.axes3[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    MySubPlot(self.axes4[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    MySubPlot(self.axes5[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    k=k+1
                else:
                    self.axes2[i,j].set_axis_off()
                    self.axes3[i,j].set_axis_off()
                    self.axes4[i,j].set_axis_off()
                    self.axes5[i,j].set_axis_off()
                    k=k+1
        # adjust space between the subplots
        self.fig2.subplots_adjust(wspace=0.4, hspace=0.9)
        self.fig3.subplots_adjust(wspace=0.4, hspace=0.9)
        self.fig4.subplots_adjust(wspace=0.4, hspace=0.9)
        self.fig5.subplots_adjust(wspace=0.4, hspace=0.9)
        # save plots to pdf
        self.fig2.savefig(filename2,dpi=1200)
        self.fig3.savefig(filename3,dpi=1200)
        self.fig4.savefig(filename4,dpi=1200)
        self.fig5.savefig(filename5,dpi=1200)
#HDispPlotwithFit(decayfitout,dispfitout,fittype,wR,MC,pcutoff,prefilename,220)

"""
# ====== Dispersion fitting and right model estimation for proton ==============
class HDispFit:
    def __init__(self,rfouttodisp,fittype,wR,B0,MC,pcutoff):
        s=rfouttodisp
        fitvalues1,fiterr1,mcfitvalues1,mcfiterr1 = [],[],[],[]
        reffitvalues1,reffiterr1,refmcfitvalues1,refmcfiterr1 = [],[],[],[]
        ftestresults,ysimref,ysim  = [],[],[]
        r1rhosimref,r1rhosim = [],[]
        for k in range(s.nres):
            # x1, x2, x3 are triplets with rf, angle and time
            x1 = s.rfangletime[k,:,0] # effrf triplet
            x2 = s.rfangletime[k,:,1]*180/pi # angle triplet
            x3 = s.rfangletime[k,:,2] # time triplet
            x4 = arange(0.2,array(s.effrf).max(),0.2)  # effrf singlet
            #x5 = array(s.angles)[:,i]*180/pi # tilt angle singlet
            x6 = empty(x4.size)
            x6.fill(90) # set all tilt angles to 90 degrees
            y  = s.normvolumes1[k,:]
            yerr = s.normnoises1[k,:]
            logy = log(abs(y))
            logyerr = abs(divide(yerr,y))
            tau0,s20,r200,a0,dwhex0 = 0.01,0.98,25,1,0.5
            taumin,s2min,r20min,a0min,dwhexmin = 0.0001,0.5,1,0.1,0.001
            taumax,s2max,r20max,a0max,dwhexmax = 0.5,1.0,200,10,20
            # No-Dispersion fitting
            poptref,pcovref = curve_fit(NoDispDecay,(x1,x2,x3),logy,p0 = [r200,a0],
                                      bounds = ([r20min,a0min], [r20max,a0max]))
            perrref = sqrt(diag(pcovref))
            reffitvalues1.append(poptref)
            reffiterr1.append(perrref)
            # Monte-corlo fitting of decays and rf fields together
            tmpfit = []
            for j in range(MC):
                tmpysimref = (NoDispDecay((x1,x2,x3),poptref[0],poptref[1]) + random.normal(0,logyerr))
                try:
                    poptsimref, pcovsimref = curve_fit(NoDispDecay,(x1,x2,x3),logy,
                                                     p0 = [r200,a0],bounds = ([r20min,a0min], [r20max,a0max]))
                    #print(poptsimref)
                    tmpfit.append(poptsimref)
                except RuntimeError:
                    print("Error - curve_fit failed for this MC iteration")
                refmcfitvalues1.append(mean(tmpfit,0))
                refmcfiterr1.append(std(tmpfit,0))
            tmpysimref = NoDispDecay((x1,x2,x3),mean(tmpfit,0)[0],mean(tmpfit,0)[1])
            ysimref.append(tmpysimref)
            #tmpr1rhosimref = NoDisp((x4,x5),mean(tmpfit,0)[0])
            tmpr1rhosimref = NoDisp((x4,x6),mean(tmpfit,0)[0])
            r1rhosimref.append(tmpr1rhosimref)

        # ===== Begining of Dispersion Fitting for proton ======================
            # Fitting for Proton NERRD
            if (fittype == "NERRD"):
                popt,pcov = curve_fit(ProtonR1rhoNERRDDecay,(x1,x2,x3),logy,p0 = [tau0,s20,r200,a0],
                                      bounds = ([taumin,s2min,r20min,a0min], [taumax,s2max,r20max,a0max]),maxfev=5000)
                #print(popt)
                perr = sqrt(diag(pcov))
                fitvalues1.append(popt)
                fiterr1.append(perr)
                tmpfit = []
                for j in range(MC):
                    tmpysim = (ProtonR1rhoNERRDDecay((x1,x2,x3),popt[0],popt[1],popt[2],popt[3])
                                            + random.normal(0,logyerr))
                    try:
                        poptsim, pcovsim = curve_fit(ProtonR1rhoNERRDDecay,(x1,x2,x3),logy,
                                                 p0 = [tau0,s20,r200,a0],
                                    bounds = ([taumin,s2min,r20min,a0min], [taumax,s2max,r20max,a0max]),maxfev=5000)
                        #print(poptsim)
                        tmpfit.append(poptsim)
                    except RuntimeError:
                        print("Error - curve_fit failed for this MC iteration")
                    mcfitvalues1.append(mean(tmpfit,0))
                    mcfiterr1.append(std(tmpfit,0))
                tmpysim = ProtonR1rhoNERRDDecay((x1,x2,x3),mean(tmpfit,0)[0],mean(tmpfit,0)[1],
                                                mean(tmpfit,0)[2],mean(tmpfit,0)[3])
                ysim.append(tmpysim)
                # Calculate fitted r1rho values
                #tmpr1rhosim = ProtonR1rhoNERRD((x4,x5),mean(tmpfit,0)[0],mean(tmpfit,0)[1],
                #                                mean(tmpfit,0)[2])
                tmpr1rhosim = ProtonR1rhoNERRD((x4,x6),mean(tmpfit,0)[0],mean(tmpfit,0)[1],
                                                mean(tmpfit,0)[2])
                r1rhosim.append(tmpr1rhosim)

            # Fitting for Proton BMRD
            if (fittype == "BMRD"):
                popt,pcov = curve_fit(ProtonR1rhoBMRDDecay,(x1,x2,x3),logy,p0 = [tau0,dwhex0,r200,a0],
                              bounds = ([taumin,dwhexmin,r20min,a0min], [taumax,dwhexmax,r20max,a0max]),maxfev=5000)
                #print(popt)
                perr = sqrt(diag(pcov))
                fitvalues1.append(popt)
                fiterr1.append(perr)
                tmpfit = []
                for j in range(MC):
                    tmpysim = (ProtonR1rhoBMRDDecay((x1,x2,x3),popt[0],popt[1],popt[2],popt[3])
                                            + random.normal(0,logyerr))
                    try:
                        poptsim, pcovsim = curve_fit(ProtonR1rhoBMRDDecay,(x1,x2,x3),logy,
                                                 p0 = [tau0,dwhex0,r200,a0], bounds =
                                        ([taumin,dwhexmin,r20min,a0min], [taumax,dwhexmax,r20max,a0max]),maxfev=5000)
                        #print(poptsim)
                        tmpfit.append(poptsim)
                    except RuntimeError:
                        print("Error - curve_fit failed for this MC iteration")
                    mcfitvalues1.append(mean(tmpfit,0))
                    mcfiterr1.append(std(tmpfit,0))
                tmpysim = ProtonR1rhoBMRDDecay((x1,x2,x3),mean(tmpfit,0)[0],mean(tmpfit,0)[1],
                                                mean(tmpfit,0)[2],mean(tmpfit,0)[3])
                ysim.append(tmpysim)
                # Calculate fitted r1rho values
                #tmpr1rhosim = ProtonR1rhoBMRD((x4,x5),mean(tmpfit,0)[0],mean(tmpfit,0)[1],
                #                                mean(tmpfit,0)[2])
                tmpr1rhosim = ProtonR1rhoBMRD((x4,x6),mean(tmpfit,0)[0],mean(tmpfit,0)[1],
                                                mean(tmpfit,0)[2])
                r1rhosim.append(tmpr1rhosim)

            # Fitting for Proton BMRD and NERRD together
            if (fittype == "BMRDNERRD"):
                popt,pcov = curve_fit(ProtonR1rhoBMRDNERRDDecay,(x1,x2,x3),logy,
                            p0=[tau0,s20,dwhex0,r200,a0],bounds = ([taumin,s2min,dwhexmin,r20min,a0min],
                                           [taumax,s2max,dwhexmax,r20max,a0max]),maxfev=5000)

                #print(popt)
                perr = sqrt(diag(pcov))
                fitvalues1.append(popt)
                fiterr1.append(perr)
                tmpfit = []
                for j in range(MC):
                    ysim = (ProtonR1rhoNERRDBMRDDecay((x1,x2,x3),popt[0],popt[1],popt[2],popt[3],popt[4])
                                            + random.normal(0,logyerr))
                    try:
                        poptsim, pcovsim = curve_fit(ProtonR1rhoBMRDDecay,(x1,x2,x3),logy,
                                                 p0 = [tau0,s20,dwhex0,r200,a0], bounds =
                             ([taumin,s2mindwhexmin,r20min,a0min], [taumax,s2max,dwhexmax,r20max,a0max]),maxfev=5000)
                        #print(poptsim)
                        tmpfit.append(poptsim)
                    except RuntimeError:
                        print("Error - curve_fit failed for this MC iteration")
                    mcfitvalues1.append(mean(tmpfit,0))
                    mcfiterr1.append(std(tmpfit,0))
                tmpysim = ProtonR1rhoNERRDBMRDDecay((x1,x2,x3),mean(tmpfit,0)[0],mean(tmpfit,0)[1],
                                                mean(tmpfit,0)[2],mean(tmpfit,0)[3],mean(tmpfit,0)[4])
                ysim.append(tmpysim)
                # Calculate fitted r1rho values
                #tmpr1rhosim = ProtonR1rhoNERRDBMRD((x4,x5),mean(tmpfit,0)[0],mean(tmpfit,0)[1],
                #                                mean(tmpfit,0)[2],mean(tmpfit,0)[3])
                tmpr1rhosim = ProtonR1rhoNERRDBMRD((x4,x6),mean(tmpfit,0)[0],mean(tmpfit,0)[1],
                                                mean(tmpfit,0)[2],mean(tmpfit,0)[3])
                r1rhosim.append(tmpr1rhosim)

            # ===== End of Dispersion Fitting for proton =======================
            # Run F-test for the fit
            tmpftestresults = FRatioToPvalue(logy,logyerr,tmpysim,tmpysimref,len(popt),
                                                 len(poptref),pcutoff)
            ftestresults.append(tmpftestresults)
        self.fitvalues1,self.fiterr1 = fitvalues1,fiterr1
        self.mcfitvalues1,self.mcfiterr1 = mcfitvalues1,mcfiterr1
        self.refmcfitvalues1,self.refmcfiterr1 = refmcfitvalues1,refmcfiterr1
        self.effrf = x4
        self.r1rhosim = array(r1rhosim)
        self.r1rhosimref = array(r1rhosimref)
        self.ftestresults = array(ftestresults)
"""

#===============================================================================
# ===== Carbon fitting functions ===============================================
#===============================================================================

# ===== Class for Carbon NERRD Dispersion Decay Fit ============================
class CNERRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,rHC,dwcsa,tauvalues,s2values,r2values,a0values,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        tau0,taumin,taumax = tauvalues[0],tauvalues[1],tauvalues[2]
        s20,s2min,s2max = s2values[0],s2values[1],s2values[2]
        r200,r20min,r20max=r2values[0],r2values[1],r2values[2]
        a0,a0min,a0max = a0values[0],a0values[1],a0values[2]
        fitmodel = lmfit.Model(CarbonR1rhoNERRDDecay1,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('tau', value = tau0, min=taumin, max=taumax)
        fitmodel.set_param_hint('S2', value = s20, min=s2min, max=s2max)
        fitmodel.set_param_hint('R20', value = r200, min=r20min, max=r20max)
        fitmodel.set_param_hint('a', value = a0, min=a0min, max=a0max)
        fitmodel.set_param_hint('B0',value=B0,vary=False)
        fitmodel.set_param_hint('wR',value=wR,vary=False)
        fitmodel.set_param_hint('rHC',value=rHC,vary=False)
        fitmodel.set_param_hint('dwC',value=dwcsa,vary=False)
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
            tmpysimref = (CarbonR1rhoNERRDDecay1(x1,x2,x3,B0,wR,rHC,dwcsa,self.popt[0],
                                self.popt[1],self.popt[2],self.popt[3]) + random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                             fitresult.params['R20'].value,fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=mean(tmpfit,0)
        self.mcfiterr2=std(tmpfit,0)
        self.tmpysim =CarbonR1rhoNERRDDecay1(x1,x2,x3,B0,wR,rHC,dwcsa,
                            mean(tmpfit,0)[0],mean(tmpfit,0)[1],mean(tmpfit,0)[2],mean(tmpfit,0)[3])
        self.tmpr1rhosim = CarbonR1rhoNERRD1(x4,x6,B0,wR,rHC,dwcsa,mean(tmpfit,0)[0],
                            mean(tmpfit,0)[1],mean(tmpfit,0)[2])

# ===== Class for Carbon BMRD Dispersion Decay Fit ============================
class CBMRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,tauvalues,dwCexvalues,r2values,a0values,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        tau0,taumin,taumax = tauvalues[0],tauvalues[1],tauvalues[2]
        dwCex0,dwCexmin,dwCexmax = dwCexvalues[0],dwCexvalues[1],dwCexvalues[2]
        r200,r20min,r20max=r2values[0],r2values[1],r2values[2]
        a0,a0min,a0max = a0values[0],a0values[1],a0values[2]
        fitmodel = lmfit.Model(CarbonR1rhoBMRDDecay1,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('tau', value = tau0, min=taumin, max=taumax)
        fitmodel.set_param_hint('dwCex', value = dwCex0, min=dwCexmin, max=dwCexmax)
        fitmodel.set_param_hint('R20', value = r200, min=r20min, max=r20max)
        fitmodel.set_param_hint('a', value = a0, min=a0min, max=a0max)
        fitmodel.set_param_hint('B0',value=B0,vary=False)
        fitmodel.set_param_hint('wR',value=wR,vary=False)
        params = fitmodel.make_params()
        #print(params)
        fitresult=fitmodel.fit(logy, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
        self.popt =[fitresult.params['tau'].value,fitresult.params['dwCex'].value,
                  fitresult.params['R20'].value,fitresult.params['a'].value]
        self.perr = [fitresult.params['tau'].stderr,fitresult.params['dwCex'].stderr,
                  fitresult.params['R20'].stderr,fitresult.params['a'].stderr]
        # Monte-corlo fitting of decays and rf fields together
        tmpfit = []
        for j in range(MC):
            tmpysimref = (CarbonR1rhoBMRDDecay1(x1,x2,x3,B0,wR,self.popt[0],
                                self.popt[1],self.popt[2],self.popt[3]) + random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['dwCex'].value,
                             fitresult.params['R20'].value,fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=mean(tmpfit,0)
        self.mcfiterr2=std(tmpfit,0)
        self.tmpysim =CarbonR1rhoBMRDDecay1(x1,x2,x3,B0,wR,
                            mean(tmpfit,0)[0],mean(tmpfit,0)[1],mean(tmpfit,0)[2],mean(tmpfit,0)[3])
        self.tmpr1rhosim = CarbonR1rhoBMRD1(x4,x6,B0,wR,mean(tmpfit,0)[0],
                            mean(tmpfit,0)[1],mean(tmpfit,0)[2])

# ===== Class for Carbon BMRD+NERRD Dispersion Decay Fit ============================
class CBMRDNERRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,rHC,dwcsa,tauvalues,s2values,dwCexvalues,r2values,
                 a0values,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        tau0,taumin,taumax = tauvalues[0],tauvalues[1],tauvalues[2]
        s20,s2min,s2max = s2values[0],s2values[1],s2values[2]
        dwCex0,dwCexmin,dwCexmax = dwCexvalues[0],dwCexvalues[1],dwCexvalues[2]
        r200,r20min,r20max=r2values[0],r2values[1],r2values[2]
        a0,a0min,a0max = a0values[0],a0values[1],a0values[2]
        fitmodel = lmfit.Model(CarbonR1rhoBMRDNERRDDecay1,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('tau', value = tau0, min=taumin, max=taumax)
        fitmodel.set_param_hint('S2', value = s20, min=s2min, max=s2max)
        fitmodel.set_param_hint('dwCex', value = dwCex0, min=dwCexmin, max=dwCexmax)
        fitmodel.set_param_hint('R20', value = r200, min=r20min, max=r20max)
        fitmodel.set_param_hint('a', value = a0, min=a0min, max=a0max)
        fitmodel.set_param_hint('B0',value=B0,vary=False)
        fitmodel.set_param_hint('wR',value=wR,vary=False)
        fitmodel.set_param_hint('rHC',value=rHC,vary=False)
        fitmodel.set_param_hint('dwC',value=dwcsa,vary=False)
        params = fitmodel.make_params()
        #print(params)
        fitresult=fitmodel.fit(logy, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
        self.popt =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                    fitresult.params['dwCex'].value,fitresult.params['R20'].value,
                    fitresult.params['a'].value]
        self.perr = [fitresult.params['tau'].stderr,fitresult.params['S2'].stderr,
                     fitresult.params['dwCex'].stderr,fitresult.params['R20'].stderr,
                     fitresult.params['a'].stderr]
        # Monte-corlo fitting of decays and rf fields together
        tmpfit = []
        for j in range(MC):
            tmpysimref = (CarbonR1rhoBMRDNERRDDecay1(x1,x2,x3,B0,wR,rHC,dwcsa,self.popt[0],
                       self.popt[1],self.popt[2],self.popt[3],self.popt[4]) + random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                          fitresult.params['dwCex'].value,fitresult.params['R20'].value,
                          fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=mean(tmpfit,0)
        self.mcfiterr2=std(tmpfit,0)
        self.tmpysim =CarbonR1rhoBMRDNERRDDecay1(x1,x2,x3,B0,wR,rHC,dwcsa,
                                    mean(tmpfit,0)[0],mean(tmpfit,0)[1],mean(tmpfit,0)[2],
                                    mean(tmpfit,0)[3],mean(tmpfit,0)[4])
        self.tmpr1rhosim = CarbonR1rhoBMRDNERRD1(x4,x6,B0,wR,rHC,dwcsa,mean(tmpfit,0)[0],
                            mean(tmpfit,0)[1],mean(tmpfit,0)[2],mean(tmpfit,0)[3])


# ====== Single residue Dispersion fitting and right model estimation for Carbon ==============
class CDispFitsr:
    def __init__(self,rfouttodisp,k,fittype,wR,B0,rHC,dwCcsa,MC,pcutoff):
        s=rfouttodisp
        # x1, x2, x3 are triplets with rf, angle and time
        x1 = s.rfangletime[k,:,0] * sin(s.rfangletime[k,:,1]) # effrf triplet
        x2 = s.rfangletime[k,:,1]*180/pi # angle triplet
        x3 = s.rfangletime[k,:,2] # time triplet
        y  = s.normvolumes1[k,:]
        yerr = s.normnoises1[k,:]
        logy = log(abs(y))
        logyerr = abs(divide(yerr,y))
        tau0,s20,r200,a0,dwCex0 = 0.01,0.98,25,1,0.5
        taumin,s2min,r20min,a0min,dwCexmin = 0.0001,0.5,1,0.1,0.001
        taumax,s2max,r20max,a0max,dwCexmax = 0.5,1.0,200,10,20
        s2= NoDispDecayFit(x1,x2,x3,logy,logyerr,r200,r20min,r20max,a0,a0min,a0max,MC)
        # NERRD fitting
        if (fittype == "NERRD"):
            s3=CNERRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,rHC,dwCcsa,[tau0,taumin,taumax],
                                  [s20,s2min,s2max],[r200,r20min,r20max],[a0,a0min,a0max],MC)
        # BMRD fitting
        if (fittype == "BMRD"):
            s3=CBMRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,[tau0,taumin,taumax],
                                 [dwCex0,dwCexmin,dwCexmax],[r200,r20min,r20max],[a0,a0min,a0max],MC)
        # NERRD+BMRD fitting
        if (fittype == "BMRDNERRD"):
            s3=CBMRDNERRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,rHC,dwCcsa,[tau0,taumin,taumax],
                                      [s20,s2min,s2max],[dwCex0,dwCexmin,dwCexmax],[r200,r20min,r20max],
                                      [a0,a0min,a0max],MC)
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

# ====== All residues Dispersion fitting and right model estimation for Carbon  ==============
class CDispFit:
    def __init__(self,rfouttodisp,fittype,wR,B0,rHC,dwCcsa,MC,pcutoff):
        s=rfouttodisp
        fitvalues1,fiterr1,mcfitvalues1,mcfiterr1 = [],[],[],[]
        reffitvalues1,reffiterr1,refmcfitvalues1,refmcfiterr1 = [],[],[],[]
        ftestresults,ysimref,ysim  = [],[],[]
        r1rhosimref,r1rhosim = [],[]
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        for k in range(s.nres):
            print(f"Residue {k} fit is running ...")
            s2=CDispFitsr(rfouttodisp,k,fittype,wR,B0,rHC,dwCcsa,MC,pcutoff)
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
        self.r1rhosim = array(r1rhosim)
        self.r1rhosimref = array(r1rhosimref)
        self.ftestresults = array(ftestresults)
#dispfitout = CDispFit(rfouttodisp,fittype,wR,B0,rHC,dwCcsa,MC,pcutoff)

# ====== Dispersion plots for Carbon ==============
class CDispPlotwithFit:
    def __init__(self,decayfitout,dispfitout,fittype,wR,MC,pcutoff,prefilename,ymax1):
        # If you don not need dispersion fit
        s1,s2 = ToDisp(decayfitout),dispfitout
        ncol  = 5
        nrows = (s1.nres//ncol)+1
        nrf = s1.nrf
        xmax = 1.15*amax(s1.effrf)
        xmin = -0.1*xmax
        if (ymax1 == -1 ):
            ymax = 1.15*amax(s1.r1rhocor,0)
        elif (ymax1 == 0):
            ymax = amax(s1.r1rhocor)
        else:
            ymax = ymax1
        ymin = -0.1*ymax
        ymax=30
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
                    self.axes2[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes3[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes4[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
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
                    tmptxt4 = r'S$^2$=' + str(round(s2.mcfitvalues1[k][1],3)) \
                              + r' $\pm$ ' + str(round(s2.mcfiterr1[k][1],3))
                    self.axes2[i,j].text(-1,0.85*ymax,tmptxt1,fontsize = 6)
                    self.axes3[i,j].text(-1,0.85*ymax,tmptxt3,fontsize = 6)
                    self.axes3[i,j].text(-1,0.70*ymax,tmptxt4,fontsize = 6)
                    self.axes3[i,j].text(-1,0.55*ymax,tmptxt2,fontsize = 6)
                    if (s2.ftestresults[k,0]==1):
                        self.axes4[i,j].plot(s2.effrf[:],s2.r1rhosim[k,:],color='green',linewidth=1.5)
                        self.axes4[i,j].text(-1,0.85*ymax,tmptxt3,fontsize = 6)
                        self.axes4[i,j].text(-1,0.70*ymax,tmptxt4,fontsize = 6)
                        self.axes4[i,j].text(-1,0.55*ymax,tmptxt2,fontsize = 6)
                    else :
                        self.axes4[i,j].plot(s2.effrf[:],s2.r1rhosimref[k,:],color='red',linewidth=1.5)
                        self.axes4[i,j].text(-1,0.85*ymax,tmptxt1,fontsize = 6)
                    self.axes2[i,j].set_xlim([-5.0,xmax])
                    self.axes3[i,j].set_xlim([-5.0,xmax])
                    self.axes4[i,j].set_xlim([-5.0,xmax])
                    self.axes2[i,j].set_ylim([ymin,ymax])
                    self.axes3[i,j].set_ylim([ymin,ymax])
                    self.axes4[i,j].set_ylim([ymin,ymax])
                    self.axes2[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes3[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes4[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes2[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    self.axes3[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    self.axes4[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    MySubPlot(self.axes2[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    MySubPlot(self.axes3[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    MySubPlot(self.axes4[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
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
#CDispPlotwithFit(decayfitout,dispfitout,fittype,wR,MC,pcutoff,prefilename,220)


#===============================================================================
# ===== Nitrogen fitting functions ===============================================
#===============================================================================

# ===== Class for Nitrogen NERRD Dispersion Decay Fit ============================
class NNERRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,rHN,dwNcsa,tauvalues,s2values,r2values,a0values,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        tau0,taumin,taumax = tauvalues[0],tauvalues[1],tauvalues[2]
        s20,s2min,s2max = s2values[0],s2values[1],s2values[2]
        r200,r20min,r20max=r2values[0],r2values[1],r2values[2]
        a0,a0min,a0max = a0values[0],a0values[1],a0values[2]
        fitmodel = lmfit.Model(NitrogenR1rhoNERRDDecay1,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('tau', value = tau0, min=taumin, max=taumax)
        fitmodel.set_param_hint('S2', value = s20, min=s2min, max=s2max)
        fitmodel.set_param_hint('R20', value = r200, min=r20min, max=r20max)
        fitmodel.set_param_hint('a', value = a0, min=a0min, max=a0max)
        fitmodel.set_param_hint('B0',value=B0,vary=False)
        fitmodel.set_param_hint('wR',value=wR,vary=False)
        fitmodel.set_param_hint('rHN',value=rHN,vary=False)
        fitmodel.set_param_hint('dwN',value=dwNcsa,vary=False)
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
            tmpysimref = (NitrogenR1rhoNERRDDecay1(x1,x2,x3,B0,wR,rHN,dwNcsa,self.popt[0],
                                self.popt[1],self.popt[2],self.popt[3]) + random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                             fitresult.params['R20'].value,fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=mean(tmpfit,0)
        self.mcfiterr2=std(tmpfit,0)
        self.tmpysim =NitrogenR1rhoNERRDDecay1(x1,x2,x3,B0,wR,rHN,dwNcsa,
                            mean(tmpfit,0)[0],mean(tmpfit,0)[1],mean(tmpfit,0)[2],mean(tmpfit,0)[3])
        self.tmpr1rhosim = NitrogenR1rhoNERRD1(x4,x6,B0,wR,rHN,dwNcsa,mean(tmpfit,0)[0],
                            mean(tmpfit,0)[1],mean(tmpfit,0)[2])

# ===== Class for Nitrogen BMRD Dispersion Decay Fit ============================
class NBMRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,tauvalues,dwNexvalues,r2values,a0values,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        tau0,taumin,taumax = tauvalues[0],tauvalues[1],tauvalues[2]
        dwCex0,dwCexmin,dwCexmax = dwCexvalues[0],dwCexvalues[1],dwCexvalues[2]
        r200,r20min,r20max=r2values[0],r2values[1],r2values[2]
        a0,a0min,a0max = a0values[0],a0values[1],a0values[2]
        fitmodel = lmfit.Model(NitrogenR1rhoBMRDDecay1,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('tau', value = tau0, min=taumin, max=taumax)
        fitmodel.set_param_hint('dwNex', value = dwNex0, min=dwNexmin, max=dwNexmax)
        fitmodel.set_param_hint('R20', value = r200, min=r20min, max=r20max)
        fitmodel.set_param_hint('a', value = a0, min=a0min, max=a0max)
        fitmodel.set_param_hint('B0',value=B0,vary=False)
        fitmodel.set_param_hint('wR',value=wR,vary=False)
        params = fitmodel.make_params()
        #print(params)
        fitresult=fitmodel.fit(logy, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
        self.popt =[fitresult.params['tau'].value,fitresult.params['dwNex'].value,
                  fitresult.params['R20'].value,fitresult.params['a'].value]
        self.perr = [fitresult.params['tau'].stderr,fitresult.params['dwNex'].stderr,
                  fitresult.params['R20'].stderr,fitresult.params['a'].stderr]
        # Monte-corlo fitting of decays and rf fields together
        tmpfit = []
        for j in range(MC):
            tmpysimref = (NitrogenR1rhoBMRDDecay1(x1,x2,x3,B0,wR,self.popt[0],
                                self.popt[1],self.popt[2],self.popt[3]) + random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['dwNex'].value,
                             fitresult.params['R20'].value,fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=mean(tmpfit,0)
        self.mcfiterr2=std(tmpfit,0)
        self.tmpysim = NitrogenR1rhoBMRDDecay1(x1,x2,x3,B0,wR,
                            mean(tmpfit,0)[0],mean(tmpfit,0)[1],mean(tmpfit,0)[2],mean(tmpfit,0)[3])
        self.tmpr1rhosim = NitrogenR1rhoBMRD1(x4,x6,B0,wR,mean(tmpfit,0)[0],
                            mean(tmpfit,0)[1],mean(tmpfit,0)[2])

# ===== Class for Nitrogen BMRD+NERRD Dispersion Decay Fit ============================
class NBMRDNERRDDispDecayFit:
    def __init__(self,x1,x2,x3,logy,logyerr,B0,wR,rHN,dwNcsa,tauvalues,s2values,dwNexvalues,r2values,
                 a0values,MC):
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        tau0,taumin,taumax = tauvalues[0],tauvalues[1],tauvalues[2]
        s20,s2min,s2max = s2values[0],s2values[1],s2values[2]
        dwNex0,dwNexmin,dwNexmax = dwNexvalues[0],dwNexvalues[1],dwNexvalues[2]
        r200,r20min,r20max=r2values[0],r2values[1],r2values[2]
        a0,a0min,a0max = a0values[0],a0values[1],a0values[2]
        fitmodel = lmfit.Model(NitrogenR1rhoBMRDNERRDDecay1,independent_vars=['w1','tilt','t'])
        fitmodel.set_param_hint('tau', value = tau0, min=taumin, max=taumax)
        fitmodel.set_param_hint('S2', value = s20, min=s2min, max=s2max)
        fitmodel.set_param_hint('dwNex', value = dwNex0, min=dwNexmin, max=dwNexmax)
        fitmodel.set_param_hint('R20', value = r200, min=r20min, max=r20max)
        fitmodel.set_param_hint('a', value = a0, min=a0min, max=a0max)
        fitmodel.set_param_hint('B0',value=B0,vary=False)
        fitmodel.set_param_hint('wR',value=wR,vary=False)
        fitmodel.set_param_hint('rHN',value=rHN,vary=False)
        fitmodel.set_param_hint('dwN',value=dwNcsa,vary=False)
        params = fitmodel.make_params()
        #print(params)
        fitresult=fitmodel.fit(logy, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
        self.popt =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                    fitresult.params['dwNex'].value,fitresult.params['R20'].value,
                    fitresult.params['a'].value]
        self.perr = [fitresult.params['tau'].stderr,fitresult.params['S2'].stderr,
                     fitresult.params['dwNex'].stderr,fitresult.params['R20'].stderr,
                     fitresult.params['a'].stderr]
        # Monte-corlo fitting of decays and rf fields together
        tmpfit = []
        for j in range(MC):
            tmpysimref = (NitrogenR1rhoBMRDNERRDDecay1(x1,x2,x3,B0,wR,rHN,dwNcsa,self.popt[0],
                       self.popt[1],self.popt[2],self.popt[3],self.popt[4]) + random.normal(0,logyerr))
            try:
                fitresult=fitmodel.fit(tmpysimref, params, w1=x1,tilt=x2,t=x3,weights=1./logyerr)
                poptsim =[fitresult.params['tau'].value,fitresult.params['S2'].value,
                          fitresult.params['dwNex'].value,fitresult.params['R20'].value,
                          fitresult.params['a'].value]
                tmpfit.append(poptsim)
            except RuntimeError:
                print("Error - curve_fit failed for this MC iteration")
        self.mcfitvalues2=mean(tmpfit,0)
        self.mcfiterr2=std(tmpfit,0)
        self.tmpysim = NitrogenR1rhoBMRDNERRDDecay1(x1,x2,x3,B0,wR,rHN,dwNcsa,
                                    mean(tmpfit,0)[0],mean(tmpfit,0)[1],mean(tmpfit,0)[2],
                                    mean(tmpfit,0)[3],mean(tmpfit,0)[4])
        self.tmpr1rhosim = NitrogenR1rhoBMRDNERRD1(x4,x6,B0,wR,rHN,dwNcsa,mean(tmpfit,0)[0],
                            mean(tmpfit,0)[1],mean(tmpfit,0)[2],mean(tmpfit,0)[3])

# ====== Single residue Dispersion fitting and right model estimation for Carbon ==============
class NDispFitsr:
    def __init__(self,rfouttodisp,k,fittype,wR,B0,rHN,dwNcsa,MC,pcutoff):
        s=rfouttodisp
        # x1, x2, x3 are triplets with rf, angle and time
        x1 = s.rfangletime[k,:,0] * sin(s.rfangletime[k,:,1])
        x2 = s.rfangletime[k,:,1]*180/pi # angle triplet
        x3 = s.rfangletime[k,:,2] # time triplet
        y  = s.normvolumes1[k,:]
        yerr = s.normnoises1[k,:]
        logy = log(abs(y))
        logyerr = abs(divide(yerr,y))
        tau0,s20,r200,a0,dwNex0 = 0.01,0.98,25,1,0.5
        taumin,s2min,r20min,a0min,dwNexmin = 0.0001,0.5,1,0.1,0.001
        taumax,s2max,r20max,a0max,dwNexmax = 0.5,1.0,200,10,20
        s2= NoDispDecayFit(x1,x2,x3,logy,logyerr,r200,r20min,r20max,a0,a0min,a0max,MC)
        # NERRD fitting
        if (fittype == "NERRD"):
            s3=NNERRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,rHN,dwNcsa,[tau0,taumin,taumax],
                                  [s20,s2min,s2max],[r200,r20min,r20max],[a0,a0min,a0max],MC)
        # BMRD fitting
        if (fittype == "BMRD"):
            s3=NBMRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,[tau0,taumin,taumax],
                                 [dwNex0,dwNexmin,dwNexmax],[r200,r20min,r20max],[a0,a0min,a0max],MC)
        # NERRD+BMRD fitting
        if (fittype == "BMRDNERRD"):
            s3=NBMRDNERRDDispDecayFit(x1,x2,x3,logy,logyerr,B0,wR,rHN,dwNcsa,[tau0,taumin,taumax],
                                      [s20,s2min,s2max],[dwNex0,dwNexmin,dwNexmax],[r200,r20min,r20max],
                                      [a0,a0min,a0max],MC)
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

# ====== All residues Dispersion fitting and right model estimation for Carbon  ==============
class NDispFit:
    def __init__(self,rfouttodisp,fittype,wR,B0,rHN,dwNcsa,MC,pcutoff):
        s=rfouttodisp
        fitvalues1,fiterr1,mcfitvalues1,mcfiterr1 = [],[],[],[]
        reffitvalues1,reffiterr1,refmcfitvalues1,refmcfiterr1 = [],[],[],[]
        ftestresults,ysimref,ysim  = [],[],[]
        r1rhosimref,r1rhosim = [],[]
        x4 = arange(0.2,150.0,0.2)  # effrf singlet
        x6 = array([90.0 for i in x4]) #angle singlet
        for k in range(s.nres):
            print(f"Residue {k} fit is running ...")
            s2=NDispFitsr(rfouttodisp,k,fittype,wR,B0,rHN,dwNcsa,MC,pcutoff)
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
        self.r1rhosim = array(r1rhosim)
        self.r1rhosimref = array(r1rhosimref)
        self.ftestresults = array(ftestresults)
#dispfitout = NDispFit(rfouttodisp,fittype,wR,B0,rHN,dwNcsa,MC,pcutoff)

# ====== Dispersion plots for Nitrogen ==============
class NDispPlotwithFit:
    def __init__(self,decayfitout,dispfitout,fittype,wR,MC,pcutoff,prefilename,ymax1):
        # If you don not need dispersion fit
        s1,s2 = ToDisp(decayfitout),dispfitout
        ncol  = 5
        nrows = (s1.nres//ncol)+1
        nrf = s1.nrf
        xmax = 1.15*amax(s1.effrf)
        xmin = -0.1*xmax
        if (ymax1 == -1 ):
            ymax = 1.15*amax(s1.r1rhocor,0)
        elif (ymax1 == 0):
            ymax = amax(s1.r1rhocor)
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
                    self.axes2[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes3[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
                                markeredgecolor = 'gold',markeredgewidth=0.5,linestyle = '',
                                markersize=markersize,ecolor='orange',capsize=capsize,lw=capsize/2)
                    self.axes4[i,j].errorbar(array(s1.effrf)[:,k],array(s1.r1rhocor)[:,k],
                                yerr=array(s1.r1rhoerrcor)[:,k],marker='o', color='darkmagenta',
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
                    tmptxt4 = r'S$^2$=' + str(round(s2.mcfitvalues1[k][1],2)) \
                              + r' $\pm$ ' + str(round(s2.mcfiterr1[k][1],2))
                    self.axes2[i,j].text(-1,0.85*ymax,tmptxt1,fontsize = 6)
                    self.axes3[i,j].text(-1,0.85*ymax,tmptxt3,fontsize = 6)
                    self.axes3[i,j].text(-1,0.70*ymax,tmptxt4,fontsize = 6)
                    self.axes3[i,j].text(-1,0.55*ymax,tmptxt2,fontsize = 6)
                    if (s2.ftestresults[k,0]==1):
                        self.axes4[i,j].plot(s2.effrf[:],s2.r1rhosim[k,:],color='green',linewidth=1.5)
                        self.axes4[i,j].text(-1,0.85*ymax,tmptxt3,fontsize = 6)
                        self.axes4[i,j].text(-1,0.70*ymax,tmptxt4,fontsize = 6)
                        self.axes4[i,j].text(-1,0.55*ymax,tmptxt2,fontsize = 6)
                    else :
                        self.axes4[i,j].plot(s2.effrf[:],s2.r1rhosimref[k,:],color='red',linewidth=1.5)
                        self.axes4[i,j].text(-1,0.85*ymax,tmptxt1,fontsize = 6)
                    self.axes2[i,j].set_xlim([-5.0,xmax])
                    self.axes3[i,j].set_xlim([-5.0,xmax])
                    self.axes4[i,j].set_xlim([-5.0,xmax])
                    self.axes2[i,j].set_ylim([ymin,ymax])
                    self.axes3[i,j].set_ylim([ymin,ymax])
                    self.axes4[i,j].set_ylim([ymin,ymax])
                    self.axes2[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes3[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes4[i,j].axvline(x=wR/2,color='dodgerblue',linestyle ='dotted',linewidth=0.5)
                    self.axes2[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    self.axes3[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    self.axes4[i,j].axvline(x=wR,color='dodgerblue',linestyle = 'dotted',linewidth=0.5)
                    MySubPlot(self.axes2[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    MySubPlot(self.axes3[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
                          s1.seqname[k])
                    MySubPlot(self.axes4[i,j],r'$\omega_{rf}$ (in kHz)',r'$R_{1rho}^{on}$ ($s^{-1}$)',
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
#NDispPlotwithFit(decayfitout,dispfitout,fittype,wR,MC,pcutoff,prefilename,220)


# ===== Settings for subplots ==================================================
def MySubPlot(ax,xlabel,ylabel,title):
    ax.set_xlabel(xlabel,fontsize=8)
    ax.set_ylabel(ylabel,fontsize=8)
    ax.set_title(title,fontsize=10)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.minorticks_on()
    return

# ===== Decay plots  ===========================================================
import numpy as np
class DecayPlots:
    def __init__(self, rfout, MC, fitplotflag, ncpucore):
        print("Decay plots are running ...")
        print(rfout, MC, fitplotflag, ncpucore)
        s = rfout
        ncol = 5
        nrows = (s.nres // ncol) + 1

        x = s.timelist
        xmax = 1.15 * np.amax(s.timelist)
        xt = np.arange(0, xmax)
        self.fig, self.axes = plt.subplots(nrows, ncol, figsize=(ncol * 2.3, 1.8 * nrows))

        if fitplotflag == 'yes':
            fitout = DecayFit(rfout, MC, ncpucore)
            self.fitvalues = fitout.fitvalues
            self.fiterr = fitout.fiterr
            self.mcfitvalues = fitout.mcfitvalues
            self.mcfiterr = fitout.mcfiterr

        max_normvolume = np.amax(s.normvolumes[0])
        min_normvolume = np.amin(s.normvolumes[0])

        for k in range(s.nres):
            i, j = divmod(k, ncol)
            self.axes[i, j].errorbar(s.timelist, s.normvolumes[k], yerr=abs(s.normnoises[k]), marker='o',
                                     color='green', linestyle='', markersize=markersize, ecolor='black',
                                     capsize=capsize, lw=capsize)
            if fitplotflag == 'yes':
                a0 = fitout.mcfitvalues[k, 0]
                R0 = fitout.mcfitvalues[k, 1]
                a0err = fitout.mcfiterr[k, 0]
                R0err = fitout.mcfiterr[k, 1]
                self.axes[i, j].plot(xt, MonoExp(xt, a0, R0 / 1000), 'r-')
                tmptxt = f"R1\u03C1={round(R0, 2)}\u00B1{round(R0err, 2)}"
                self.axes[i, j].text(0.2 * xmax, 1.3, tmptxt, fontsize=7)

            self.axes[i, j].set_xlim([-10.0, 1.2 * xmax])
            self.axes[i, j].set_ylim([0, 1.5 * max_normvolume])
            MySubPlot(self.axes[i, j], 'Delay(in ms.)', 'rel.int.', s.seqname[k])

        # Turn off any remaining empty subplots
        for k in range(s.nres, nrows * ncol):
            i, j = divmod(k, ncol)
            self.axes[i, j].set_axis_off()

        self.fig.subplots_adjust(wspace=0.4, hspace=0.9)
        self.fig.savefig(s.prefilename + 'decay.pdf')
        plt.close()
        #                 self.axes[i].errorbar(s.timelist,s.normvolumes[k], yerr=abs(s.normnoises[k]),marker='o',
        #                 color ='green',linestyle = '',markersize=markersize,ecolor = 'black',
        #                 capsize =capsize,lw=capsize)
        #                 if (fitplotflag == 'yes'):
        #                     a0 = fitout.mcfitvalues[k,0]
        #                     R0 = fitout.mcfitvalues[k,1]
        #                     #print(a0,R0)
        #                     a0err=fitout.mcfiterr[k,0]
        #                     R0err=fitout.mcfiterr[k,1]
        #                     self.axes[i].plot(xt,MonoExp(xt,a0,R0/1000),'r-')
        #                 self.axes[i].set_xlim([-10.0,1.2*amax(s.timelist)])
        #                 self.axes[i].set_ylim([0.0,1.5*amax(s.normvolumes[0])])
        #                 MySubPlot(self.axes[i],'Delay(in ms.)','rel.int.',s.seqname[k])
        #                 if (fitplotflag == 'yes'):
        #                     tmptxt = "R1\u03C1="+str(round(R0,2))+"\u00B1"+str(round(R0err,2))
        #                     self.axes[i].text(0.2*xmax,1.3,tmptxt,fontsize = 7)
        #                 k=k+1
        #         else:
        #             self.axes[i].set_axis_off()
        #             k=k+1
        #             self.fig.subplots_adjust(wspace=0.4, hspace=0.9)
        #             self.fig.savefig(s.prefilename+'decay.pdf')
        # plt.close()

