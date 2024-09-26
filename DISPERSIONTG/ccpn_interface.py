import os
import xml.etree.ElementTree as ET
import numpy as np
import re
#===== Begin: Definition of functions for CCPNMR loading =======================
#ccpnmrprojdir = '/Users/Suresh/Documents/home_workspace/SMILE/SH3_NUSrelax.ccpn'
#filename = os.listdir(ccpnmrprojdir+'/ccpnv3/ccp/nmr/Nmr/')[0]

class CCPNMRProject:
    def __init__(self, ccpnmrprojdir: str):
        """
        A class used to load CCPNMR projects.
        
        Args:
        -----
        ccpnmrprojdir (str): Path to the CCPNMR project directory.

        Attributes:
        ------------
        peaklistnames (list): List of peaklist names.
        spectranames (list): List of spectra names.
        heights (list): List of arrays with peak heights.
        volumes (list): List of arrays with peak volumes.
        seqname (list): List of sequence names.
        ndim (dict): Number of dimensions per spectrum.
        scale (dict): Scale factors per spectrum.
        sf1, sf2, sf3 (dict): Spectrometer frequencies for dimensions 1, 2, and 3.
        si1, si2, si3 (dict): Spectral size (number of points) for dimensions 1, 2, and 3.
        dw1, dw2, dw3 (dict): Spectral widths per point (dw) for dimensions 1, 2, and 3.
        ref1, ref2, ref3 (dict): Reference values for dimensions 1, 2, and 3.
        refpos1, refpos2, refpos3 (dict): Reference positions for dimensions 1, 2, and 3.
        assign1, assign2, assign3 (list): Peak assignments for dimensions 1, 2, and 3.
        pos1, pos2, pos3 (list): Peak positions for dimensions 1, 2, and 3.
        nspectra (int): Number of spectra in the project.
        """

        # Initialize empty lists and dictionaries for storing project information
        peaklistnames, spectranames, heights, volumes, seqname, residues = [], [], [], [], [], []
        si1, dw1, ref1, refpos1, sf1 = {}, {}, {}, {}, {}
        si2, dw2, ref2, refpos2, sf2 = {}, {}, {}, {}, {}
        si3, dw3, ref3, refpos3, sf3 = {}, {}, {}, {}, {}
        atomcodes1, atomcodes2, ndim, scale = {}, {}, {}, {}
        assign1, assign2, assign3, pos1, pos2, pos3 = [], [], [], [], [], []

        # Load the first NMR project file in the directory (assumes single file)
        filename = os.listdir(ccpnmrprojdir+'/ccpnv3/ccp/nmr/Nmr/')[0]
        tree = ET.parse(ccpnmrprojdir+'/ccpnv3/ccp/nmr/Nmr/' + filename)
        root = tree.getroot()

        # Define the XML paths for extracting the necessary information
        rootname1 = "NMR.NmrProject/NMR.NmrProject.experiments/NMR.Experiment"
        rootname2 = rootname1 + "/NMR.Experiment.dataSources/NMR.DataSource"

        # Find and store the spectrometer frequencies (sf1, sf2, sf3) for all spectra
        for experiment_elem in root.findall(rootname1):
            try:
                # Extract spectra name
                spectraname = experiment_elem.find("NMR.Experiment.dataSources/NMR.DataSource/NMR.DataSource.name/IMPL.Line").text
                ndim[spectraname] = int(experiment_elem.find("NMR.Experiment.dataSources/NMR.DataSource").attrib["numDim"])

                # Extract the spectrometer frequencies for each dimension
                for expdim_elem in experiment_elem.findall('NMR.Experiment.expDims/NMR.ExpDim'):
                    if expdim_elem.attrib["dim"] == '1':
                        sf1[spectraname] = float(expdim_elem.find("NMR.ExpDim.expDimRefs/NMR.ExpDimRef").attrib["sf"])
                    if expdim_elem.attrib["dim"] == '2':
                        sf2[spectraname] = float(expdim_elem.find("NMR.ExpDim.expDimRefs/NMR.ExpDimRef").attrib["sf"])
                    if expdim_elem.attrib["dim"] == '3':
                        sf3[spectraname] = float(expdim_elem.find("NMR.ExpDim.expDimRefs/NMR.ExpDimRef").attrib["sf"])
            except AttributeError:
                pass

        # Find residue names and atoms (residue and atom codes)
        residue_code_counter = 999
        for resonance_group_elem in root[0].find('NMR.NmrProject.resonanceGroups'):
            try:
                # Construct a residue identifier
                residue_id = resonance_group_elem.attrib['seqCode'] + resonance_group_elem.attrib['residueType']
            except KeyError:
                residue_id = str(residue_code_counter) + 'NAN'
                residue_code_counter -= 1
            residues.append(residue_id)

            # Get resonances for the residue
            try:
                resonance_ids = resonance_group_elem.find('NMR.ResonanceGroup.resonances').text.split()
            except AttributeError:
                resonance_ids = ['-999']
            resonance_to_residue_map = dict.fromkeys(resonance_ids, residue_id)
            atomcodes1.update(resonance_to_residue_map)

        # Find atom types and resonance peak assignments
        try:
            for resonance_elem in root[0].find('NMR.NmrProject.resonances'):
                resonance_id = resonance_elem.attrib['_ID']
                try:
                    atom_id = atomcodes1[resonance_id] + '-' + resonance_elem.attrib['implName']
                    peakdim_contrib_ids = resonance_elem.find('NMR.Resonance.peakDimContribs').text.split()
                except KeyError:
                    atom_id = atomcodes1[resonance_id] + '-' + resonance_elem.attrib['isotopeCode'][-1] + '?'
                except AttributeError:
                    peakdim_contrib_ids = ['-999']
                contrib_to_atom_map = dict.fromkeys(peakdim_contrib_ids, atom_id)
                atomcodes2.update(contrib_to_atom_map)
        except TypeError:
            atomcodes2.update({'-999': 'NAN'})

        # Find spectra and peaks
        for spectrum_elem in root.findall(rootname2):
            # Each spectrum_elem represents a spectrum in the CCPNMR project
            spectraname = spectrum_elem.find("NMR.DataSource.name/IMPL.Line").text
            spectranames.append(spectraname)

            # Get scale factors (if available)
            try:
                scale[spectraname] = float(spectrum_elem.attrib["scale"])
            except KeyError:
                scale[spectraname] = 1.0

            # Get frequency data (dw, si, ref values for dimensions 1, 2, 3)
            for freq_data_dim_elem in spectrum_elem.findall("NMR.DataSource.dataDims/NMR.FreqDataDim"):
                ref_value = freq_data_dim_elem.find("NMR.FreqDataDim.dataDimRefs/NMR.DataDimRef").attrib["refValue"]
                ref_point = freq_data_dim_elem.find("NMR.FreqDataDim.dataDimRefs/NMR.DataDimRef").attrib["refPoint"]

                # Store spectral information for each dimension
                if freq_data_dim_elem.attrib["dim"] == "1":
                    dw1[spectraname] = float(freq_data_dim_elem.attrib["valuePerPoint"])
                    si1[spectraname] = float(freq_data_dim_elem.attrib["numPoints"])
                    ref1[spectraname] = float(ref_value)
                    refpos1[spectraname] = float(ref_point)
                if freq_data_dim_elem.attrib["dim"] == "2":
                    dw2[spectraname] = float(freq_data_dim_elem.attrib["valuePerPoint"])
                    si2[spectraname] = float(freq_data_dim_elem.attrib["numPoints"])
                    ref2[spectraname] = float(ref_value)
                    refpos2[spectraname] = float(ref_point)
                if freq_data_dim_elem.attrib["dim"] == "3":
                    dw3[spectraname] = float(freq_data_dim_elem.attrib["valuePerPoint"])
                    si3[spectraname] = float(freq_data_dim_elem.attrib["numPoints"])
                    ref3[spectraname] = float(ref_value)
                    refpos3[spectraname] = float(ref_point)

            # For each spectrum, find peak lists and load peak data
            for peak_list_elem in spectrum_elem.findall("NMR.DataSource.peakLists/NMR.PeakList"):
                peaklist_num = peak_list_elem.attrib["serial"]
                peak_heights, peak_volumes = [], []
                peak_pos_dim1, peak_pos_dim2, peak_pos_dim3 = [], [], []
                peak_assign_dim1, peak_assign_dim2, peak_assign_dim3 = [], [], []

                # For each peak, extract height, volume, positions, and assignments
                for peak_elem in peak_list_elem.iter("NMR.Peak"):
                    peak_heights.append(float(peak_elem.attrib["height"]))

                    # Volume may not always be present, use default if not
                    try:
                        peak_volumes.append(float(peak_elem.attrib["volume"]))
                    except KeyError:
                        peak_volumes.append(-999.999)

                    # Extract peak positions and assignments for each dimension
                    for peak_dim_elem in peak_elem.findall("NMR.Peak.peakDims/NMR.PeakDim"):
                        if peak_dim_elem.attrib["dim"] == '1':
                            peak_pos_dim1.append(float(peak_dim_elem.attrib["position"]))
                            try:
                                peak_assign_dim1.append(atomcodes2[peak_dim_elem.find(
                                    'NMR.PeakDim.peakDimContribs/NMR.PeakDimContrib').attrib['_ID']])
                            except AttributeError:
                                peak_assign_dim1.append(atomcodes2.setdefault('-999', '-999'))
                        if peak_dim_elem.attrib["dim"] == '2':
                            peak_pos_dim2.append(float(peak_dim_elem.attrib["position"]))
                            try:
                                peak_assign_dim2.append(atomcodes2[peak_dim_elem.find(
                                    'NMR.PeakDim.peakDimContribs/NMR.PeakDimContrib').attrib['_ID']])
                            except AttributeError:
                                peak_assign_dim2.append(atomcodes2.setdefault('-999', '-999'))
                        if peak_dim_elem.attrib["dim"] == '3':
                            peak_pos_dim3.append(float(peak_dim_elem.attrib["position"]))
                            try:
                                peakdim_contribs = peak_dim_elem.findall('NMR.PeakDim.peakDimContribs/NMR.PeakDimContrib')
                                atom_assignments = ' '.join([atomcodes2[contrib.attrib['_ID']] for contrib in peakdim_contribs])
                                peak_assign_dim3.append(atom_assignments)
                            except AttributeError:
                                peak_assign_dim3.append(atomcodes2['-999'])

                # Store peak data for the current peak list
                peaklistnames.append(spectraname + "." + peaklist_num)
                heights.append(np.array(peak_heights))
                volumes.append(np.array(peak_volumes))
                assign1.append(peak_assign_dim1)
                assign2.append(peak_assign_dim2)
                assign3.append(peak_assign_dim3)
                pos1.append(np.array(peak_pos_dim1))
                pos2.append(np.array(peak_pos_dim2))
                pos3.append(np.array(peak_pos_dim3))

        # Set class attributes with the loaded data
        self.peaklistnames, self.spectranames = peaklistnames, spectranames
        self.heights = heights
        self.volumes = volumes
        self.seqname = seqname
        self.ndim, self.scale = ndim, scale
        self.ref1, self.dw1, self.si1, self.pos1, self.sf1, self.assign1 = ref1, dw1, si1, pos1, sf1, assign1
        self.ref2, self.dw2, self.si2, self.pos2, self.sf2, self.assign2 = ref2, dw2, si2, pos2, sf2, assign2
        self.ref3, self.dw3, self.si3, self.pos3, self.sf3, self.assign3 = ref3, dw3, si3, pos3, sf3, assign3
        self.refpos1, self.refpos2, self.refpos3 = refpos1, refpos2, refpos3
        self.nspectra = len(peaklistnames)

    def __str__(self):
        return f"CCPNMRProject Object with: {self.nspectra} spectra"


class CCPNMRPeakList:
    """
    A class to load, process, and organize peak data from a CCPNMR project for a specific peak list.

    This class retrieves peak assignments, volumes, heights, noise levels, and frequency offsets 
    from a CCPNMR dataset for a given spectrum file. It also supports multiple dimensions and spin systems.

    Attributes:
    -----------
    seqname (list): 
        List of sequence names corresponding to the peaks.
    seqnum (list): 
        List of sequence numbers corresponding to the peaks.
    noises (list): 
        List of noise levels for each residue.
    heights (list): 
        List of peak heights.
    offsets (list): 
        List of frequency offsets for each peak.
    timevalue (float): 
        Decay time value corresponding to the peak list.
    spectraname (str): 
        Name of the spectrum associated with the peak list.
    ndim (int): 
        Number of dimensions of the spectra.
    volumes (array): 
        Array of normalized peak volumes.
    freq1, freq2, freq3 (array): 
        Arrays containing frequencies for each dimension (only freq1 and freq2 for 2D, all three for 3D).
    assign1, assign2, assign3 (array): 
        Peak assignments for each dimension (only assign1 and assign2 for 2D, all three for 3D).
    nres (int): 
        Number of residues (peak sequences) in the data.

    Methods:
    --------
    __init__(self, ccpnmrdata, filename, timevalue, noisevalue, nscans, centerppm, spin):
        Initializes the class by loading the peak list from the CCPNMR data and extracting peak information.

    Parameters:
    -----------
    ccpnmrdata : object
        Data object containing the CCPNMR project information.
    filename : str
        Name of the peak list file to be loaded.
    timevalue : float
        Decay time value corresponding to the peak list.
    noisevalue : float
        Noise level corresponding to the peak list.
    nscans : int
        Number of scans used in the experiment.
    centerppm : float
        Center ppm value for the peak assignments.
    spin : str
        Spin system, either '1H' or another nucleus (e.g., '13C', '15N').

    Raises:
    -------
    IndexError:
        If the peak list does not exist in the CCPNMR project.

    Example:
    --------
    >>> ccpnmrdata = loadccpnmrproject("project_directory")
    >>> filename = "spectrum1.xml"
    >>> timevalue = 0.1
    >>> noisevalue = 0.02
    >>> nscans = 16
    >>> centerppm = 4.7
    >>> spin = '1H'
    >>> peaklist = CCPNMRPeakList(ccpnmrdata, filename, timevalue, noisevalue, nscans, centerppm, spin)
    """

    def __init__(self, ccpnmrdata, filename, timevalue, noisevalue, nscans, centerppm, spin):
        self.seqname = []
        self.seqnum = []
        self.noises = []
        self.heights = []
        self.offsets = []
        self.timevalue = timevalue
        try:
            # Locate the index of the peak list in the project data
            peaklist_index = [ccpnmrdata.peaklistnames.index(i)
                              for i in ccpnmrdata.peaklistnames if filename in i][0]
            print(filename)
            self.spectraname = filename[:-2]  # Removing the file extension
            self.ndim = ccpnmrdata.ndim[self.spectraname]
            
            # Get sequence names based on spin system
            if spin == '1H':
                sequence_names = np.array(ccpnmrdata.assign1[peaklist_index])
            else:
                sequence_names = np.array(ccpnmrdata.assign2[peaklist_index])
            
            # Extract sequence numbers from sequence names
            try:
                sequence_numbers = np.array([int(re.findall(r'\d+', name)[0]) for name in sequence_names])
            except IndexError:
                sequence_numbers = np.array([i for i in range(len(sequence_names))])

            # Sort the sequences by sequence numbers
            sorted_indices = sequence_numbers.argsort()
            self.seqname = sequence_names[sorted_indices]
            self.nres = len(self.seqname)
            self.seqnum = sequence_numbers[sorted_indices]

            # Normalize volumes and heights by the number of scans
            self.volumes = (ccpnmrdata.scale[self.spectraname]) * (1 / nscans) * \
                           (ccpnmrdata.volumes[peaklist_index][sorted_indices])
            self.heights = (ccpnmrdata.scale[self.spectraname]) * (1 / nscans) * \
                           (ccpnmrdata.heights[peaklist_index][sorted_indices])
            self.noises = [noisevalue / nscans for _ in range(self.nres)]

            # Calculate frequency shifts (freq1 and freq2 for 2D, add freq3 for 3D)
            dw1 = ccpnmrdata.dw1[self.spectraname] / ccpnmrdata.sf1[self.spectraname]
            dw2 = ccpnmrdata.dw2[self.spectraname] / ccpnmrdata.sf2[self.spectraname]

            # Calculate reference points for frequency axes
            if ccpnmrdata.refpos1[self.spectraname] == 1:
                highref1 = ccpnmrdata.ref1[self.spectraname]
            else:
                highref1 = (ccpnmrdata.si1[self.spectraname] * dw1 / 2) + ccpnmrdata.ref1[self.spectraname]
            if ccpnmrdata.refpos2[self.spectraname] == 1:
                highref2 = ccpnmrdata.ref2[self.spectraname]
            else:
                highref2 = (ccpnmrdata.si2[self.spectraname] * dw2 / 2) + ccpnmrdata.ref2[self.spectraname]
            
            # Calculate frequencies for the first and second dimensions
            self.freq1 = highref1 - (ccpnmrdata.pos1[peaklist_index][sorted_indices] * dw1)
            self.freq2 = highref2 - (ccpnmrdata.pos2[peaklist_index][sorted_indices] * dw2)

            # Store peak assignments
            self.assign1 = np.array(ccpnmrdata.assign1[peaklist_index])[sorted_indices]
            self.assign2 = np.array(ccpnmrdata.assign2[peaklist_index])[sorted_indices]

            # Handle 3D spectra if applicable
            if self.ndim == 3:
                dw3 = ccpnmrdata.dw3[self.spectraname] / ccpnmrdata.sf3[self.spectraname]
                if ccpnmrdata.refpos3[self.spectraname] == 1:
                    highref3 = ccpnmrdata.ref3[self.spectraname]
                else:
                    highref3 = (ccpnmrdata.si3[self.spectraname] * dw3 / 2) + ccpnmrdata.ref3[self.spectraname]
                self.freq3 = highref3 - (ccpnmrdata.pos3[peaklist_index][sorted_indices] * dw3)
                self.assign3 = np.array(ccpnmrdata.assign3[peaklist_index])[sorted_indices]

            # Calculate offsets based on spin system and frequency shifts
            if spin == '1H':
                self.offsets = 0.001 * ccpnmrdata.sf1[self.spectraname] * (self.freq1 - centerppm)
            else:
                self.offsets = 0.001 * ccpnmrdata.sf2[self.spectraname] * (self.freq2 - centerppm)
        except IndexError:
            print(f"Peaklist = '{filename}' does not exist in the CCPNMR project: Please Check")

    def __str__(self):
        return f"CCPNMRPeakList Object with: {self.nres} residues, {self.ndim} dimensions, and {self.timevalue} ms decay time"
    def __repr__(self):
        return self.__str__()




# ===== Class definition for multiple spectra peaks into a single object =======
class CCPNMRPeaksListGroup:
    """
    A class to load multiple CCPNMR peak lists, process and normalize the data.

    This class handles the loading of peak data from multiple files, processes volumes, 
    noise, and other spectral data, and computes derived properties such as normalized volumes, 
    effective RF frequency, and angles based on the given RF and offset values.

    Attributes:
    -----------
    rfkHz (float): 
        Radio frequency in kHz used during the NMR experiment.
    offset (float): 
        Frequency offset in kHz.
    prefilename (str): 
        The base filename, formatted with RF and offset values.
    volumes (list of arrays): 
        A list of arrays where each array contains the peak volumes from a peak list.
    noises (list of arrays): 
        A list of arrays where each array contains noise values from a peak list.
    assign1, assign2, assign3 (list): 
        Lists containing peak assignments for dimensions 1, 2, and 3 (if applicable).
    freq1, freq2, freq3 (list): 
        Lists containing frequencies for dimensions 1, 2, and 3 (if applicable).
    ndim (int): 
        Number of dimensions in the spectra.
    seqname (list): 
        Sequence names corresponding to the peaks.
    seqnum (list): 
        Sequence numbers corresponding to the peaks.
    offsets (array): 
        Array of offset values adjusted with the input offset.
    timelist (array): 
        Array of decay times used in the experiment.
    nres (int): 
        Number of residues (peak sequences) in the data.
    ntimes (int): 
        Number of time points (decay times).
    angles (array): 
        Array of angles computed based on the RF and offset values.
    normvolumes (array): 
        Normalized volumes computed as a ratio of volumes to the first time point.
    normnoises (array): 
        Normalized noise computed as a ratio of noise values to the first time point.
    effrfkhz (array): 
        Effective RF frequencies computed based on the RF and offset values.

    Methods:
    --------
    __init__(self, ccpnmrdata, filenamelist, rfkHz, offset, timelist, noisevaluelist, 
             nscanlist, centerppm, prefilename, spin):
        Initializes the class by loading and processing multiple CCPNMR peak lists.

    Parameters:
    -----------
    ccpnmrdata : object
        Data object containing the CCPNMR project information.
    filenamelist : list of str
        List of filenames corresponding to the peak lists to be loaded.
    rfkHz : float
        The radio frequency in kHz used in the experiment.
    offset : float
        The offset frequency in kHz applied during the experiment.
    timelist : list of float
        List of decay times corresponding to each peak list.
    noisevaluelist : list of float
        List of noise values for each peak list.
    nscanlist : list of int
        List of scan numbers corresponding to each peak list.
    centerppm : float
        Center ppm value for the peak assignments.
    prefilename : str
        Base prefix for filenames, to which RF and offset values are appended.
    spin : object
        Spin object containing spin system information.

    Raises:
    -------
    IndexError:
        If the number of files does not match the number of times, noise values, 
        or scan values, an IndexError is raised with an appropriate error message.

    Example:
    --------
    >>> ccpnmrdata = loadccpnmrproject("project_directory")
    >>> filenamelist = ["file1.xml", "file2.xml"]
    >>> rfkHz = 500.0
    >>> offset = 100.0
    >>> timelist = [0.1, 0.2]
    >>> noisevaluelist = [0.01, 0.02]
    >>> nscanlist = [16, 16]
    >>> centerppm = 4.7
    >>> prefilename = "experiment"
    >>> spin = spin_object
    >>> peakslist = loadccpnmrpeakslist(ccpnmrdata, filenamelist, rfkHz, offset, timelist,
                                        noisevaluelist, nscanlist, centerppm, prefilename, spin)
    """
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
                peaklist = CCPNMRPeakList(ccpnmrdata,filenamelist[i],timelist[i],
                                       noisevaluelist[i],nscanlist[i],centerppm,spin)
                self.volumes.append(peaklist.heights)
                self.noises.append(peaklist.noises)
            # for sub_array in self.volumes:
            #     print(len(sub_array))
            self.assign1 = peaklist.assign1
            self.assign2 = peaklist.assign2
            self.freq1 = peaklist.freq1
            self.freq2 = peaklist.freq2
            self.ndim = peaklist.ndim
            if (self.ndim == 3):
                self.assign3 = peaklist.assign3
                self.freq3 = peaklist.freq3
            self.seqname= peaklist.seqname
            self.seqnum = peaklist.seqnum
            self.offsets = peaklist.offsets + offset
            self.timelist = np.array(timelist)
            self.volumes = np.array(self.volumes)
            self.noises  = np.array(self.noises)
            self.nres    = len(self.seqnum)
            self.ntimes  = ntimes
            self.angles  = np.arctan(np.divide(rfkHz, self.offsets))
            self.normvolumes = np.transpose(np.divide(self.volumes,self.volumes[0]))
            self.normnoises  = np.transpose(np.divide(self.noises,self.volumes[0]))
            self.effrfkhz = np.sqrt(((rfkHz)**2) + ((self.offsets)**2))
        except IndexError:
            if not (nfiles == ntimes):
                print("Total no.of peaklist doesnot match with total no. of decay times: Please Check")
            if not (nfiles == len(noisevaluelist)):
                print("Total no.of peaklist doesnot match with total no. of noise values: Please Check")
            if not (nfiles == len(nscanlist)):
                print("Total no.of peaklist doesnot match with total no. of nscan values: Please Check")
    # function for printing
    def __str__(self):
        return f"CCPNMRPeaksListGroup Object with: Effective field {self.effrfkhz[0]:.2f}, Number of residues: {self.nres}, Number of time points: {self.ntimes}"
    

    def __repr__(self):
        return self.__str__()