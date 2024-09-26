import os
import re

def parse_title_file(folder_num):
    title_file_path = f"{folder_num}/pdata/1/title"
    acqus_file_path = f"{folder_num}/acqus"
    
    if not os.path.isfile(title_file_path) or not os.path.isfile(acqus_file_path):
        return None

    with open(title_file_path, 'r') as file:
        content = file.read()

    # Extracting data using regular expressions
    match = re.search(r'(\d+C-.+?)\s+with negative offresonance angle = ([\d.]+)', content)
    if match:
        experiment = match.group(1)
        offres_angle = float(match.group(2))
    else:
        return None

    match = re.search(r'pldb11 = ([\d.-]+) in dB', content)
    pldb11 = float(match.group(1)) if match else None

    match = re.search(r'cnst11 = ([\d.-]+) in Hz', content)
    cnst11 = float(match.group(1)) if match else None

    match = re.search(r'p11 = ([\d.]+) in us', content)
    p11 = float(match.group(1)) / 1000.0 if match else None  # Convert to milliseconds

    match = re.search(r'spinlock rfvalue = ([\d.]+) in kHz', content)
    spinlock_rfvalue = float(match.group(1)) if match else None

    match = re.search(r'spinlock eff.rf value = ([\d.]+) \(in kHz\)', content)
    eff_rf_value = float(match.group(1)) if match else None

    # Extracting nscans value from the acqus file
    with open(acqus_file_path, 'r') as file:
        acqus_content = file.read()

    match = re.search(r'##\$NS= (\d+)', acqus_content)
    nscans = int(match.group(1)) if match else None

    return {
        'folder_num': folder_num,
        'experiment': experiment,
        'offres_angle': offres_angle,
        'pldb11': pldb11,
        'cnst11': cnst11,
        'p11': p11,
        'spinlock_rfvalue': spinlock_rfvalue,
        'eff_rf_value': eff_rf_value,
        'nscans': nscans
    }

# Dictionary to hold data grouped by RF
data_by_rf = {}

# Process each folder
experiment_data = [range(12, 81), range(87, 104)]
noisevolumelist=[3.92e+06,7.96e+06, 8.59e+06 ,2.58e+07,1.59e+07,3.92e+06,8.47e+06,7.96e+06, 3.92e+06]
for i in experiment_data:
    data = parse_title_file(i)
    if data:
        eff_rf_key = f"RF_{str(data['eff_rf_value']).replace('.', '_')}"
        if eff_rf_key not in data_by_rf:
            data_by_rf[eff_rf_key] = {
                'filenamelist': [],
                'timelist': [],
                'noisevolumelist': noisevolumelist, # Placeholder value
                'nscanlist': [],
                'spinlockrf': data['spinlock_rfvalue'],
                'eff_rf_value': data['eff_rf_value'],
                'offset': data['cnst11']  # Placeholder value
            }
        data_by_rf[eff_rf_key]['filenamelist'].append(f"{data['folder_num']}-1.4")
        data_by_rf[eff_rf_key]['timelist'].append(data['p11'])
        # data_by_rf[eff_rf_key]['noisevolumelist'].append(7.96e+06)  # Placeholder value
        data_by_rf[eff_rf_key]['nscanlist'].append(data['nscans'])

# Print the results
counter = 1
for eff_rf_key, rf_data in data_by_rf.items():
    eff_rf_value = rf_data['eff_rf_value']
    print(f"# RF {counter} ({eff_rf_value} kHz)")
    print(f"filenamelist{counter} = {rf_data['filenamelist']}")
    print(f"timelist{counter} = {rf_data['timelist']}")
    print(f"noisevolumelist{counter} = {rf_data['noisevolumelist']}")
    print(f"nscanlist{counter} = {rf_data['nscanlist']}")
    print(f"spinlockrf{counter} = {rf_data['spinlockrf']}")
    print(f"offset{counter} = {rf_data['offset']}")
    print()
    counter += 1