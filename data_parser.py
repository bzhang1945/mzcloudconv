import xml.etree.ElementTree as ET
from pyteomics import mzml
from xml.dom import minidom
from lxml import etree

def load_massbank(filename):
    """
    Loads and processes the massbank MSP file. 
    Returns a dict of compounds- keys: compound name, value: dict of its specs
    """
    compounds = {}
    current_compound = None

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()

            if key == "Name":
                if current_compound:
                    compounds[current_compound["Name"]] = current_compound
                current_compound = {"Name": value}
            elif key == "Num Peaks":
                num_peaks = int(value)
                current_compound[key] = num_peaks
                peaks = [tuple(map(float, file.readline().strip().split())) for _ in range(num_peaks)]
                current_compound["Peaks"] = peaks
            else:
                current_compound[key] = value

    if current_compound:
        compounds[current_compound["Name"]] = current_compound

    return compounds


def parse_to_xml(file_path, output_path):
    """
    Parses the mzML file at 'file_path' into an output.xml file at 'output_path'
    Returns a list containing each element pertaining to the instrument parameters.
    """
    instrument_params = []
    with mzml.read(file_path, use_index=False) as reader:
        root = ET.Element('mzml_data')
        
        tree = etree.parse(file_path)
        ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}
        instrument_elem = ET.SubElement(root, 'instrument_information')
        for elem in tree.iterfind('*/mzml:referenceableParamGroupList/mzml:referenceableParamGroup', ns):
            ref_param_group = ET.SubElement(instrument_elem, 'referenceableParamGroup')
            ref_param_group.set('id', elem.get('id'))
            

            for cv in elem.iterfind('mzml:cvParam', ns):
                cv_param = ET.SubElement(ref_param_group, 'cvParam')
                cv_param.set('name', cv.get('name'))
                cv_param.set('value', cv.get('value'))
                instrument_params.append((cv.get('name'), cv.get('value')))
        
        scans_elem = ET.SubElement(root, 'scans')
        for spectrum in reader:
            scan = ET.SubElement(scans_elem, 'scan')
            try:
                scan.set('ms_level', str(spectrum['ms level']))
            except KeyError:
                break
            scan.set('mz_start', str(spectrum['lowest observed m/z']))
            scan.set('mz_end', str(spectrum['highest observed m/z']))
    
    with open(output_path, 'w', encoding='utf-8') as output_file:
        output_file.write(prettify(root))

    return instrument_params

def prettify(elem):
    """
    XML indentation helper function
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


def run_parser(file_path):
    """
    Returns a list containing each element pertaining to the instrument parameters.
    """
    output_path = '/tmp/output.xml'
    return parse_to_xml(file_path, output_path)


if __name__ == '__main__':
    """
    Development testing
    """
    print(run_parser())
