import sys
import pymzml
from data_parser import load_massbank
import numpy as np

compounds = {}

def find_peaks(mzml_file, min_mz_diff=0.1, top_n=10):
    """
    Given the path of an mzml file, return a list of the top_n peaks.
    To prevent clustering, a minimum m/z difference is enforced, and
    each selected peak must be a local maximum to get a more diverse
    selection of readings on the spectrum.
    """
    ms_reader = pymzml.run.Reader(mzml_file)
    all_peaks = []

    for spectrum in ms_reader:
        if spectrum.ms_level == 1:  # MS1 spectra scans
            peaks = np.array(spectrum.peaks("centroided"))

            if len(peaks) == 0:
                continue

            mz, intensity = peaks[:, 0], peaks[:, 1]

            # Check if current peak is a local maximum
            local_maxima = np.r_[True, intensity[1:] > intensity[:-1]] & np.r_[intensity[:-1] > intensity[1:], True]
            local_max_peaks = peaks[local_maxima]

            all_peaks.extend(local_max_peaks)

    # Sort the peaks based on intensity
    all_peaks = np.array(sorted(all_peaks, key=lambda x: x[1], reverse=True))

    # Select the top N peaks while enforcing a minimum m/z difference
    selected_peaks = []
    for mz, intensity in all_peaks:
        if not any(abs(mz - mz2) < min_mz_diff for mz2, _ in selected_peaks):
            selected_peaks.append((mz, intensity))
            if len(selected_peaks) >= top_n:
                break

    # Normalize the intensities (0-1000) and sort by m/z
    selected_peaks = np.array(selected_peaks)
    max_intensity = np.max(selected_peaks[:, 1])
    normalized_peaks = np.column_stack((selected_peaks[:, 0], (selected_peaks[:, 1] / max_intensity) * 1000))
    normalized_peaks = normalized_peaks[np.argsort(normalized_peaks[:, 0])]
    
    return normalized_peaks.tolist()


def find_best_match(peaks, compounds):
    """
    Finds the compound with the highest matching similarity score given the peaks.
    Returns the matching compound's name, formula, molecular weight, and relative intensity.
    """
    best_score = 0
    best_match = None
    best_molecular_weight = 0
    best_relative_intensity = 0
    best_formula = ''

    for compound_name, compound_data in compounds.items():
        reference_peaks = compound_data['Peaks']
        score = cosine_similarity(peaks, reference_peaks)

        if score > best_score:
            best_score = score
            best_match = compound_name
            best_molecular_weight = float(compound_data['MW'])
            best_relative_intensity = max(peak[1] for peak in reference_peaks)
            best_formula = compound_data['Formula']

    return best_match, best_formula, best_molecular_weight, best_relative_intensity


def cosine_similarity(peaks1, peaks2, mz_tolerance=0.01, intensity_threshold=0.1):
    """
    Algorithm used for finding the similarity score between two sets of peaks
    """
    peaks1 = np.array([p for p in peaks1 if p[1] >= intensity_threshold])
    peaks2 = np.array([p for p in peaks2 if p[1] >= intensity_threshold])

    shared_intensity = 0.0
    intensity1_squared = 0.0
    intensity2_squared = 0.0

    i, j = 0, 0
    while i < len(peaks1) and j < len(peaks2):
        mz_diff = peaks1[i][0] - peaks2[j][0]
        if abs(mz_diff) <= mz_tolerance:
            shared_intensity += peaks1[i][1] * peaks2[j][1]
            intensity1_squared += peaks1[i][1] ** 2
            intensity2_squared += peaks2[j][1] ** 2
            i += 1
            j += 1
        elif mz_diff < 0:
            intensity1_squared += peaks1[i][1] ** 2
            i += 1
        else:
            intensity2_squared += peaks2[j][1] ** 2
            j += 1

    intensity1_squared += np.sum(peaks1[i:, 1] ** 2)
    intensity2_squared += np.sum(peaks2[j:, 1] ** 2)

    if intensity1_squared == 0 or intensity2_squared == 0:
        return 0.0

    return shared_intensity / ((intensity1_squared * intensity2_squared) ** 0.5)


def load_compounds():
    global compounds
    massbank_dir = "MassBank_NIST.msp"
    compounds = load_massbank(massbank_dir)


def main(mzml):
    global compounds
    peaks = find_peaks(mzml)
    # for peak in peaks:
    #     print(f'{peak[0]} {peak[1]}')
    return find_best_match(peaks, compounds)  


if __name__ == "__main__":
    """
    (Development Testing)
    """
    mzml_file = sys.argv[1]
    peaks = find_peaks(mzml_file)
    # for peak in peaks:
    #     print(f'{peak[0]} {peak[1]}')
    massbank_dir = "MassBank_NIST.msp"
    compounds = load_massbank(massbank_dir)
    print(find_best_match(peaks, compounds))
