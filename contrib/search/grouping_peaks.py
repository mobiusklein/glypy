

class DeconRow(object):
    def __init__(self, scan_number=0, charge=0, abundance=0.0, mz=0.0, fit=0.0,
                 average_mass_weight=0.0, monoisotopic_mass_weight=0.0, most_abundant_mass_weight=0.0,
                 full_width_half_height=0.0, signal_to_noise=0.0, monoisotopic_abundance=0.0,
                 monoisotopic_plus_2_abundance=0.0, flag=False, inference_score=0.0):
        self.scan_number = int(scan_number)
        self.charge = int(charge)
        self.abundance = abundance
        self.mz = mz
        self.fit = fit

        self.average_mass_weight = average_mass_weight
        self.monoisotopic_mass_weight = monoisotopic_mass_weight
        self.most_abundant_mass_weight = most_abundant_mass_weight

        self.full_width_half_height = full_width_half_height
        self.signal_to_noise = signal_to_noise

        self.monoisotopic_abundance = monoisotopic_abundance
        self.monoisotopic_plus_2_abundance = monoisotopic_plus_2_abundance

        self.flag = bool(flag)

        self.inference_score = inference_score

    @classmethod
    def from_csv(stream, noise_threshold=0):
        for line in stream:
            row = DeconRow(*map(float, line))
            yield row


class ResultsGroup(object):
    def __init__(self):
        pass
