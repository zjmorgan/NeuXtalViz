from mantid.kernel import Atom

from NeuXtalViz.config.atoms import isotopes, names


class PeriodicTableModel:
    def __init__(self, atom):
        self.value = atom

    def get_atom_model(self, atm):
        return AtomModel(atm)


class AtomModel:
    def __init__(self, atm="H"):
        self.atm, self.isotopes = atm, isotopes.get(atm)

        self.name = names[atm]

    def get_symbol_name(self):
        return self.atm, self.name

    def get_isotope_numbers(self):
        return self.isotopes

    def generate_data(self, iso):
        atom = Atom(self.atm, iso)

        self.atom_dict = {
            "mass_number": atom.a_number,
            "abundance": atom.abundance,
            "mass": atom.mass,
            "z": atom.z_number,
        }

        neutron = atom.neutron()

        self.neutron_dict = {
            "sigma_coh": neutron["coh_scatt_xs"],
            "sigma_inc": neutron["inc_scatt_xs"],
            "sigma_tot": neutron["tot_scatt_xs"],
            "sigma_abs": neutron["abs_xs"],
            "b_coh_re": neutron["coh_scatt_length_real"],
            "b_coh_im": neutron["coh_scatt_length_img"],
            "b_inc_re": neutron["inc_scatt_length_real"],
            "b_inc_im": neutron["inc_scatt_length_img"],
        }
