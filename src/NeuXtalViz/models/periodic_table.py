from mantid.kernel import Atom

from NeuXtalViz.config.atoms import isoptopes, names

class PeriodicTableModel():

    def __init__(self):

        pass


class AtomModel():

    def __init__(self, atm='H'):

        self.atm, self.isotopes = atm, isoptopes[atm]

        self.name = names[atm]

    def get_symbol(self):

        return self.atm

    def get_isotope_numbers(self):

        return self.isotopes

    def generate_data(self, iso):

        atom = Atom(self.atm, iso)

        self.atom_dict = {'mass_number': atom.a_number,
                          'abundace': atom.abundace,
                          'mass': atom.mass,
                          'mass_density': atom.mass_density,
                          'number_density': atom.number_density,
                          'z': atom.z_number}

        neutron = atom.neutron()

        self.neutron_dict = {'sigma_coh' : neutron['coh_scatt_xs'],
                             'sigma_inc' : neutron['inc_scatt_xs'],
                             'sigma_tot' : neutron['tot_scatt_xs'],
                             'sigma_abs' : neutron['abs_xs'],
                             'b_coh_re' : neutron['coh_scatt_length_real'],
                             'b_coh_im' : neutron['coh_scatt_length_img'],
                             'b_inc_re' : neutron['inc_scatt_length_real'],
                             'b_inc_im' : neutron['inc_scatt_length_img']}
