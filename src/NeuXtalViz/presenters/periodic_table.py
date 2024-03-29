class PeriodicTable:

    def __init__(self, view, model):

        self.view = view
        self.model = model

        self.view.connect_atoms(self.show_atom_dialog)

    def connect_atom_model_view(self, atom):

        view = self.view.get_atom_view()
        model = self.model.get_atom_model(atom)

        self.atom = Atom(view, model)
        self.view.value = self.model.value

    def show_atom_dialog(self, atom):

        self.connect_atom_model_view(atom)

        self.atom.view.connect_selected(self.update_selection)
        self.atom.view.show()

    def update_selection(self, data):

        self.model.value = data
        self.view.value = data

class Atom:

    def __init__(self, view, model):

        self.view = view
        self.model = model

        self.view.set_symbol_name(*self.model.get_symbol_name())
        self.view.set_isotope_numbers(self.model.get_isotope_numbers())
        self.set_info()

        self.view.connect_isotopes(self.set_info)
        self.view.connect_selection(self.set_isotope)

    def set_info(self):

        self.model.generate_data(self.view.get_isotope())
        self.view.set_atom_parameters(self.model.atom_dict,
                                      self.model.neutron_dict)

    def set_isotope(self):

        self.view.close()