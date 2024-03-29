class PeriodicTable:

    def __init__(self, view, model):

        self.view = view
        self.model = model

        self.view.connect_atoms(self.show_atom_dialog)

    def connect_atom_model_view(self):
        
        view = self.view.get_atom_view() 
        model = self.model.get_atom_model()

        self.atom = Atom(view, model)

    def show_atom_dialog(self, atom):

        print(atom)        

        self.connect_atom_model_view()

        self.atom.view.show()

class Atom:

    def __init__(self, view, model):

        self.view = view
        self.model = model