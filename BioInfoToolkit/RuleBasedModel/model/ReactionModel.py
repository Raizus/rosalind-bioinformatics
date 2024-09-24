from collections import defaultdict
from functools import reduce
from itertools import accumulate, product
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Observable import Observable
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import BidirectionalReaction, ReactionRule
from BioInfoToolkit.RuleBasedModel.model.Pattern import ComplexReactant, Molecule, generate_species
import numpy as np


class Parameter:
    declaration: str
    name: str
    value: int|float

    def __init__(self, name: str, declaration: str, value: int|float) -> None:
        self.name = name
        self.declaration = declaration
        self.value = value

    def __repr__(self) -> str:
        out = f"{self.name} {self.declaration}"
        return out


class ReactionModel:
    molecules: dict[str, MoleculeType] = dict()
    defined_reactions: dict[str, ReactionRule] = dict()

    # this dictionary includes each reaction as a unidirectional reaction,
    # and also splits reactions by state, for example:
    # if we have L(t) and T(l,Phos~U~P) then the reaction
    #   LigandReceptor: L(t) + T(l) <-> L(t!1).T(l!1) k_lr_bind, k_lr_dis
    # Should get split into
    #   LigandReceptor_forward_U: L(t) + T(l,Phos~U) -> L(t!1).T(l!1,Phos~U) k_lr_bind
    #   LigandReceptor_forward_P: L(t) + T(l,Phos~P) -> L(t!1).T(l!1,Phos~P) k_lr_bind
    #   LigandReceptor_reverse_U: L(t!1).T(l!1,Phos~U) -> L(t) + T(l,Phos~U) k_lr_dis
    #   LigandReceptor_reverse_P: L(t!1).T(l!1,Phos~P) -> L(t) + T(l,Phos~P) k_lr_dis
    _reactions: defaultdict[str, ReactionRule] = defaultdict()
    _species: list[Molecule] = []
    _ligands: list[ComplexReactant] = []

    parameters_dict: dict[str, Parameter] = dict()
    rate_constants: set[str] = set()
    species_init_concentrations: dict[str, Parameter] = dict()
    observables: list[Observable] = []

    def add_molecule_type(self, declaration: str):
        molecule = MoleculeType(declaration)
        if molecule.name in self.molecules:
            raise ValueError(f"Molecule with name '{molecule.name}' already declared.")

        self.molecules[molecule.name] = molecule
        for molecule in self.molecules.values():
            for specie in generate_species(molecule):
                self._species.append(specie)

    def add_molecules(self, declarations: list[str]):
        for declaration in declarations:
            self.add_molecule_type(declaration)


    def add_observable(self, declaration: str):
        observable = Observable(declaration)
        matches1 = [r for r in self._species if r.match_declaration(observable.pattern)]
        matches2 = [r for r in self._ligands if r.match_declaration(
            observable.pattern)]
        self.observables.append(observable)


    def add_observables(self, declarations: list[str]):
        for declaration in declarations:
            self.add_observable(declaration)

    def get_reagents_from_name(self, name: str):
        match = [
            reagent for reagent in self._species if reagent.match_declaration(name)]
        return match

    def add_ligand(self, ligand: ComplexReactant):
        parts: list[list[Molecule]] = []
        for part in ligand.parts:
            aux = [reagent2.copy()
                   for reagent2 in self._species if reagent2.match(part)]
            for reagent2 in aux:
                reagent2.bonds = part.bonds.copy()
            parts.append(aux)                

        for prod in product(*parts):
            ligand_str = '.'.join([str(a) for a in prod])
            ligand2 = ComplexReactant(ligand_str)
            if not any(ligand2.match_ligand(ligand3) for ligand3 in self._ligands):
                self._ligands.append(ligand2)

    def _add_hidden_reaction(self, reaction: ReactionRule):
        reactions: list[ReactionRule] = []
        if isinstance(reaction, BidirectionalReaction):
            r_forward = reaction.get_forward()
            r_reverse = reaction.get_reverse()
            reactions = [r_forward, r_reverse]
        else:
            reactions.append(reaction)

        # make sure that for each reaction added to self._reactions, every reagent has their states descriminated. For example:
        # if we have L(t) and T(l,Phos~U~P) and the user supplied reaction is
        #   LigandReceptor: L(t) + T(l) <-> L(t!1).T(l!1) k_lr_bind, k_lr_dis
        # Should get split into four reactions
        #   LigandReceptor_forward_T_PHOS_U: L(t) + T(l,Phos~U) -> L(t!1).T(l!1,Phos~U) k_lr_bind
        #   LigandReceptor_forward_T_PHOS_P: L(t) + T(l,Phos~P) -> L(t!1).T(l!1,Phos~P) k_lr_bind
        #   LigandReceptor_reverse_T_PHOS_U: L(t!1).T(l!1,Phos~U) -> L(t) + T(l,Phos~U) k_lr_dis
        #   LigandReceptor_reverse_T_PHOS_P: L(t!1).T(l!1,Phos~P) -> L(t) + T(l,Phos~P) k_lr_dis

        # Note that there might be more than one Molecule with their states not explicit,
        # and Molecules may have more than one state not explicit.
        # Note also that if theres a molecule in both sides of the reaction with implicit states,
        # then for each hidden reaction their states should match (like the example above).

        # If every reagent in a reaction has no implicit states, then we can add the reaction to _reactions, directly

        for reaction in reactions:
            # TODO: finish this
            pass
        

    def add_reaction_rule(self, reaction_str: str):
        dict_aux = {key: val.value for key,
                    val in self.parameters_dict.items()}
        
        # create reaction
        try:
            reaction = ReactionRule.from_declaration(reaction_str)
            value = eval(reaction.forward_rate, {"__builtins__": None}, dict_aux)
            reaction.forward_rate_value = value
        except ValueError:
            reaction = BidirectionalReaction.from_declaration(reaction_str)
            value_forward = eval(reaction.forward_rate, {"__builtins__": None}, dict_aux)
            value_reverse = eval(reaction.reverse_rate, {"__builtins__": None}, dict_aux)
            reaction.forward_rate_value = value_forward
            reaction.reverse_rate_value = value_reverse

        # Check that the reaction name is unique
        if reaction.name in self.defined_reactions or reaction.name in self._reactions:
            raise ValueError(f"Reaction name {reaction.name} already exists.")

        # Check that each reagent has been declared as a molecule
        for reagent in (reaction.left_reactants + reaction.right_reactants):

            # check if it is declared as molecule / is in the reagent list
            if isinstance(reagent, Molecule):
                matching_reagents = [
                    reagent2 for reagent2 in self._species if reagent2.match(reagent)]
                if not len(matching_reagents):
                    raise ValueError(f"Reagent {reagent} has not been declared as a molecule.")

            else: # is ligand
                # Check if ligand reagents are valid (formed from declared unbound reagents) 
                for part in reagent.parts:
                    matching_reagents = [
                        reagent2 for reagent2 in self._species if reagent2.match(part)]
                    if not len(matching_reagents):
                        raise ValueError(
                            f"Reagent {reagent} has not been declared as a molecule.")

                self.add_ligand(reagent)

        # check for conservation of mass (this may not necessarily apply)
        # if not reaction.is_mass_conserved():
        #     raise ValueError(f"Conservation of mass is violated for reaction {reaction}")
        
        self.defined_reactions[reaction.name] = reaction
        self._add_hidden_reaction(reaction)
                    

    def add_reaction_rules(self, reactions: list[str]):
        for reaction in reactions:
            self.add_reaction_rule(reaction)

    def add_parameter(self, name: str, declaration: str):
        dict_aux = {key: val.value for key, val in self.parameters_dict.items()}
        value = eval(declaration, {"__builtins__": None}, dict_aux)
        parameter = Parameter(name, declaration, value)
        self.parameters_dict[name] = parameter

    def add_parameters(self, parameters: list[tuple[str, str]]):
        for name, val in parameters:
            self.add_parameter(name, val)

    def add_specie_concentration(self, name: str, value_str: str):
        # see if name matches existing reagent
        match = self.get_reagents_from_name(name)
        if not match:
            raise ValueError(f"Could not match specie with name {name} against a reagent.")
        dict_aux = {key: val.value for key,
                    val in self.parameters_dict.items()}
        value = eval(value_str, {"__builtins__": None}, dict_aux)
        parameter = Parameter(name, value_str, value)
        self.species_init_concentrations[name] = parameter

    def add_species_concentration(self, species: list[tuple[str, str]]):
        for name, declaration in species:
            self.add_specie_concentration(name, declaration)

    def __repr__(self) -> str:
        out = "Molecules:\n"
        for molecule in self.molecules:
            out += f"\t{molecule}\n"

        out += "Reactions:\n"
        for reaction in self.defined_reactions:
            out += f"\t{reaction}\n"

        out += "Parameters:\n"
        for parameter in self.parameters_dict.values():
            out += f"\t{parameter.name} {parameter.declaration}\n"

        out += "Species:\n"
        for specie in self.species_init_concentrations.values():
            out += f"\t{specie.name} {specie.declaration}\n"

        return out

    def simulate(self, total_time: float):

        # compute all the initial concentration values for each possible individual molecule
        # and also all possible complexes
        concentrations: dict[str, int] = dict()
        for reagent in self._species:
            match = [value for name, value in self.species_init_concentrations.items() if reagent.match_declaration(name)]
            if len(match) == 0:
                concentrations[str(reagent)] = 0
            elif len(match) == 1:
                concentrations[str(reagent)] = int(match[0].value)
            else:
                raise ValueError("More than one match.")
        
        for ligand in self._ligands:
            match = [value for name, value in self.species_init_concentrations.items()
                     if ligand.match_declaration(name)]
            if len(match) == 0:
                concentrations[str(ligand)] = 0
            elif len(match) == 1:
                concentrations[str(ligand)] = int(match[0].value)
            else:
                raise ValueError("More than one match.")
            
        time = 0.0
        times: list[float] = [time]
        rates: defaultdict[str, float] = defaultdict()

        observables_map: dict[str, set[str]] = dict()
        observables: dict[str, list[int]] = dict()
        for observable in self.observables:
            label = observable.label
            reagent_str = observable.pattern
            matches1 = [r for r in self._species if r.match_declaration(reagent_str)]
            matches2 = [r for r in self._ligands if r.match_declaration(reagent_str)]
            aux: set[str] = set()
            for match in matches1 + matches2:
                aux.add(str(match))
            observables_map[label] = aux

        for label, reagents in observables_map.items():
            concentration = sum([concentrations[reagent]
                    for reagent in reagents])
            observables[label] = [concentration]

        # apply gillespie algorithm
        while time < total_time:
            # compute the rates for each reaction
            # the rate is equal to the product of the concentrations of the left side reagents and the reaction rate constant
            for name, reaction in self._reactions.items():
                rate_constant = reaction.forward_rate_value
                concentrations_left: list[int] = []
                for reagent in reaction.reactants:
                    if str(reagent) in concentrations:
                        concentrations_left.append(concentrations[str(reagent)])
                    else:
                        raise KeyError(f"{str(reagent)} not in the concentrations dictionary.")
                rate = rate_constant * reduce(lambda x,y: x*y, concentrations_left)
                rates[name] = rate

            total_rate = sum(rates.values())
            if total_rate == 0:
                break

            # Time until the next reaction
            tau = np.random.exponential(1 / total_rate)
            time += tau

            # Determine which reaction occurs
            reaction_choice = np.random.rand() * total_rate
            # select the reaction by computing the comulative distribution function from the rates
            cumulative_rates = list(accumulate(rates.values()))
            chosen_reaction = next(i for i, rate in enumerate(
                cumulative_rates) if rate > reaction_choice)

            # Update concentrations based on the chosen reaction
            reaction_name = list(self._reactions.keys())[chosen_reaction]
            reaction = self._reactions[reaction_name]
            for reagent in reaction.reactants:
                concentrations[str(reagent)] -= 1
            for reagent in reaction.products:
                concentrations[str(reagent)] += 1

            times.append(time)
            for label, reagents in observables_map.items():
                concentration = sum([concentrations[reagent]
                                    for reagent in reagents])
                observables[label].append(concentration)

def test1():
    molecules = ["L(t)",
                 "T(l)"]
    reactions = ["LigandReceptor: L(t) + T(l) <-> L(t!1).T(l!1) k_lr_bind, k_lr_dis"]
    parameters = [("NaV", "6.02e8"),
                  ("L0", "1e4"),
                  ("T0", "7000"),

                  ("k_lr_bind", "8.8e6/NaV"),
                  ("k_lr_dis", "35")]

    species = [("L(t)", "L0"),
               ("T(l)", "T0")]

    model = ReactionModel()
    model.add_molecules(molecules)
    model.add_reaction_rules(reactions)
    model.add_parameters(parameters)
    model.add_species_concentration(species)

    print(model)


def test2():
    molecules = ["L(t)",
                "T(l,Phos~U~P)",
                "CheY(Phos~U~P)",
                "CheZ()"]
    
    reactions = ["LigandReceptor: L(t) + T(l) <-> L(t!1).T(l!1) k_lr_bind, k_lr_dis",
                 "FreeTP: T(l,Phos~U) -> T(l,Phos~P) k_T_phos",
                 "BoundTP: L(t!1).T(l!1,Phos~U) -> L(t!1).T(l!1,Phos~P) k_T_phos*0.2",
                 "YP: T(Phos~P) + CheY(Phos~U) -> T(Phos~U) + CheY(Phos~P) k_Y_phos",
                 "YDep: CheZ() + CheY(Phos~P) -> CheZ() + CheY(Phos~U) k_Y_dephos"]

    parameters = [("NaV", "6.02e8"), 
                  ("L0", "1e4"), 
                  ("T0", "7000"),
                  ("CheY0", "20000"),
                  ("CheZ0", "6000"),

                  ("k_lr_bind", "8.8e6/NaV"), 
                  ("k_lr_dis", "35"),
                  ("k_T_phos", "15"),
                  ("k_Y_phos", "3.8e6/NaV"),
                  ("k_Y_dephos", "8.6e5/NaV")]

    species = [("L(t)", "L0"), 
               ("T(l,Phos~U)", "T0*0.8"),
               ("T(l,Phos~P)", "T0*0.2"),
               ("CheY(Phos~U)", "CheY0*0.5"),
               ("CheY(Phos~P)", "CheY0*0.5"),
               ("CheZ()", "CheZ0")]
    
    observables = [
        "Molecules phosphorylated_CheY CheY(Phos~P)", 
        "Molecules phosphorylated_CheA T(Phos~P)",
        "Molecules bound_ligand L(t!1).T(l!1)"]

    model = ReactionModel()
    model.add_molecules(molecules)
    model.add_parameters(parameters)
    model.add_reaction_rules(reactions)
    model.add_species_concentration(species)
    model.add_observables(observables)
    model.simulate(3.0)

    print(model)


def test3():
    molecules = ["L(t)",
                 "T(l,r,Phos~U~P,Meth~A~B~C)",
                 "CheY(Phos~U~P)",
                 "CheZ()",
                 "CheB(Phos~U~P)",
                 "CheR(t)"]

    reactions = ["LigandReceptor: L(t) + T(l) <-> L(t!1).T(l!1) k_lr_bind, k_lr_dis",

                 # Receptor complex (specifically CheA) autophosphorylation
                 # Rate dependent on methylation and binding states
                 # Also on free vs. bound with ligand
                 "TaUnboundP: T(l,Meth~A,Phos~U) -> T(l,Meth~A,Phos~P) k_TaUnbound_phos",
                 "TbUnboundP: T(l,Meth~B,Phos~U) -> T(l,Meth~B,Phos~P) k_TaUnbound_phos*1.1",
                 "TcUnboundP: T(l,Meth~C,Phos~U) -> T(l,Meth~C,Phos~P) k_TaUnbound_phos*2.8",
                 "TaLigandP: L(t!1).T(l!1,Meth~A,Phos~U) -> L(t!1).T(l!1,Meth~A,Phos~P) 0",
                 "TbLigandP: L(t!1).T(l!1,Meth~B,Phos~U) -> L(t!1).T(l!1,Meth~B,Phos~P) k_TaUnbound_phos*0.8",
                 "TcLigandP: L(t!1).T(l!1,Meth~C,Phos~U) -> L(t!1).T(l!1,Meth~C,Phos~P) k_TaUnbound_phos*1.6",

                 # CheY phosphorylation by T and dephosphorylation by CheZ
                 "YP: T(Phos~P) + CheY(Phos~U) -> T(Phos~U) + CheY(Phos~P) k_Y_phos",
                 "YDep: CheZ() + CheY(Phos~P) -> CheZ() + CheY(Phos~U) k_Y_dephos",

                 # CheR binds to and methylates receptor complex
                 # Rate dependent on methylation states and ligand binding
                 "TRBind: T(r) + CheR(t) <-> T(r!2).CheR(t!2) k_TR_bind, k_TR_dis",
                 "TaRUnboundMeth: T(r!2,l,Meth~A).CheR(t!2) -> T(r,l,Meth~A) + CheR(t) k_TaR_meth",
                 "TbRUnboundMeth: T(r!2,l,Meth~B).CheR(t!2) -> T(r,l,Meth~B) + CheR(t) k_TaR_meth*0.1",
                 "TaRLigandMeth: T(r!2,l!1,Meth~A).L(t!1).CheR(t!2) -> T(r,l!1,Meth~A).L(t!1) + CheR(t) k_TaR_meth*30",
                 "TbRLigandMeth: T(r!2,l!1,Meth~B).L(t!1).CheR(t!2) -> T(r,l!1,Meth~B).L(t!1) + CheR(t) k_TaR_meth*3",

                 # CheB is phosphorylated by receptor complex, and autodephosphorylates
                 "CheBphos: T(Phos~P) + CheB(Phos~U) -> T(Phos~U) + CheB(Phos~P) k_B_phos",
                 "CheBdephos: CheB(Phos~P) -> CheB(Phos~U) k_B_dephos",

                 # CheB demethylates receptor complex
                 # Rate dependent on methyaltion states
                 "TbDemeth: T(Meth~B) + CheB(Phos~P) -> T(Meth~A) + CheB(Phos~P) k_Tb_demeth",
                 "TcDemeth: T(Meth~C) + CheB(Phos~P) -> T(Meth~B) + CheB(Phos~P) k_Tc_demeth",

                 ]

    parameters = [("NaV", "6.02e8"),
                  ("miu", "1e-6"),
                  ("L0", "1e4"),
                  ("T0", "7000"),
                  ("CheY0", "20000"),
                  ("CheZ0", "6000"),
                  ("CheR0", "120"),
                  ("CheB0", "250"),

                  ("k_lr_bind", "8.8e6/NaV"),
                  ("k_lr_dis", "35"),

                  ("k_TaUnbound_phos", "7.5"),

                  ("k_Y_phos", "3.8e6/NaV"),
                  ("k_Y_dephos", "8.6e5/NaV"),

                  ("k_TR_bind", "2e7/NaV"),
                  ("k_TR_dis", "1"),
                  ("k_TaR_meth", "0.08"),

                  ("k_B_phos", "1e5/NaV"),
                  ("k_B_dephos", "0.17"),

                  ("k_Tb_demeth", "5e4/NaV"),
                  ("k_Tc_demeth", "2e4/NaV"),
                  ]

    species = [("L(t)", "L0"),

               ("T(l,r,Meth~A,Phos~U)", "T0*0.84*0.9"),
               ("T(l,r,Meth~B,Phos~U)", "T0*0.15*0.9"),
               ("T(l,r,Meth~C,Phos~U)", "T0*0.01*0.9"),
               ("T(l,r,Meth~A,Phos~P)", "T0*0.84*0.1"),
               ("T(l,r,Meth~B,Phos~P)", "T0*0.15*0.1"),
               ("T(l,r,Meth~C,Phos~P)", "T0*0.01*0.1"),

               ("CheY(Phos~U)", "CheY0*0.71"),
               ("CheY(Phos~P)", "CheY0*0.29"),
               ("CheZ()", "CheZ0"),

               ("CheB(Phos~U)", "CheB0*0.62"),
               ("CheB(Phos~P)", "CheB0*0.38"),

               ("CheR(t)", "CheR0")
               ]

    model = ReactionModel()
    model.add_molecules(molecules)
    model.add_parameters(parameters)
    model.add_species_concentration(species)
    model.add_reaction_rules(reactions)
    model.simulate(3.0)

    print(model)


def test4():
    molecules = ["L(rec)"
                 "R(lig, dim, y1068~0~p, y1173~0~p)"]

    reactions = [
        # Dimer dependent phosphorilation
        "R(dim!+,y1068~0) -> R(dim!+,y1068~p) k_ph1",
        "R(dim!+,y1173~0) -> R(dim!+,y1173~p) k_ph2",

        # Receptor background dephosphorilation
        "R(y1068~p) -> R(y1068~0) k_deph",
        "R(y1173~p) -> R(y1173~0) k_deph",

        # Ligand-dependent Dimerization
        "R(lig,dim) + L(rec) <-> R(lig!0,dim).L(rec!0) k_f, k_r"
        "R(lig!+,dim) + R(lig!+,dim) <-> R(lig!+,dim!0).R(lig!+,dim!0) k_f_d, k_r_d"

        # Ligand association
        "R(lig) + L(rec) -> R(lig!1).L(rec!1) k_f",

        # Ligand dissociation for monomer, singly-ligated dimer and doubly-ligated dimer
        "R(dim,lig!1).L(rec!1) -> R(dim,lig) + L(rec) k_r",
        "R(dim!2,lig).R(dim!2,lig!1).L(rec!1) -> R(dim!2,lig).R(dim!2,lig) + L(rec) k_r/alpha",
        "R(dim!2,lig!+).R(dim!2,lig!1).L(rec!1) -> R(dim!2,lig!+).R(dim!2,lig) + L(rec) k_r/beta",

        # Dimer association
        "R(dim) + R(dim) -> R(dim!1).R(dim!1) k_f_d",

        # Dimer dissociation for unligated/singly-ligated/doubly-ligated dimers
        "R(lig,dim!1).R(lig,dim!1) -> R(lig,dim) + R(lig,dim) k_r_d",
        "R(lig!+,dim!1).R(lig,dim!1) -> R(lig!+,dim) + R(lig,dim) k_r_d/alpha",
        "R(lig!+,dim!1).R(lig!+,dim!1) -> R(lig!+,dim) + R(lig!+,dim) k_r_d/(alpha*beta)",
    ]

    parameters = [("V_ext", "1.6e-9"),  # liters
                  ("N_Avo", "6.022e23"),  # molecule number per mole
                  ("R0", "1e5"),  # molecule number per cell
                  ("L0", "L_conc*V_ext*N_Avo") # molecule number
                  ]

    species = [("(R(lig, dim, y1068~0, y1173~0)", "R0"),
               ("L(rec)", "L0")
               ]

    observables = [
        ("Molecules", "BoundLigand", "L(rec!+)"),
        ("Molecules", "BoundReceptor", "R(lig!+)"),
        ("Species", "Dimer", "R(dim!0).R(dim!0)"),
        ("Species", "UnligatedSpecies", "R(lig,dim), R(lig,dim!0).R(lig,dim!0)"),
        ("Species", "PhosphSpecies", "R(y1068~p), R(y1173~p)"),
    ]

    model = ReactionModel()
    model.add_molecules(molecules)
    model.add_parameters(parameters)
    model.add_species_concentration(species)
    model.add_reaction_rules(reactions)
    model.simulate(3.0)

    print(model)


if __name__ == "__main__":
    test2()
