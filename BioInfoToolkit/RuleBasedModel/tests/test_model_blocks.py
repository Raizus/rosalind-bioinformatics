
from BioInfoToolkit.RuleBasedModel.model.ModelBlocks import MoleculeTypesBlock, ObservablesBlock, ParametersBlock, ReactionRulesBlock, SeedSpeciesBlock
from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Observable import Observable
from BioInfoToolkit.RuleBasedModel.model.Parameter import Parameter
from BioInfoToolkit.RuleBasedModel.model.Parsers import parse_parameter
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import ReactionRule
from BioInfoToolkit.RuleBasedModel.model.Species import Species


class TestMoleculeTypesBlock():
    mol_type_str_list = [
        "L(t)",
        "T(l,Phos~U~P)",
        "CheY(Phos~U~P)",
        "CheZ()",
    ]

    def test_valid(self):
        block = MoleculeTypesBlock()
        for mol_type_str in self.mol_type_str_list:
            mol_type = MoleculeType.from_declaration(mol_type_str)
            block.add_molecule_type(mol_type)
        assert isinstance(block, MoleculeTypesBlock)


class TestParametersBlock():
    parameters_list = [
        "NaV 6.02e8",
        "Vcell 1000",
        "Vec   1000*Vcell",
        "d_pm  0.01",
        "Acell 1000",
        "Vpm   Acell*d_pm",
        "lig_conc 1e-8",
        "L0 lig_conc*NaV*Vec",
        "R0 10000",
        "kp1 1e6/NaV",
        "km1 0.01 ",
    ]

    def test_valid(self):
        block = ParametersBlock()
        for param_str in self.parameters_list:
            parsed_param = parse_parameter(param_str)
            parameter = Parameter(parsed_param["name"], parsed_param["expression"])
            block.add_parameter(parameter)
        assert isinstance(block, ParametersBlock)


class TestObservablesBlock():
    molecules: dict[str, MoleculeType] = {
        "A": MoleculeType.from_declaration("A(b,b,c)"),
        "B": MoleculeType.from_declaration("B(a)"),
        "C": MoleculeType.from_declaration("C(a)"),
    }

    observables_list = [
        "Molecules A A(b)",
        "Species B B(a)",
        "Molecules C C(a!?)",
        "Species test_2 A(b!+)",
    ]

    def test_valid(self):
        block = ObservablesBlock()
        for obs_str in self.observables_list:
            observable = Observable.from_declaration(obs_str, self.molecules)
            block.add_observable(observable)
        assert isinstance(block, ObservablesBlock)


class TestSpeciesBlock():
    molecules = {
        "L": MoleculeType.from_declaration("L(t)"),
        "T": MoleculeType.from_declaration("T(l,r,Phos~U~P,Meth~A~B~C)"),
        "CheY": MoleculeType.from_declaration("CheY(Phos~U~P)"),
        "CheZ": MoleculeType.from_declaration("CheZ()"),
        "CheB": MoleculeType.from_declaration("CheB(Phos~U~P)"),
        "CheR": MoleculeType.from_declaration("CheR(t)")
    }

    species_str_list = [
        "L(t) L0",
        "T(l,r,Meth~A,Phos~U) T0*0.84*0.9",
        "T(l,r,Meth~B,Phos~U) T0*0.15*0.9",
        "T(l,r,Meth~C,Phos~U) T0*0.01*0.9",
        "T(l,r,Meth~A,Phos~P) T0*0.84*0.1",
        "T(l,r,Meth~B,Phos~P) T0*0.15*0.1",
        "T(l,r,Meth~C,Phos~P) T0*0.01*0.1",
        "CheY(Phos~U) CheY0*0.71",
        "CheY(Phos~P) CheY0*0.29",
        "CheB(Phos~U) CheB0*0.62",
        "CheB(Phos~P) CheB0*0.38",
    ]

    def test_valid(self):
        block = SeedSpeciesBlock()
        for species_str in self.species_str_list:
            specie = Species.from_declaration(species_str, self.molecules)
            block.add_species(specie)
        assert isinstance(block, SeedSpeciesBlock)


class TestReactionRulesBlock():
    molecules = {
        "L": MoleculeType.from_declaration("L(t)"),
        "T": MoleculeType.from_declaration("T(l,r,Phos~U~P,Meth~A~B~C)"),
        "CheY": MoleculeType.from_declaration("CheY(Phos~U~P)"),
        "CheZ": MoleculeType.from_declaration("CheZ()"),
        "CheB": MoleculeType.from_declaration("CheB(Phos~U~P)"),
        "CheR": MoleculeType.from_declaration("CheR(t)")
    }

    rules_str_list = [
        "LigandReceptor: L(t) + T(l) <-> L(t!1).T(l!1) k_lr_bind, k_lr_dis",

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

    def test_valid(self):
        block = ReactionRulesBlock()
        for rule_str in self.rules_str_list:
            rule = ReactionRule.from_declaration(rule_str, self.molecules)
            block.add_rule(rule)
        assert isinstance(block, ReactionRulesBlock)
