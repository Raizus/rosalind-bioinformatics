import pytest

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern, match_pattern_specie
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import ReactionRule, \
    decompose_reaction, verify_transformations
from BioInfoToolkit.RuleBasedModel.model.ReactionTransformations import apply_transforms
from BioInfoToolkit.RuleBasedModel.model.chemical_array import compare_chemical_array_graphs
from BioInfoToolkit.RuleBasedModel.model.chemical_array import build_chemical_array_graph


def test_valid_1():
    molecules = {
        "R": MoleculeType.from_declaration("R(lig,ch~open~closed)"),
        "L": MoleculeType.from_declaration("L(rec)")
    }

    react_declaration = "LR: R(lig) + L(rec) -> R(lig!0).L(rec!0) k"
    reaction = ReactionRule.from_declaration(react_declaration, molecules)
    assert isinstance(reaction, ReactionRule)


class TestReaction:
    molecules = {
        "L": MoleculeType.from_declaration("L(t)"),
        "T": MoleculeType.from_declaration("T(l,r,Phos~U~P,Meth~A~B~C)"),
        "CheY": MoleculeType.from_declaration("CheY(Phos~U~P)"),
        "CheZ": MoleculeType.from_declaration("CheZ()"),
        "CheB": MoleculeType.from_declaration("CheB(Phos~U~P)"),
        "CheR": MoleculeType.from_declaration("CheR(t)")
    }

    @pytest.mark.parametrize("declaration", [
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
    ])
    def test_valid_1(self, declaration: str):
        reaction = ReactionRule.from_declaration(declaration, self.molecules)
        assert isinstance(reaction, ReactionRule)

    @pytest.mark.parametrize("declaration", [
        "LigandReceptor: L(t!1).T(l!1) -> L(t) + T(l) k_lr_dis",
        "LigandReceptor: L(t) + T(l) -> L(t!1).T(l!1) k_lr_bind",
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
    ])
    def test_decompose(self, declaration: str):
        reaction = ReactionRule.from_declaration(declaration, self.molecules)
        # reaction_graph_node_matchings(reaction)
        # reactants_to_products_node_map(reaction)
        transformations = decompose_reaction(reaction)
        g_equal = verify_transformations(reaction.reactants, transformations,
                                         reaction.products_graph)
        assert g_equal is True


class TestReaction2:
    molecules = {
        "L": MoleculeType.from_declaration("L(t)"),
        "T": MoleculeType.from_declaration("T(l,Phos~U~P)"),
        "CheY": MoleculeType.from_declaration("CheY(Phos~U~P)"),
        "CheZ": MoleculeType.from_declaration("CheZ()")
    }

    @pytest.fixture
    def species(self):
        out = [
            Pattern.from_declaration("L(t)", self.molecules),
            Pattern.from_declaration("T(l,Phos~U)", self.molecules),
            Pattern.from_declaration("T(l,Phos~P)", self.molecules),
            Pattern.from_declaration("CheY(Phos~U)", self.molecules),
            Pattern.from_declaration("CheY(Phos~P)", self.molecules),
            Pattern.from_declaration("CheZ()", self.molecules),
            Pattern.from_declaration("L(t!1).T(l!1,Phos~U)", self.molecules),
            Pattern.from_declaration("L(t!1).T(l!1,Phos~P)", self.molecules),
        ]
        return out

    @pytest.mark.parametrize("declaration", [
        "LigandReceptor: L(t) + T(l) -> L(t!1).T(l!1) k_lr_bind",
        "LigandReceptor_r: L(t!1).T(l!1) -> L(t) + T(l) k_lr_dis",
        "FreeTP: T(l,Phos~U) -> T(l,Phos~P) k_T_phos",
        "BoundTP: L(t!1).T(l!1,Phos~U) -> L(t!1).T(l!1,Phos~P) k_T_phos*0.2",
        "YP: T(Phos~P) + CheY(Phos~U) -> T(Phos~U) + CheY(Phos~P) k_Y_phos",
        "YDep: CheZ() + CheY(Phos~P) -> CheZ() + CheY(Phos~U) k_Y_dephos"
    ])
    def test_valid_1(self, declaration: str):
        reaction = ReactionRule.from_declaration(declaration, self.molecules)
        transformations = decompose_reaction(reaction)
        reactants = reaction.reactants
        products_res = apply_transforms(reactants, transformations)
        g2 = build_chemical_array_graph(products_res)
        equal = compare_chemical_array_graphs(reaction.products_graph, g2)
        assert equal is True


class TestReaction3:
    molecules = {
        "Lig": MoleculeType.from_declaration("Lig(l,l)"),
        "Lyn": MoleculeType.from_declaration("Lyn(U,SH2)"),
        "Syk": MoleculeType.from_declaration("Syk(tSH2,l~Y~pY,a~Y~pY)"),
        "Rec": MoleculeType.from_declaration("Rec(a,b~Y~pY,g~Y~pY)")
    }

    @pytest.mark.parametrize("declaration", [
        "R2: Rec(a) + Lig(l,l!+) <-> Rec(a!2).Lig(l!2,l!+)  kp2, km2",
    ])
    def test_valid_1(self, declaration: str):
        reaction = ReactionRule.from_declaration(declaration, self.molecules)
        transformations = decompose_reaction(reaction)
        reactants = reaction.reactants
        products_res = apply_transforms(reactants, transformations)
        g2 = build_chemical_array_graph(products_res)
        equal = compare_chemical_array_graphs(reaction.products_graph, g2)
        assert equal is True


class TestReaction4:
    molecules = {
        "L": MoleculeType.from_declaration("L(l,l)"),
        "R": MoleculeType.from_declaration("R(r,r)"),
    }

    @pytest.fixture
    def species(self):
        out = [
        Pattern.from_declaration("R(r,r)", self.molecules),
        Pattern.from_declaration("L(l,l)", self.molecules),
        Pattern.from_declaration("L(l,l!1).R(r,r!1)", self.molecules),
        Pattern.from_declaration(
                "L(l,l!1).L(l,l!2).R(r!2,r!1)", self.molecules),
        Pattern.from_declaration("L(l!1,l!2).R(r,r!2).R(r,r!1)", self.molecules),
        Pattern.from_declaration("L(l!1,l!2).L(l,l!3).R(r!1,r!3).R(r,r!2)", self.molecules),
        Pattern.from_declaration(
            "L(l!1,l!2).L(l,l!3).L(l,l!4).R(r!4,r!1).R(r!3,r!2)", self.molecules),
        Pattern.from_declaration(
            "L(l!1,l!2).L(l!3,l!4).R(r!4,r!2).R(r,r!3).R(r,r!1)", self.molecules),
        Pattern.from_declaration(
            "L(l!1,l!2).L(l!3,l!4).L(l,l!5).R(r!1,r!3).R(r!2,r!5).R(r,r!4)", self.molecules),
        ]
        return out

    @pytest.fixture
    def reactions(self):
        out = [
            ReactionRule.from_declaration(
                "R(r) + L(l,l) -> R(r!1).L(l!1,l) kp1", self.molecules),
            ReactionRule.from_declaration(
                "R(r!1).L(l!1,l) -> R(r) + L(l,l) km1", self.molecules),
            ReactionRule.from_declaration(
                "R(r) + L(l,l!+) -> R(r!1).L(l!1,l!+) kp2", self.molecules),
            ReactionRule.from_declaration(
                "R(r!1).L(l!1,l!+) -> R(r) + L(l,l!+) km2", self.molecules)
        ]
        return out

    @pytest.mark.parametrize("sel_r, sel_sp, prod_sp", [
        (0, [0, 1], [2]),
        (0, [2, 1], [3]),
        (1, [2], [0, 1]),
        (2, [0, 2], [4]),
        (2, [2, 2], [5]),
        (0, [4, 1], [5]),
        (0, [5, 1], [6]),
        (1, [3], [1, 2]),
        (1, [5], [1, 4]),
        (2, [0, 3], [5]),
        (2, [0, 5], [7]),
        (2, [2, 3], [6]),
        (2, [2, 5], [8]),
        (2, [4, 2], [7]),
        (2, [4, 3], [8]),
    ])
    def test_valid(self, species, reactions, sel_r: int, sel_sp: list[int], prod_sp: list[int]):
        reaction = reactions[sel_r]
        transformations = decompose_reaction(reaction)
        # reactants = reaction.reactants
        selected_species = [species[i] for i in sel_sp]
        assert len(selected_species) == len(reaction.reactants)
        for sp, react in zip(selected_species, reaction.reactants):
            assert match_pattern_specie(react, sp) is 1

        products_species = apply_transforms(selected_species, transformations)
        assert len(products_species) == len(prod_sp)
        expected_products = [species[i] for i in prod_sp]

