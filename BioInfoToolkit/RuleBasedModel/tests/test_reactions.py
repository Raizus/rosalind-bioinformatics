import pytest

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.model.ReactionRule import ReactionRule, build_chemical_array_graph, compare_chemical_array_graphs, decompose_reaction, verify_transformations
from BioInfoToolkit.RuleBasedModel.model.ReactionTransformations import apply_transforms


def test_valid_1():
    molecules = {
        "R": MoleculeType("R(lig,ch~open~closed)"),
        "L": MoleculeType("L(rec)")
    }

    react_declaration = "LR: R(lig) + L(rec) -> R(lig!0).L(rec!0) k"
    reaction = ReactionRule.from_declaration(react_declaration, molecules)
    assert isinstance(reaction, ReactionRule)


class TestReaction:
    molecules = {
        "L": MoleculeType("L(t)"),
        "T": MoleculeType("T(l,r,Phos~U~P,Meth~A~B~C)"),
        "CheY": MoleculeType("CheY(Phos~U~P)"),
        "CheZ": MoleculeType("CheZ()"),
        "CheB": MoleculeType("CheB(Phos~U~P)"),
        "CheR": MoleculeType("CheR(t)")
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
        transformations = decompose_reaction(reaction)
        g_equal = verify_transformations(reaction.reactants, transformations,
                                         reaction.products_graph)
        assert g_equal is True


class TestReaction2:
    molecules = {
        "A": MoleculeType("A()"),
        "B": MoleculeType("B()"),
        "C": MoleculeType("C()"),
        "S": MoleculeType("S(x~u~p)"),
    }

    @pytest.mark.parametrize("declaration", [
        "R1: S(x~u) + A() + B() + C() -> S(x~p) + A() + B() + C() k1"
    ])
    def test_valid_1(self, declaration: str):
        reaction = ReactionRule.from_declaration(declaration, self.molecules)
        assert isinstance(reaction, ReactionRule)
