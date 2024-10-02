
from BioInfoToolkit.RuleBasedModel.model.ModelBlocks import MoleculeTypesBlock, ObservablesBlock, ParametersBlock, ReactionRulesBlock, SeedSpeciesBlock


class Model:
    parameters_block: ParametersBlock
    observables_block: ObservablesBlock
    molecule_types_block: MoleculeTypesBlock
    species_block: SeedSpeciesBlock
    reaction_rules_block: ReactionRulesBlock

    def __init__(self) -> None:
        self.parameters_block = ParametersBlock()
        self.observables_block = ObservablesBlock()
        self.molecule_types_block = MoleculeTypesBlock()
        self.species_block = SeedSpeciesBlock()
        self.reaction_rules_block = ReactionRulesBlock()
