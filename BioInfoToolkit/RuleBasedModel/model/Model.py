
from BioInfoToolkit.RuleBasedModel.model.ModelBlocks import CompartmentsBlock, \
    MoleculeTypesBlock, ObservablesBlock, ParametersBlock, ReactionRulesBlock, SeedSpeciesBlock


class Model:
    molecule_types_block: MoleculeTypesBlock
    parameters_block: ParametersBlock
    reaction_rules_block: ReactionRulesBlock
    observables_block: ObservablesBlock
    compartments_block: CompartmentsBlock
    species_block: SeedSpeciesBlock

    def __init__(self) -> None:
        self.parameters_block = ParametersBlock()
        self.observables_block = ObservablesBlock()
        self.molecule_types_block = MoleculeTypesBlock()
        self.species_block = SeedSpeciesBlock()
        self.reaction_rules_block = ReactionRulesBlock()
        self.compartments_block = CompartmentsBlock()

    def as_string(self) -> str:
        out = ''
        molecule_types_block_str = self.molecule_types_block.gen_string()
        parameters_block_str = self.parameters_block.gen_string()
        compartments_block_str = self.compartments_block.gen_string()
        species_block_str = self.species_block.gen_string()
        reactions_block_str = self.reaction_rules_block.gen_string()
        observables_block_str = self.observables_block.gen_string()

        strings = [molecule_types_block_str,
                   parameters_block_str,
                   compartments_block_str,
                   species_block_str,
                   reactions_block_str,
                   observables_block_str]

        for string in strings:
            if string:
                out += string

        return out
