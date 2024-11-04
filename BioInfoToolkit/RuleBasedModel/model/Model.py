
from BioInfoToolkit.RuleBasedModel.model.ModelBlocks import CompartmentsBlock, InvalidModelBlockError, \
    MoleculeTypesBlock, ObservablesBlock, ParametersBlock, ReactionRulesBlock, SeedSpeciesBlock

class Model:
    molecule_types_block: MoleculeTypesBlock
    parameters_block: ParametersBlock
    reaction_rules_block: ReactionRulesBlock
    observables_block: ObservablesBlock
    compartments_block: CompartmentsBlock
    species_block: SeedSpeciesBlock
    filename: str

    def __init__(self) -> None:
        self.parameters_block = ParametersBlock()
        self.observables_block = ObservablesBlock()
        self.molecule_types_block = MoleculeTypesBlock()
        self.species_block = SeedSpeciesBlock()
        self.reaction_rules_block = ReactionRulesBlock()
        self.compartments_block = CompartmentsBlock()
        self.filename = "model.bngl"

    def validate(self) -> bool:
        molecule_types = self.molecule_types_block.items

        # evaluate parameters
        try:
            self.parameters_block.evaluate_parameters()
        except (TypeError, ValueError) as exc:
            msg = "Parameters block is not valid."
            raise InvalidModelBlockError(msg) from exc
        variables = self.parameters_block.evaluated_params

        # validate compartments
        try:
            self.compartments_block.validate()
        except (TypeError, ValueError) as exc:
            msg = "Compartments block is not valid."
            raise InvalidModelBlockError(msg) from exc

        # if the model is compartmentalized, then we need to validate the definition
        # of observables, species and rules with the compartments
        compartments = self.compartments_block.as_compartments()

        #validate observables
        if not self.observables_block.validate(molecule_types, compartments):
            msg = "Observables block is not valid."
            raise InvalidModelBlockError(msg)

        # validate species
        if not self.species_block.validate_species(molecule_types, compartments):
            msg = "Species block is not valid, due to invalid species."
            raise InvalidModelBlockError(msg)
        if not self.species_block.validate_expressions(variables):
            msg = "Species block is not valid, due to invalid expression."
            raise InvalidModelBlockError(msg)

        # validate rules
        if not self.reaction_rules_block.validate_reactants(molecule_types, compartments):
            msg = "Reaction rules block is not valid, due to invalid pattern."
            raise InvalidModelBlockError(msg)
        if not self.reaction_rules_block.validate_rates(variables):
            msg = "Reaction rules block is not valid, due to invalid rate expression."
            raise InvalidModelBlockError(msg)

        # check if reactions can be decomposed into simple transformations
        if self.is_compartmentalized():
            msg = "Reaction decomposition for compartmentalized model is not implemented."
            raise NotImplementedError(msg)
        try:
            self.reaction_rules_block.decompose_reactions()
        except ValueError as exc:
            msg = ("Reaction rules block is not valid, could not decompose " +
                   "all rules into simple transformations.")
            raise InvalidModelBlockError(msg) from exc

        return True

    def is_compartmentalized(self) -> bool:
        return len(self.compartments_block.items) > 0

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
