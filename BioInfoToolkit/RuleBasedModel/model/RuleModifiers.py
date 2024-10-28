
from enum import Enum
from typing import TypedDict

from BioInfoToolkit.RuleBasedModel.model.MoleculeType import MoleculeType
from BioInfoToolkit.RuleBasedModel.utils.parsing_utils import ModifiersDict, MoleculeDict
from BioInfoToolkit.RuleBasedModel.model.Pattern import Pattern


class AllowedModTypes(Enum):
    """Allowed modification types enumeration

    Args:
        Enum (_type_): _description_
    """
    DELETE_MOLECULES = "DeleteMolecules"
    INCLUDE_REACTANTS = "include_reactants"
    INCLUDE_PRODUCTS = "include_products"
    EXCLUDE_REACTANTS = "exclude_reactants"
    EXCLUDE_PRODUCTS = "exclude_products"


class RuleModTypedDict(TypedDict):
    """Rule modifier typed dictionary for parsing

    Args:
        TypedDict (_type_): _description_
    """
    # allowed types: DeleteMolecules, include_[reactants/products],  exclude_[reactants/products]
    type: str
    idx: int
    excluded: list[list[MoleculeDict]]


def mod_type_from_string(type_str: str) -> AllowedModTypes:
    """Return the AllowedModTypes corresponding to the string. 
    If the string does not match raises a value error.

    Args:
        type_str (str): _description_

    Raises:
        ValueError: error if the string is not valid

    Returns:
        AllowedModTypes: _description_
    """
    match type_str:
        case "DeleteMolecules":
            return AllowedModTypes.DELETE_MOLECULES
        case "include_reactants":
            return AllowedModTypes.INCLUDE_REACTANTS
        case "include_products":
            return AllowedModTypes.INCLUDE_PRODUCTS
        case "exclude_reactants":
            return AllowedModTypes.EXCLUDE_REACTANTS
        case "exclude_products":
            return AllowedModTypes.EXCLUDE_PRODUCTS
        case _:
            raise ValueError(f"{type_str} must match the one of the allowed keywords.")


class RuleModifiers:
    delete_molecules: bool

    def __init__(self, mod_dict: ModifiersDict) -> None:
        self.delete_molecules = mod_dict['delete_molecules']
