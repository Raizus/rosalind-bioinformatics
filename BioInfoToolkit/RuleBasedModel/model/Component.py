from typing import Any


class MoleculeTypeComponent:
    _name: str
    _states: set[str] | None

    def __init__(self, name: str, states: set[str] | None) -> None:
        self._name = name
        self._states = states

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, MoleculeTypeComponent):
            return False

        return self._name == other._name and self._states == other._states

    def __repr__(self) -> str:
        out = f"{self._name}"
        if self._states:
            out = out + '~' + "~".join(sorted(self._states))
        return out

    @property
    def name(self):
        return self._name

    @property
    def states(self):
        return self._states

    def is_stateless(self) -> bool:
        states = self._states
        out = states is None or len(states) == 0
        return out

    def matches_component(self, other: Any) -> bool:
        if not isinstance(other, MoleculeTypeComponent):
            return False

        if self._name != other.name:
            return False

        if ((self._states is not None)
            and (other.states is not None)
            and (self._states != other.states)):
            return False

        return True

    def copy(self):
        states = self._states.copy() if self._states is not None else None
        new = MoleculeTypeComponent(self._name, states)
        return new


class Component(MoleculeTypeComponent):
    _bond: str = ''
    _state: str = ''

    def __init__(self,
                 name: str,
                 states: set[str] | None,
                 state: str = '',
                 bond: str = ''
                 ) -> None:
        super().__init__(name, states)

        if (self._states is not None) and state and (state not in self._states):
            msg = ( f"Reagent component {name} cannot have the state {state}," +
                    f"it is not in the set of allowed states ({self._states})")
            raise ValueError(msg)

        self._bond = bond
        self._state = state

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, state: str):
        if self._states is None:
            self._state = state
        elif state in self._states:
            self._state = state
        else:
            raise ValueError(
                f"state {state} must be in the set of allowed states ({self._states}).")

    @property
    def bond(self):
        return self._bond

    @bond.setter
    def bond(self, bond: str):
        self._bond = bond

    def is_bonded(self) -> bool:
        return len(self.bond) > 0 and self.bond != '?'

    def is_bond_wildcard(self) -> bool:
        return self.bond in ('?', '+')

    def __repr__(self) -> str:
        out = f"{self._name}"
        if self._state:
            out += f"~{self._state}"
        if self.is_bonded():
            out += f"!{self._bond}"

        return out

    def copy(self):
        states = self._states.copy() if self._states is not None else None
        new = Component(self._name, states,
                        self._state, self._bond)
        return new

    def as_tuple(self):
        return (self._name, self._state, self._bond)

    def matches_component(self, other: Any) -> bool:
        matches = super().matches_component(other)
        if not matches:
            return False



        return True


def sort_components(components: list[Component]) -> list[Component]:
    sorted_comps = sorted(components, key=lambda x: (x.name, x.state, x.bond))
    return sorted_comps


def validate_component(component: Component, component_types: list[MoleculeTypeComponent]) -> bool:
    valid = any(component.matches_component(comp2)
                for comp2 in component_types)
    return valid


def components_gen(component: MoleculeTypeComponent | Component):
    bond = ""
    state = ''
    if isinstance(component, Component):
        bond = component.bond
        state = component.state

    # generate a component for each possible state and maintain the bond
    if component.states and len(state) == 0:
        for state in component.states:
            component2 = Component(component.name, component.states, state, bond)
            yield component2
    else:
        yield Component(component.name, component.states, state, bond)


def components_all_equal(components: list[MoleculeTypeComponent] | list[Component]):
    if len(components) <= 1:
        return True
    return all(components[0].matches_component(comp2) for comp2 in components[1:])
