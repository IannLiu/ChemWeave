"""Rule-based reaction network generation for evolving chemical mechanisms."""

from .generator import GenerationConfig, Generator
from .models import Network, Reaction, Species
from .pathways import find_pathways
from .rules import ReactionRule, load_rules
from .single_event import annotate_single_events, reaction_center, register_single_event_calculator
from .single_event_calculators import methane_h_abstraction as _methane_h_abstraction
from .path_degeneracy import reaction_path_degeneracy as _reaction_path_degeneracy, path_family_rules, calculate_path_factor

register_single_event_calculator("methane_h_abstraction_v1", _methane_h_abstraction)
register_single_event_calculator("path_degeneracy_v1", _reaction_path_degeneracy)
register_single_event_calculator("mapped_event_v1", _reaction_path_degeneracy)

__version__ = "0.6.0"
__all__ = ["Species", "Reaction", "ReactionRule", "Network", "GenerationConfig", "Generator", "load_rules", "find_pathways", "annotate_single_events", "reaction_center", "register_single_event_calculator", "path_family_rules", "calculate_path_factor"]
