"""Rule-based reaction network generation for evolving chemical mechanisms."""

from .generator import GenerationConfig, Generator
from .models import Network, Reaction, Species
from .pathways import find_pathways
from .rules import ReactionRule, load_rules
from .single_event import annotate_single_events, reaction_center

__version__ = "0.6.0"
__all__ = ["Species", "Reaction", "ReactionRule", "Network", "GenerationConfig", "Generator", "load_rules", "find_pathways", "annotate_single_events", "reaction_center"]
