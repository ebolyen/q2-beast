from qiime2.plugin import SemanticType

from q2_types.tree import Phylogeny

Chain = SemanticType('Chain', field_names='posterior')
BEAST = SemanticType('BEAST', variant_of=Chain.field['posterior'])
MCC = SemanticType('MCC', variant_of=Phylogeny.field['type'])
