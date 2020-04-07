from qiime2.plugin import (
    Plugin, Properties, MetadataColumn, Numeric, Int, Range, Bool, List)

from q2_types.feature_data import FeatureData, AlignedSequence
from q2_types.tree import Phylogeny

import q2_beast
from q2_beast.methods import (
    site_heterogeneous_hky, merge_chains, maximum_clade_credibility)
from q2_beast.types import Chain, BEAST, MCC
from q2_beast.formats import (
    PosteriorLogFormat, NexusFormat, BEASTControlFileFormat,
    BEASTOpsFileFormat, BEASTPosteriorDirFmt, NexusDirFmt)

plugin = Plugin(
    name='beast',
    website='',
    package='q2_beast',
    version=q2_beast.__version__,
    description='',
    short_description=''
)


plugin.register_formats(
    PosteriorLogFormat, NexusFormat, BEASTControlFileFormat,
    BEASTOpsFileFormat, BEASTPosteriorDirFmt, NexusDirFmt)

plugin.register_semantic_types(Chain, BEAST, MCC)
plugin.register_semantic_type_to_format(
    Chain[BEAST], artifact_format=BEASTPosteriorDirFmt)
plugin.register_semantic_type_to_format(
    Phylogeny[MCC], artifact_format=NexusDirFmt)

NONZERO_INT = Int % Range(1, None)
NONNEGATIVE_INT = Int % Range(0, None)

plugin.methods.register_function(
    function=site_heterogeneous_hky,
    inputs={
        'coding_regions': FeatureData[AlignedSequence % Properties('ORF')],
        'noncoding_regions': FeatureData[AlignedSequence % Properties('NC')]},
    parameters={'time': MetadataColumn[Numeric],
                'time_uncertainty': MetadataColumn[Numeric],
                'n_generations': NONZERO_INT,
                'sample_every': NONZERO_INT,
                'print_every': NONZERO_INT,
                'use_gpu': Bool,
                'n_threads': NONZERO_INT},
    outputs=[('chain', Chain[BEAST])],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name='',
    description='')

plugin.methods.register_function(
    function=merge_chains,
    inputs={'chains': List[Chain[BEAST]]},
    parameters={'burn_in': List[NONNEGATIVE_INT],
                'resample': NONNEGATIVE_INT},
    outputs=[('posterior', Chain[BEAST])],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name='',
    description='')

plugin.methods.register_function(
    function=maximum_clade_credibility,
    inputs={'posterior': Chain[BEAST]},
    parameters={'burn_in': NONNEGATIVE_INT},
    outputs=[('tree', Phylogeny[MCC])],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name='',
    description='')
