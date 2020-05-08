import importlib

from qiime2.plugin import (
    Plugin, MetadataColumn, Numeric, Int, Range, Bool, List, Str, Choices,
    Float)

from q2_types.feature_data import FeatureData, AlignedSequence
from q2_types.tree import Phylogeny

import q2_beast
from q2_beast.methods import (
    site_heterogeneous_hky, merge_chains, maximum_clade_credibility,
    gtr_single_partition)
from q2_beast.visualizations import traceplot, auspice
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


importlib.import_module('q2_beast.transformers')

NONZERO_INT = Int % Range(1, None)
NONNEGATIVE_INT = Int % Range(0, None)

plugin.methods.register_function(
    function=gtr_single_partition,
    inputs={
        'alignment': FeatureData[AlignedSequence]},
    parameters={'time': MetadataColumn[Numeric],
                'n_generations': NONZERO_INT,
                'sample_every': NONZERO_INT,
                'time_uncertainty': MetadataColumn[Numeric],
                'base_freq': Str % Choices("estimated", "empirical"),
                'site_gamma': Int % Range(0, 10, inclusive_end=True),
                'site_invariant': Bool,
                'clock': Str % Choices("ucln", "strict"),
                'coalescent_model': Str % Choices("skygrid", "constant",
                                                  "exponential"),
                'skygrid_intervals': NONZERO_INT,
                'skygrid_duration': Float % Range(0, None,
                                                  inclusive_start=False),
                'print_every': NONZERO_INT,
                'use_gpu': Bool,
                'n_threads': NONZERO_INT},
    outputs=[('chain', Chain[BEAST])],
    input_descriptions={
        'alignment': 'The alignment to construct a tree with.',
    },
    parameter_descriptions={
        'time': 'The decimal date for when that sequence was collected.',
        'time_uncertainty': 'Uncertainty in the collection time,'
                            ' this should be in decimal years.',
        'n_generations': 'The number of generations (or iterations) to run the'
                         ' MCMC procedure for. Higher values are more likely'
                         ' to result in samples from the posterior'
                         ' distribution. Typical values are on the order of'
                         ' tens of millions of generations.',
        'sample_every': 'How many generations should occur between samples'
                        ' which will form the chain. This is a thinning '
                        ' parameter, and can be used to reduce autocorrelation'
                        ' increasing your effective sample size.',
        'base_freq': '',
        'site_gamma': '',
        'site_invariant': '',
        'clock': '',
        'coalescent_model': '',
        'skygrid_intervals': '',
        'skygrid_duration': '',
        'print_every': 'How many generations should occur before printing to'
                       ' stdout. This is a cosmetic feature, and by default'
                       ' will match `sample_every`.',
        'use_gpu': 'Whether to perform MCMC on a CUDA enabled GPU.',
        'n_threads': 'The number of threads to use, TODO: this is not quite'
                     ' accurate, as some extra math happens with partitions'
    },
    output_descriptions={
        'chain': 'An output chain of (ideally) the posterior distribution for'
                 ' the phylogenetic analysis. Multiple chains should be'
                 ' analyzed to ensure that each chain has converged on the'
                 ' posterior distribution.'
    },
    name='',
    description='')

plugin.methods.register_function(
    function=site_heterogeneous_hky,
    inputs={
        'coding_regions': FeatureData[AlignedSequence],
        'noncoding_regions': FeatureData[AlignedSequence]},
    parameters={'time': MetadataColumn[Numeric],
                'time_uncertainty': MetadataColumn[Numeric],
                'n_generations': NONZERO_INT,
                'sample_every': NONZERO_INT,
                'print_every': NONZERO_INT,
                'use_gpu': Bool,
                'n_threads': NONZERO_INT},
    outputs=[('chain', Chain[BEAST])],
    input_descriptions={
        'coding_regions': 'An alignment of concatenated open reading frames.',
        'noncoding_regions': 'An alignment of concatenated non-coding regions.'
    },
    parameter_descriptions={
        'time': 'The decimal date for when that sequence was collected.',
        'time_uncertainty': 'Uncertainty in the collection time,'
                            ' this should be in decimal years.',
        'n_generations': 'The number of generations (or iterations) to run the'
                         ' MCMC procedure for. Higher values are more likely'
                         ' to result in samples from the posterior'
                         ' distribution. Typical values are on the order of'
                         ' tens of millions of generations.',
        'sample_every': 'How many generations should occur between samples'
                        ' which will form the chain. This is a thinning '
                        ' parameter, and can be used to reduce autocorrelation'
                        ' increasing your effective sample size.',
        'print_every': 'How many generations should occur before printing to'
                       ' stdout. This is a cosmetic feature, and by default'
                       ' will match `sample_every`.',
        'use_gpu': 'Whether to perform MCMC on a CUDA enabled GPU.',
        'n_threads': 'The number of threads to use, TODO: this is not quite'
                     ' accurate, as some extra math happens with partitions'
    },
    output_descriptions={
        'chain': 'An output chain of (ideally) the posterior distribution for'
                 ' the phylogenetic analysis. Multiple chains should be'
                 ' analyzed to ensure that each chain has converged on the'
                 ' posterior distribution.'
    },
    name='',
    description='')

plugin.methods.register_function(
    function=merge_chains,
    inputs={'chains': List[Chain[BEAST]]},
    parameters={'burn_in': List[NONNEGATIVE_INT],
                'resample': NONNEGATIVE_INT},
    outputs=[('posterior', Chain[BEAST])],
    input_descriptions={
        'chains': 'A list of BEAST chains to merge together.'
    },
    parameter_descriptions={
        'burn_in': 'The number of generations (not samples!) to treat as the'
                   ' warmup period for the MCMC procedure. The samples'
                   ' collected from generations prior to this `burn_in` will'
                   ' be discarded as they are not representative of the'
                   ' posterior distribution. If a single value is given, then'
                   ' all chains will be given the same burn-in.',
        'resample': 'Will preform additional thinning on each chain before'
                    ' merging. This value is in generations (not samples!)'
                    ' and must be an even multiple of the original sampling'
                    ' rate.'  # why can't BEAST just use iter and thin?
    },
    output_descriptions={
        'posterior': 'A merged chain of posterior samples.'},
    name='Merge multiple posterior chains, remove burn-in, and thin.',
    description='Merge multiple posterior chains, remove burn-in, and thin.')

plugin.methods.register_function(
    function=maximum_clade_credibility,
    inputs={'posterior': Chain[BEAST]},
    parameters={'burn_in': NONNEGATIVE_INT},
    outputs=[('tree', Phylogeny[MCC])],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name='Create a Maximum Clade Credibility tree from BEAST.',
    description='Calculate the Maximum Clade Credibility tree from a BEAST'
                ' posterior distribution. Ensure that the chain used has'
                ' properly converged before running, or this will fail to'
                ' produce meaningful results.')

plugin.visualizers.register_function(
    function=traceplot,
    inputs={'chains': List[Chain[BEAST]]},
    parameters={'params': List[Str]},
    input_descriptions={},
    parameter_descriptions={
        'params': 'Additional parameter traces to plot. By default only the'
                  ' log-likelihood is plotted.'},
    name='Create traceplots of BEAST chains.',
    description=''
)


plugin.visualizers.register_function(
    function=auspice,
    inputs={'mcc': Phylogeny[MCC]},
    parameters={'time': MetadataColumn[Numeric]},
    input_descriptions={},
    parameter_descriptions={},
    name='',
    description=''
)


def not_real(output_dir: str, nope: int = None):
    pass


plugin.visualizers.register_function(
    function=not_real,
    inputs={},
    parameters={'nope': Int},
    input_descriptions={},
    parameter_descriptions={},
    name='',
    description='')
