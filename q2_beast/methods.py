import pkg_resources
import subprocess

import jinja2
import pandas as pd

import qiime2

from q2_beast.formats import (BEASTPosteriorDirFmt, NexusFormat,
                              PosteriorLogFormat)


def _get_template(name):
    path = pkg_resources.resource_filename('q2_beast',
                                           'xml-templates')
    loader = jinja2.FileSystemLoader(searchpath=path)
    env = jinja2.Environment(loader=loader)
    return env.get_template(name)


def gtr_single_partition(
        alignment: qiime2.Metadata,
        time: qiime2.NumericMetadataColumn,
        n_generations: int,
        sample_every: int,
        time_uncertainty: qiime2.NumericMetadataColumn = None,
        base_freq: str = "estimated",
        site_gamma: int = 4,
        site_invariant: bool = True,
        clock: str = 'ucln',
        coalescent_model: str = 'skygrid',
        skygrid_intervals: int = None,
        skygrid_duration: float = None,
        print_every: int = None,
        use_gpu: bool = False,
        n_threads: int = 1) -> BEASTPosteriorDirFmt:

    if coalescent_model == 'skygrid':
        if skygrid_duration is None or skygrid_intervals is None:
            raise ValueError("skygrid not parameterized (TODO: better error)")

    # Parallelization options
    beast_call = ['beast']
    if use_gpu:
        if n_threads != 1:
            raise ValueError
        beast_call += ['-beagle_GPU', '-beagle_cuda', '-beagle_instances', '1']
    else:
        beast_call += ['-beagle_CPU', '-beagle_SSE',
                       '-beagle_instances', str(n_threads)]

    # Set up directory format where BEAST will write everything
    result = BEASTPosteriorDirFmt()
    control_file = str(result.control.path_maker())

    ops_file = str(result.ops.path_maker().relative_to(result.path))
    log_file = str(result.log.path_maker().relative_to(result.path))
    trees_file = str(result.trees.path_maker().relative_to(result.path))

    # Setup up samples for templating into control file
    seq_series = alignment.get_column('Sequence').to_series()
    time_series = time.to_series()

    if time_uncertainty is not None:
        uncertainty_series = time_uncertainty.to_series()
    else:
        uncertainty_series = time_series.copy()
        uncertainty_series[...] = None

    samples_df = pd.concat([seq_series, time_series, uncertainty_series],
                           axis='columns', join='inner')
    samples_df.index.name = 'id'
    samples_df.columns = ['seq', 'time', 'time_uncertainty']
    samples_df = samples_df.replace({pd.np.nan: None})
    samples = list(samples_df.itertuples(index=True))

    # Default print behavior
    if print_every is None:
        print_every = sample_every

    # Generate control file for BEAST
    template_kwargs = dict(trees_file=trees_file, ops_file=ops_file,
                           log_file=log_file, sample_every=sample_every,
                           print_every=print_every,
                           n_generations=n_generations, time_unit='years',
                           samples=samples, base_freq=base_freq,
                           site_gamma=site_gamma,
                           site_invariant=site_invariant, clock=clock,
                           coalescent_model=coalescent_model,
                           skygrid_duration=skygrid_duration,
                           skygrid_intervals=skygrid_intervals)

    template = _get_template("gtr_single_partition.xml")
    template.stream(**template_kwargs).dump(control_file)

    beast_call += [str(control_file)]

    # Execute
    subprocess.run(beast_call, check=True, cwd=result.path)

    return result


def site_heterogeneous_hky(
        coding_regions: qiime2.Metadata,
        noncoding_regions: qiime2.Metadata,
        time: qiime2.NumericMetadataColumn,
        n_generations: int,
        sample_every: int,
        print_every: int = None,
        time_uncertainty: qiime2.NumericMetadataColumn = None,
        use_gpu: bool = False,
        n_threads: int = 1) -> BEASTPosteriorDirFmt:

    # Parallelization options
    beast_call = ['beast']
    if use_gpu:
        if n_threads != 1:
            raise ValueError
        beast_call += ['-beagle_GPU', '-beagle_cuda', '-beagle_instances', '1']
    else:
        beast_call += ['-beagle_CPU', '-beagle_SSE',
                       '-beagle_instances', str(n_threads)]

    # Set up directory format where BEAST will write everything
    result = BEASTPosteriorDirFmt()
    control_file = str(result.control.path_maker())

    ops_file = str(result.ops.path_maker().relative_to(result.path))
    log_file = str(result.log.path_maker().relative_to(result.path))
    trees_file = str(result.trees.path_maker().relative_to(result.path))

    # Setup up samples for templating into control file
    orf_series = coding_regions.get_column('Sequence').to_series()
    nc_series = noncoding_regions.get_column('Sequence').to_series()
    time_series = time.to_series()
    uncertainty_series = time_uncertainty.to_series()

    samples_df = pd.concat([orf_series, nc_series, time_series,
                            uncertainty_series], axis='columns', join='inner')
    samples_df.index.name = 'id'
    samples_df.columns = ['seq_orf', 'seq_nc', 'time', 'time_uncertainty']
    samples_df = samples_df.replace({pd.np.nan: None})
    samples = list(samples_df.itertuples(index=True))

    # Default print behavior
    if print_every is None:
        print_every = sample_every

    # Generate control file for BEAST
    template_kwargs = dict(trees_file=trees_file, ops_file=ops_file,
                           log_file=log_file, sample_every=sample_every,
                           print_every=print_every,
                           n_generations=n_generations, time_unit='years',
                           samples=samples)
    template = _get_template("orf_and_nc.xml")
    template.stream(**template_kwargs).dump(control_file)

    beast_call += [str(control_file)]

    # Execute
    subprocess.run(beast_call, check=True, cwd=result.path)

    return result


def _log_combiner(files, out, burn_in, is_tree, resample=None):
    combiner_call = ['logcombiner', '-burnin', str(burn_in)]
    if is_tree:
        combiner_call += ['-trees']
    if resample is not None:
        combiner_call += ['-resample', str(resample)]
    combiner_call += list(map(str, files))
    combiner_call += [str(out)]
    subprocess.run(combiner_call, check=True)


def merge_chains(chains: BEASTPosteriorDirFmt, burn_in: int,
                 resample: int = None) -> BEASTPosteriorDirFmt:
    if len(burn_in) > 1 and len(burn_in) != len(chains):
        raise ValueError("burn_in")

    CONTROL_FMT = chains[0].control.format
    md5sums = {c.control.view(CONTROL_FMT).md5sum() for c in chains}
    if len(md5sums) > 1:
        raise ValueError("Chains do not share a posterior distribution as they"
                         " were generated with different inputs/parameters/"
                         "priors, so they cannot be merged.")

    if len(burn_in) > 1:
        logs_to_merge = [PosteriorLogFormat() for _ in chains]
        trees_to_merge = [NexusFormat() for _ in chains]
        # remove the burn_in first as the CLI doesn't allow these to differ
        for chain, single_burn_in, out_trees, out_log in zip(
                chains, burn_in, trees_to_merge, logs_to_merge):
            _log_combiner([chain.log.view(chain.log.format)],
                          out=out_log, burn_in=single_burn_in, is_tree=False)
            _log_combiner([chain.trees.view(chain.trees.format)],
                          out=out_trees, burn_in=single_burn_in, is_tree=True)
        burn_in = 0  # disable global burn-in
    else:
        logs_to_merge = [c.log.view(c.log.format) for c in chains]
        trees_to_merge = [c.trees.view(c.trees.format) for c in chains]
        burn_in = burn_in[0]

    result = BEASTPosteriorDirFmt()
    result.control.write_data(chains[0].control.view(CONTROL_FMT),
                              view_type=CONTROL_FMT)
    with result.ops.path_maker().open('w') as fh:
        fh.write('')  # intentionally empty file

    _log_combiner(logs_to_merge, out=result.log.path_maker(), burn_in=burn_in,
                  is_tree=False, resample=resample)
    _log_combiner(trees_to_merge, out=result.trees.path_maker(),
                  burn_in=burn_in, is_tree=True, resample=resample)

    return result


def maximum_clade_credibility(posterior: BEASTPosteriorDirFmt,
                              burn_in: int = 0) -> NexusFormat:
    result = NexusFormat()

    trees = posterior.trees.view(posterior.trees.format)
    annotator_call = ['treeannotator', '-burnin', str(burn_in),
                      str(trees), str(result)]
    subprocess.run(annotator_call, check=True)

    return result
