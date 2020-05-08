import os
import shutil
import tempfile
import subprocess
import pkg_resources

import pandas as pd
import altair as alt
from altair import expr, datum

import qiime2

from q2_beast.formats import BEASTPosteriorDirFmt, NexusFormat


def traceplot(output_dir: str, chains: BEASTPosteriorDirFmt,
              params: str = None):
    CONTROL_FMT = chains[0].control.format
    md5sums = {c.control.view(CONTROL_FMT).md5sum() for c in chains}
    if len(md5sums) > 1:
        raise ValueError("Chains do not share a posterior distribution as they"
                         " were generated with different inputs/parameters/"
                         "priors, so they cannot be visualized together.")
    if params is None:
        params = []
    params = list(reversed(params)) + ['likelihood']
    dfs = []
    for idx, chain in enumerate(chains, 1):
        df = chain.log.view(pd.DataFrame)[['state'] + params]
        df['CHAIN'] = 'Chain %d' % idx
        dfs.append(df)
    M = idx

    data = pd.concat(dfs)

    url = 'data.json'
    gen_end = data['state'].iloc[-1]
    gen_step = gen_end - data['state'].iloc[-2]

    slider = alt.binding_range(min=0, max=gen_end, step=gen_step,
                               name='Burn-in: ')
    selector = alt.selection_single(name="BurnIn", fields=['burnin'],
                                    bind=slider, init={'burnin': 0})
    traceplots = []
    for param in params:
        def P(NAME):
            return '_'.join([param, NAME])

        line = alt.Chart(url).mark_line(
            interpolate='step-after',
            opacity=0.8
        ).encode(
            x=alt.X('state:Q', title='Generation'),
            y=alt.Y(param, type='quantitative', scale=alt.Scale(zero=False)),
            color='CHAIN:N'
        ).add_selection(
            selector
        ).transform_filter(
            alt.datum.state >= selector.burnin
        ).properties(
            width=700,
            title='Posterior Trace Plot: ' + param
        ).interactive(bind_y=False)

        hist = alt.Chart(url).mark_bar().encode(
            x=alt.X('count()', title='Frequency'),
            y=alt.Y(param, type='quantitative', bin=alt.Bin(), title=None),
            color='CHAIN:N'
        ).transform_filter(
            alt.datum.state >= selector.burnin
        ).properties(width=200, title='Posterior Hist.')

        table = alt.Chart(url).mark_text().encode(
            y=alt.Y('row_number:O', axis=None),
        ).transform_window(
            row_number='row_number()'
        ).transform_filter(
            alt.datum.state >= selector.burnin
        ).transform_joinaggregate(
            [alt.AggregatedFieldDef(field=param, op='mean',
                                    **{'as': P('OMEAN')})]
        )

        mean = table.properties(title='Mean').encode(
            text=P('OMEAN:N')
        ).transform_aggregate(
            [alt.AggregatedFieldDef(field=param, op='mean',
                                    **{'as': P('OMEAN')})]
        )


        rhat = table.properties(title='PSRF').transform_joinaggregate(
            [alt.AggregatedFieldDef(field=param, op='mean',
                                    **{'as': P('MEAN')}),
             alt.AggregatedFieldDef(field=param, op='variance',
                                    **{'as': P('VAR')}),
             alt.AggregatedFieldDef(op='count', **{'as': 'N'}),
             alt.AggregatedFieldDef(field=P('OMEAN'), op='min',
                                    **{'as': P('OMEAN')})],
            groupby=['CHAIN']
        ).transform_calculate(
            **{P('MSE'): expr.pow(datum[P('MEAN')] - datum[P('OMEAN')], 2)}
        ).transform_aggregate(
            [alt.AggregatedFieldDef(field=P('MSE'), op='sum',
                                    **{'as': P('MSSE')}),
             alt.AggregatedFieldDef(field=P('VAR'), op='sum',
                                    **{'as': P('SVAR')}),
             alt.AggregatedFieldDef(field='N', op='min',
                                    **{'as': 'N'})]
        ).transform_calculate(
            **{P('B'): datum.N / (M - 1) * datum[P('MSSE')],
               P('W'): 1 / M * datum[P("SVAR")]}
        ).transform_calculate(
            **{P('V'): (datum.N - 1) / datum.N * datum[P('W')]
                       + (M + 1) / (M * datum.N) * datum[P('B')]}
        ).transform_calculate(
            **{P('PSRF'): expr.sqrt(datum[P('V')] / datum[P('W')])}
        ).encode(
            text=P('PSRF:N')
        )

        diag = alt.vconcat(mean, rhat)

        traceplot = alt.hconcat(
            diag,
            alt.hconcat(line, hist).resolve_scale(y='shared')
        )
        traceplots.append(traceplot)

    dash = alt.vconcat(*traceplots).configure_view(strokeWidth=0)

    dash.save(os.path.join(output_dir, 'index.html'))
    data.to_json(os.path.join(output_dir, url), orient='records')

def _get_auspice(name):
    path = pkg_resources.resource_filename('q2_beast', 'auspice')
    return os.path.join(path, name)

def auspice(output_dir: str,
            mcc: NexusFormat,
            time: qiime2.NumericMetadataColumn):

    charon = os.path.join(output_dir, 'charon')
    dist = os.path.join(output_dir, 'dist')
    os.mkdir(charon)

    tip_date = time.to_series().max(skipna=True)
    # TODO verify the ids match the data
    with tempfile.TemporaryDirectory(prefix="q2auspice") as tmp:
        out_tree = os.path.join(tmp, 'tree.new')
        out_beast_nd = os.path.join(tmp, 'beast_data.json')
        config = _get_auspice('config.json')

        subprocess.run(['augur', 'import', 'beast', '--mcc', str(mcc),
                        '--output-tree', out_tree, '--output-node-data',
                        out_beast_nd, '--most-recent-tip-date', str(tip_date)],
                       check=True)

        subprocess.run(['augur', 'export', 'v2', '--tree', out_tree,
                        '--node-data', out_beast_nd, '--auspice-config',
                        config, '--output', os.path.join(charon, 'getDataset')
                       ], check=True)

    shutil.copyfile(_get_auspice('index.html'),
                    os.path.join(output_dir, 'index.html'))

    shutil.copyfile(_get_auspice('getAvailable.json'),
                    os.path.join(charon, 'getAvailable'))

    shutil.copytree(_get_auspice('dist'), dist)
