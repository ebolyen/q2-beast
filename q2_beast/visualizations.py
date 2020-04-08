import os

import pandas as pd
import altair as alt

from q2_beast.formats import BEASTPosteriorDirFmt


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
        ).properties(width=800).interactive(bind_y=False)

        hist = alt.Chart(url).mark_bar().encode(
            x=alt.X('count()', title='Frequency'),
            y=alt.Y(param, type='quantitative', bin=alt.Bin(), title=None),
            color='CHAIN:N'
        ).transform_filter(
            alt.datum.state >= selector.burnin
        ).properties(width=200)

        traceplot = alt.hconcat(line, hist).resolve_scale(y='shared')
        traceplots.append(traceplot)

    dash = alt.vconcat(*traceplots)

    dash.save(os.path.join(output_dir, 'index.html'))
    data.to_json(os.path.join(output_dir, url), orient='records')
