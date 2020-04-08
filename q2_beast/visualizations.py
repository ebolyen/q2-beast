import os

import pandas as pd
import altair as alt

from q2_beast.formats import BEASTPosteriorDirFmt


def traceplot(output_dir: str, chains: BEASTPosteriorDirFmt):
    CONTROL_FMT = chains[0].control.format
    md5sums = {c.control.view(CONTROL_FMT).md5sum() for c in chains}
    if len(md5sums) > 1:
        raise ValueError("Chains do not share a posterior distribution as they"
                         " were generated with different inputs/parameters/"
                         "priors, so they cannot be visualized together.")

    log1 = chains[0].log.view(pd.DataFrame)
    url = 'data.json'
    gen_end = log1['state'].iloc[-1]
    gen_step = gen_end - log1['state'].iloc[-2]

    slider = alt.binding_range(min=0, max=gen_end, step=gen_step,
                               name='Burn-in: ')
    selector = alt.selection_single(name="BurnIn", fields=['burnin'],
                                    bind=slider, init={'burnin': 0})
    line = alt.Chart(url).mark_line(interpolate='step').encode(
        x=alt.X('state:Q', title='Generation'),
        y=alt.Y('likelihood:Q', scale=alt.Scale(zero=False))
    ).add_selection(
        selector
    ).transform_filter(
        alt.datum.state >= selector.burnin
    ).properties(width=800).interactive(bind_y=False)

    hist = alt.Chart(url).mark_bar().encode(
        x=alt.X('count()', title='Frequency'),
        y=alt.Y('likelihood:Q', bin=alt.Bin(), title=None)
    ).transform_filter(alt.datum.state >= selector.burnin).properties(
        width=200)

    dash = alt.hconcat(
        line,
        hist
    ).resolve_scale(
        y='shared'
    )

    dash.save(os.path.join(output_dir, 'index.html'))
    log1.to_json(os.path.join(output_dir, url), orient='records')
