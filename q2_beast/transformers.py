import pandas as pd

from q2_beast.plugin_setup import plugin
from q2_beast.formats import PosteriorLogFormat


@plugin.register_transformer
def _1(ff: PosteriorLogFormat) -> pd.DataFrame:
    return pd.read_csv(str(ff), sep='\t', skip_blank_lines=True , comment='#')
