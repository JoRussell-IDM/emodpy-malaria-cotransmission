import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform
from idmtools.entities import IAnalyzer

class TimeseriesAnalyzer(IAnalyzer):
    data_group_names = ['group', 'sim_id', 'channel']
    ordered_levels = ['channel', 'group', 'sim_id']
    output_file = 'sim_id_mapping.csv'

    def __init__(self, filenames=[os.path.join('output', 'InsetChart.json')], channels=('Daily EIR',
                                                                                        'PfHRP2 Prevalence',
                                                                                        'Statistical Population'),
                 save_output=True):

        super(TimeseriesAnalyzer, self).__init__(filenames=filenames)
        self.channels = set(channels)
        self.save_output = save_output

    def initialize(self):
        if not os.path.exists(os.path.join(self.working_dir, "output")):
            os.mkdir(os.path.join(self.working_dir, "output"))

    def default_select_fn(self, ts):
        return pd.Series(ts)

    def default_group_fn(self, k, v):
        return k

    def default_plot_fn(self, df, ax):
        grouped = df.groupby(level=['group'], axis=1)
        m = grouped.mean()
        m.plot(ax=ax, legend=False, alpha = 0.2)

    def default_filter_fn(self, md):
        return True

    def filter(self, simulation):
        return self.default_filter_fn(simulation.tags)

    def get_channel_data(self, data_by_channel, selected_channels):
        channel_series = [self.default_select_fn(data_by_channel[channel]["Data"]) for channel in selected_channels]
        return pd.concat(channel_series, axis=1, keys=selected_channels)

    def map(self, data, simulation):
        cdata = data[self.filenames[0]]['Channels']
        selected_channels = self.channels.intersection(cdata.keys()) if self.channels else cdata.keys()
        return self.get_channel_data(cdata, selected_channels)

    def plot_by_channel(self, channels, plot_fn):

        import matplotlib.pyplot as plt

        ncol = int(1 + len(channels) / 4)
        nrow = int(np.ceil(float(len(channels)) / ncol))

        fig, axs = plt.subplots(figsize=(max(6, min(8, 4 * ncol)), min(6, 3 * nrow)), nrows=nrow, ncols=ncol,
                                sharex=True)

        flat_axes = [axs] if ncol * nrow == 1 else axs.flat
        for (channel, ax) in zip(channels, flat_axes):
            ax.set_title(channel)
            plot_fn(channel, ax)

    def reduce(self, all_data):
        output_dir = os.path.join(self.working_dir, "output")
        selected = []
        all_df = pd.DataFrame()
        for sim, data in all_data.items():
            # Enrich the data with info
            sim_id = sim.uid
            aEIR = np.mean(data['Daily EIR'][(25*365):(26*365)])*365
            dict = {'sim_id':[sim_id], 'aEIR':[aEIR]}
            # sim_df = pd.DataFrame(dict)
            # all_df = pd.concat([all_df,sim_df])
            data.group = self.default_group_fn(sim.tags['habitat_multiplier'],sim.uid)
            data.sim_id = sim.uid
            selected.append(data)
        all_df.to_csv(os.path.join(output_dir, self.output_file))

        if len(selected) == 0:
            print("\n No data have been returned... Exiting...")
            return


        # Combining selected data...
        combined = pd.concat(selected, axis=1,
                             keys=[(d.group, d.sim_id) for d in selected],
                             names=self.data_group_names)

        # Re-ordering multi-index levels...
        data = combined.reorder_levels(self.ordered_levels, axis=1).sort_index(axis=1)

        if self.save_output:
            data.to_csv(os.path.join(output_dir, self.output_file))

        def plot_fn(channel, ax):
            self.default_plot_fn(data[channel].dropna(), ax)

        channels = data.columns.levels[0]
        self.plot_by_channel(channels, plot_fn)

        plt.show()
        plt.savefig(os.path.join(output_dir, 'timeseries.png'))


if __name__ == "__main__":
    platform = Platform('Calculon')

    sim_id = '1b2eae64-8120-ec11-9ecd-9440c9bee941'  # comps2 exp_id
    exp_id = 'c66fa3d2-435c-ec11-a9f1-9440c9be2c51'
    # filenames = ['output/ReportSimpleMalariaTransmissionJSON.json']
    filenames = ['output/InsetChart.json']
    analyzers = [TimeseriesAnalyzer(filenames=filenames)]

    manager = AnalyzeManager(platform=platform, ids=[(exp_id, ItemType.EXPERIMENT)], analyzers=analyzers)
    manager.analyze()