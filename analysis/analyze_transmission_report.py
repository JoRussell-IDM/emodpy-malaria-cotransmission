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
    output_file = 'timeseries.csv'

    def __init__(self, filenames=[os.path.join('output', 'ReportSimpleMalariaTransmissionJSON.json')],
                 save_output=True):

        super(TimeseriesAnalyzer, self).__init__(filenames=filenames)
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
        m.plot(ax=ax, legend=False)

    def default_filter_fn(self, md):
        return True

    def filter(self, simulation):
        return self.default_filter_fn(simulation.tags)

    def get_transmission_data(self, data):
        transmission_df = pd.DataFrame(data)
        return transmission_df

    def map(self, data, simulation):
        tdata = data[self.filenames[0]]['transmissions']
        return self.get_transmission_data(tdata)

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

    def plot_biting_timeseries(self,vector_to_human,human_to_vector,roots):
        fig, axs = plt.subplots(3, 1, figsize=(12, 6), sharex=True, sharey=True)
        vector_to_human.acquireTime.hist(bins=np.arange(25 * 365, 25 * 365 + 100, 1), alpha=1, color='#800080',
                                         ax=axs[0])
        human_to_vector.transmitTime.hist(bins=np.arange(25 * 365, 25 * 365 + 100, 1), alpha=1, color='#f2bf33',
                                          ax=axs[1])
        roots.acquireTime.hist(bins=np.arange(25 * 365, 25 * 365 + 100, 1), alpha=1, color='k', ax=axs[2])
        axs[0].set(title='vector to human')
        axs[1].set(title='human to vector', ylabel='counts')
        axs[2].set(title='new introduced roots', xlabel='simulation time')

        plt.savefig(r'C:\git\emodpy-malaria-cotransmission\analysis\output\biting_timeseries.png')
        plt.show()

    def plot_gametocyte_density_by_progeny(self,vector_to_human):
        transmitters = [val for sublist in vector_to_human['transmitInfectionIds'] for val in sublist]
        acquired = [val for sublist in vector_to_human['acquireInfectionIds'] for val in sublist]

        has_progeny = np.unique([inf for inf in acquired if inf in transmitters])
        no_progeny = np.unique([inf for inf in acquired if inf not in transmitters])

        has_progeny_gams = []
        no_progeny_gams = []
        sum_has_gams = []
        sum_no_gams = []

        # for i, row in vector_to_human.iterrows():
        #     transmit_ids = row.transmitInfectionIds
        #     transmit_gams = row.transmitGametocyteDensities
        #     for j, inf in enumerate(transmit_ids):
        #         if inf in has_progeny:
        #             has_progeny_gams.append(transmit_gams[j])
        #         else:
        #             no_progeny_gams.append(transmit_gams[j]

        for i, row in vector_to_human.iterrows():
            transmit_ids = row.transmitInfectionIds
            transmit_gams = row.transmitGametocyteDensities

            for j, inf in enumerate(transmit_ids):
                if inf in has_progeny:
                    has_progeny_gams.append(transmit_gams[j])
                    # sum_has_gams.append(sum(transmit_gams))
                    # break
                else:
                    no_progeny_gams.append(transmit_gams[j])
                    # sum_no_gams.append(sum(transmit_gams))
                    # break
        from scipy.stats.mstats import gmean
        fig, axs = plt.subplots(2, 1, figsize=(12, 6), sharex=True, sharey=True)
        axs[0].hist(has_progeny_gams, bins=np.logspace(-5, 4, 100))
        axs[0].set(xscale='log', xlabel='gametocytes/uL')
        axs[1].hist(no_progeny_gams, bins=np.logspace(-5, 4, 100))
        axs[1].set(xscale='log', xlabel='gametocytes/uL')

        success_mean = round(gmean(has_progeny_gams), 2)
        success_std = round(np.std(sum_has_gams), 2)
        no_mean = round(gmean(no_progeny_gams), 2)
        no_std = round(np.std(sum_no_gams), 2)

        plt.yscale('log')
        axs[0].set(title=f'successful progeny (gmean:{success_mean})')
        axs[1].set(title=f'no progeny (gmean:{no_mean})', xlabel='gametocyte density')
        plt.savefig(r'C:\git\emodpy-malaria-cotransmission\analysis\output\gametocyte_breakdown_by_progeny_gmean_importation.png')

    def reduce(self, all_data):
        output_dir = os.path.join(self.working_dir, "output")
        selected = []
        for sim, data in all_data.items():
            # Break down the data by transmission types
            # Human_to_Vector transmission events are those where acquireIndividualId = 0 and transmitTime == acquireTime
            human_to_vector = data[(data['acquireIndividualId'] == 0) & (data['transmitTime']==data['acquireTime'])]
            vector_to_human = data[(data['acquireIndividualId'] != 0) & (data['vectorId'] != 0)]
            roots = data[data['vectorId'] == 0]

            transmitters = [val for sublist in vector_to_human['transmitInfectionIds'] for val in sublist]
            acquired = [val for sublist in vector_to_human['acquireInfectionIds'] for val in sublist]
            sampled_by_vectors = [val for sublist in human_to_vector['transmitInfectionIds'] for val in sublist]

            has_progeny = np.unique([inf for inf in acquired if inf in transmitters])
            no_progeny = np.unique([inf for inf in acquired if inf not in transmitters])

            made_it_to_a_vector = np.unique([inf for inf in sampled_by_vectors if inf in no_progeny])

            #add additional resolution on whether the cumulative gam density is dragging some along

            self.plot_biting_timeseries(vector_to_human, human_to_vector, roots)
            self.plot_gametocyte_density_by_progeny(vector_to_human)


            plt.show()

            all_gams = [val for sublist in data['transmitGametocyteDensities'] for val in sublist]
            h2v_gams = [val for sublist in human_to_vector['transmitGametocyteDensities'] for val in sublist]
            v2h_gams = [val for sublist in vector_to_human['transmitGametocyteDensities'] for val in sublist]


            data.group = self.default_group_fn(sim.uid, sim.tags)
            data.sim_id = sim.uid
            selected.append(data)

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

        plt.legend()
        # plt.show()
        plt.savefig(os.path.join(output_dir, 'timeseries.png'))


if __name__ == "__main__":
    platform = Platform('Calculon')

    # sim_id = '1b2eae64-8120-ec11-9ecd-9440c9bee941'
    sim_id = 'b1680973-b420-ec11-9ecd-9440c9bee941'
    filenames = ['output/ReportSimpleMalariaTransmissionJSON.json']
    # filenames = ['output/InsetChart.json']
    # import json
    # with open(f'outputs/{sim_id}/ReportSimpleMalariaTransmissionJSON.json') as json_file:
    #     data = json.load(json_file)


    analyzers = [TimeseriesAnalyzer(filenames=filenames)]

    manager = AnalyzeManager(platform=platform, ids=[(sim_id, ItemType.SIMULATION)], analyzers=analyzers)
    manager.analyze()