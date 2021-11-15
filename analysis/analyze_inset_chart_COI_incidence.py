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

    def __init__(self,
                 filenames=[os.path.join('output', 'InsetChart.json')],
                 channels=('Daily EIR',
                        'New Clinical Cases',
                        'Statistical Population',
                        'True Prevalence'),
                 year = 1,
                 save_output=True,
                 af = 0.5,
                 vars = 24):

        super(TimeseriesAnalyzer, self).__init__(filenames=filenames)
        self.channels = set(channels)
        self.save_output = save_output
        self.year = year
        self.vars = vars
        self.af = af

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

    def get_channel_data(self, data_by_channel, selected_channels,tags):
        channel_series = [self.default_select_fn(data_by_channel[channel]["Data"]) for channel in selected_channels]
        df = pd.concat(channel_series, axis=1, keys=selected_channels)
        year_to_report = self.year
        df_sub = df.iloc[year_to_report*365:(year_to_report+1)*365,:]
        df_sub['annual EIR'] = [x * 365 for x in df_sub['Daily EIR']]
        df_sub['LHM'] = tags['larval_habitat_multiplier']
        return df_sub

    def map(self, data, simulation):
        cdata = data[self.filenames[0]]['Channels']
        selected_channels = self.channels.intersection(cdata.keys()) if self.channels else cdata.keys()
        tags = simulation.tags
        return self.get_channel_data(cdata, selected_channels, tags)

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
        for sim, data in all_data.items():
            # Enrich the data with info


            data.group = self.default_group_fn(sim.uid, sim.tags)
            data.sim_id = sim.uid
            selected.append(data)

        if len(selected) == 0:
            print("\n No data have been returned... Exiting...")
            return

        af = self.af
        vars = self.vars


        # true_COI_filepath = r'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\genomics\synthetic_genomes\eir_sweep_210910\year_allCOI_summaryStats_newRootMoi.txt'
        true_COI_filepath = r'C:\Users\jorussell\Downloads\eir_sweep_210910\year_allCOI_summaryStats_newRootMoi_1025.txt'
        data = pd.read_csv(true_COI_filepath, sep=r'\t', header=0)
        # data_year_one = data[(data['year']== year) & (data['af'] == af) & (data['variants']== vars)]
        data_year_one = data[(data['year'] == self.year)]
        data_year_one.reset_index(inplace=True)
        data_year_one['tag'] = data_year_one.index
        data_year_one['trim'] = [float(str(x)[0:5]) for x in data_year_one['multiplier']]

        df = pd.DataFrame()

        for i, sim_ in enumerate(selected):
            aEIR = np.mean(sim_['annual EIR'])
            CI = np.sum(sim_['New Clinical Cases']) / np.mean(sim_['Statistical Population'])
            sub_df = pd.DataFrame(
                {'index': [i], 'aEIR': aEIR, 'CI': CI, 'LHM': round(float(str(sim_['LHM'].values[0])[0:6]), 3)})
            df = pd.concat([df, sub_df])

        # do inner join on LHM and year = 0 with fields of mean,median,sd
        # fig_name = f'Clinical Incidence v. fraction monogenomic (truth)_year{self.year}'
        # fig_name = f'Clinical Incidence v. fraction monogenomic (af = {af}, variants = {vars}, year = {self.year})'
        fig_name = f"Clinical Incidence v. monogenomic fraction (true)"
        merged = df.merge(data_year_one, left_on='LHM', right_on='trim')
        # merged.to_csv(fr'C:\test\merged_coi_true_{self.year}.csv')
        import seaborn as sns
        # merged = merged[merged['af'] == 'true']
        # merged_subset = merged[(merged['af'] == self.af) & (merged['variants']== str(self.vars))]
        # merged_subset_24_100 = merged[(merged['variants'] == str(24)) | (merged['variants'] == str(100))]
        merged_subset_true = merged[(merged['variants'] == str('true'))]

        sns.lmplot(x='percent_mono', y='CI', data=merged_subset_true, order=2,hue = 'variants',palette='k', scatter_kws={'alpha':0.5,'color':'k'})


        # plt.ylim(3.5, 4.2)
        # plt.xlim(0, 0.3)
        plt.title(fig_name)

        plt.xlabel('fraction monogenomic')
        plt.ylabel('average clinical incidence')
        plt.savefig(
            rf'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\example_with_cotransmission\analyzer_output\{fig_name}.png',bbox_inches='tight')

        plt.savefig(
            rf'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\example_with_cotransmission\analyzer_output\{fig_name}.eps',bbox_inches='tight')

        # plt.show()
        plt.close()
        print('success')




if __name__ == "__main__":
    platform = Platform('Calculon')

    exp_id = 'a09d4076-b235-ec11-9ecd-9440c9bee941'  # comps2 exp_id
    years = [2]
    variants = [24]
    afs = [0.5, 0.1]
    # filenames = ['output/ReportSimpleMalariaTransmissionJSON.json']
    filenames = ['output/InsetChart.json']

    # for year in years:
    #     for v in variants:
            # for af in afs:
            #     analyzers = [TimeseriesAnalyzer(filenames=filenames, year = year, vars = v, af = af)]
    analyzers = [TimeseriesAnalyzer(filenames=filenames, year=2, vars=24)]

    manager = AnalyzeManager(platform=platform, ids=[(exp_id, ItemType.EXPERIMENT)], analyzers=analyzers)
    manager.analyze()