import json
import os
from typing import Any

from idmtools.core.interfaces.iitem import IItem
from idmtools.entities.ianalyzer import IAnalyzer as BaseAnalyzer
import matplotlib as mpl

mpl.use('Agg')


class DownloadAnalyzer(BaseAnalyzer):
    """
       This analyzer is based on the DownloadAnalyzer and allows the download of files based on tags
    """

    def __init__(self, filenames, name='idm'):
        super().__init__(filenames=[filenames])
        print(name)

    def initialize(self):
        if not os.path.exists(os.path.join(self.working_dir, "output")):
            os.mkdir(os.path.join(self.working_dir, "output"))

        # output_dir = r'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\example_with_cotransmission\transmission_report_outputs\genepi_test_suite_EIR_sweep_08112021'

    def map(self, data: Any, item: IItem) -> Any:

        return data[self.filenames[0]]["Channels"]["Statistical Population"]["Data"]

    def reduce(self, all_data: dict) -> Any:
        output_dir = os.path.join(self.working_dir, "output")


if __name__ == "__main__":
    platform = Platform('COMPS2')

    filenames = ['output/InsetChart.json']
    analyzers = [PopulationAnalyzer(working_dir=".")]

    exp_id = '8bb8ae8f-793c-ea11-a2be-f0921c167861'  # comps2 exp_id
    manager = AnalyzeManager(platform=platform, ids=[(exp_id, ItemType.EXPERIMENT)], analyzers=analyzers)
    manager.analyze()