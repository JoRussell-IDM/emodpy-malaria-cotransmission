# Example DownloadAnalyzer for EMOD Experiment
# In this example, we will demonstrate how to create an DownloadAnalyzer to download simulation output files locally

# First, import some necessary system and idmtools packages.
from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.analysis.download_analyzer import DownloadAnalyzer
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform

if __name__ == '__main__':

    # Set the platform where you want to run your analysis
    # In this case we are running in BELEGOST, but this can be changed to run 'Local'
    with Platform('Calculon') as platform:

        # Arg option for analyzer init are uid, working_dir, data in the method map (aka select_simulation_data),
        # and filenames
        # In this case, we want to provide a filename to analyze
        filenames = ['output/ReportSimpleMalariaTransmissionJSON.json','output/InsetChart.json']
        # Initialize the analyser class with the path of the output files to download
        analyzers = [DownloadAnalyzer(filenames=filenames, output_path='download')]

        # Set the experiment you want to analyze
        experiment_id = 'c66fa3d2-435c-ec11-a9f1-9440c9be2c51'  # comps exp id

        # Specify the id Type, in this case an Experiment
        manager = AnalyzeManager(ids=[(experiment_id, ItemType.EXPERIMENT)],
                                 analyzers=analyzers)
        manager.analyze()