2/1/2022
Jon Russell

Documentation for running epidemiological layer of EMOD-GenEpi layered model

Directory structure:
/analysis: analyzer scripts for download and plotting of inset chart and transmission report channels

/download: Eradication binaries and schema used to run model in emodpy (passed from Svetlana or DanB in place of ideal behavior using get_model_files, not yet functional as I understand it)

/inputs: demographics files read in for EMOD nodes

/outputs: a dir that gets overwritten with file outputs from run_sweep_and_download or similar scripts that inherit from download analyzer


key files:
run_sweep_and_download.py: commissions, runs sims based on configuration specified within body of script, downloads to outputs dir using in line download analyzer

Steps:
Download and install emodpy-malaria: https://github.com/InstituteforDiseaseModeling/emodpy-malaria

Setup appropriate virtualenv in Pycharm (or similar IDE) using emodpy-malaria

-------------
Commissioning
-------------

Change parameters in run_sweep_and_download.py for particular experiment setup 
	- Ln 36: seasonality = dict to decribe different seasonal splines as basis for tuning seasonal larval habitat abundnace
	- Ln 93: adds outbreak to establish initial prevalence in sim
	- Ln 99: adds recurring outbreaks (local importation) to sustain transmission and prevent low transmission from eliminating
	- Ln 135: generic treatment seekking behavior based on parasite densities triggering pyrogenic threshold
	- Ln 137: defining length of sim in years
	- Ln 170: Function that does rescaling of larval habitat spline to achieve range of EIR
	- Ln 186: Adding local migration for importation from external (travel) node
	- Ln 191: Scaling the importation as reference above in Ln 99 for sustaining transmission
	- Ln 197: For sweeping over run_number, i.e. new random seeds for DTK replicates
	- Ln 219: the main run function for executing a sim (as referenced in if name == main)

	- Ln 229: Submission details (i.e setting priority and number of cores on Calculon/Belegost)
	- Ln 233: Set number of random seeds
	- Ln 234: set the first day for reports to start on (i.e. starting Transmission Report on day 25*365 in a 30 year sim gets you last 5 years of transmission events)
	- Ln 236: an array of values for sweeping human population size
	- Ln 237: the importation scalar determining demographic coverage of the importation used in Ln 99 adding outbreaks to combat elimination, this could be a sweep of values in the list
	- Ln 238: habitat_multipliers is the sweep space used to determine a range of EIR by scaling the larval habitat spline
	- Ln 239: a selection of particular habitat scaling values selected for manuscript
	- Ln 240: where to specify the larval habitat spline shape/seasonality from above dict
	- Ln 242: where to name the experiment

	- Ln 298: where the call to transmission report gets generated from
	- Ln 304-309: Where the sweeps that define the experiment get generated, using calls to appropritae functions and invoking args defined above
	- Ln 332-347: Download Analyzer for downloading the result of the sim, alternatively can use ln 359-368 (commented out) to call download analyzer on particular experiment, or Download Analyzer in /analysis
	- Ln 371: Calls the main run function!

After successful run and download from COMPS, move transmission report output from specified /outputs dir to Dropbox for use as input by rest of GenEpi pipeline (parse_dtk)

--------------
Analyzers
--------------

analyze_inset_chart.py for plotting epi channels containined in InsetChart (population transmission characteristics and demographic information)
analyze_inset_chart_COI_incidence.py for mapping incidence to COI metric
analyze_inset_chart_EIR_sim_id_mapping.py for making a csv output that maps sim_ids to EIR values
analyze_transmission_report.py a starting point for analyzing transmission characteristics from transmission report (also see Edward's blog post)
download_analyzer.py for downloading outputs like InsetCharts or transmission reports from single sims or experiments
analyze_patient_report.py a starting point for looking at MalariaPatientReport output if available for investigating individual properties like gam densities

--------------
Inputs
--------------
For information about setup of demographics files see EMOD documentation or other malaria team workflow examples:
https://docs.idmod.org/projects/emod-malaria/en/latest/index.html
