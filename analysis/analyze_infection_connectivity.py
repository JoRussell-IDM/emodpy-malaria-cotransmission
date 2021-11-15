

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

if __name__ == "__main__":
    report_cotransmission_contagion = pd.read_csv(r'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\Archive_CoTransmission_Development\HumanToVectorTransmissionInReport\ReportCoTransmissionContagionDeposited\ReportCoTransmissionContagionDeposited.csv')
    sub_df = report_cotransmission_contagion[report_cotransmission_contagion['Time']<9225]

    sub_df['log_gam'] = np.log10(sub_df['GametocyteDensity'])
    sub_df['log_infectiousness'] = np.log10(sub_df['Infectiousness'])
    sub_df['infected'] = sub_df['InfectedVectorIDs'].notnull()
    sub_df['length'] = sub_df['InfectedVectorIDs'].apply(lambda x:len(str(x).split("-")))

    infection_events = sub_df[sub_df['InfectedVectorIDs'].notnull()]

    sns.scatterplot(x='log_gam', y='length', data=sub_df, hue='infected')
    plt.title('Number of vectors infected by infected hosts')
    plt.savefig(
        r'C:\git\emodpy-malaria-cotransmission\analysis\output\Number of vectors infected by infected hosts.png')
    plt.show()
    print('Success')