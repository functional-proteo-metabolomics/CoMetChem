import pandas as pd
import os

# Path where the folders containing the data are stored
PATH = "D:Data\\20210609_results\\results"

REPLICATES = {
        "191207":1,
        "191214":2,
        "191220":3
    }

TIMEPOINTS = {
    "T0" : 0,
    "T05": 0.5,
    "T1" : 1,
    "T15": 1.5,
    "T2" : 2,
    "T4" : 4,
    "T6" : 6,
    "T8" : 8,
    "T24": 24
}

def get_fitted_area(df, fragment_ion, low = True):
    """Return the fitted area of a given fragment ion.
    Set low to false if the fitted area of the heavy fragment is desired.

    Parameters
    ----------
    df : Dataframe
        Dataframe containing the fitted areas of the fragment ions from one precursor m/z.
        Must contain the columns "fragment_ion", "fragment_mz" and "fitted_area"

    fragment_ion : str
        Desired fragment ion (e.g. "y7")

    low : bool
        True if the fitted area of the light fragment ion is desired.
        False for the heavy fragment ion.

    Returns
    -------
    float
        Fitted area for the given fragment ion
    """
    fitted_areas = df[df.fragment_ion == fragment_ion].sort_values(by=["fragment_mz"]).fitted_area.values
    if low:
        return fitted_areas[0]
    else:
        return fitted_areas[1]

def get_percentage(df, fragment_ion, weight):
    """Return the relative abundance of a given fragment ion.

    Parameters
    ----------
    df : Dataframe
        Dataframe containing the information of the fitted areas of a given fragment ion

    fragment_ion : str
        Desired fragment ion (e.g. "y7")

    weight : str
        Choose weather the light or heavy fragment ion should be selected

    Returns
    -------
    float
        Relative abundance (percentage) of the chosen fragment ion
    """
    return (df[f"fitted_area_{fragment_ion}_{weight}"] / (df[f"fitted_area_{fragment_ion}_light"]+df[f"fitted_area_{fragment_ion}_heavy"])) * 100

cols = [
    "precursor_mz",
    "Time in h",
    "Exp",
    "Rep",
    "fitted_area_b2_light",
    "fitted_area_b2_heavy",
    "fitted_area_y7_light",
    "fitted_area_y7_heavy",   
]

output = pd.DataFrame(columns=cols)

# name of folders:
# replicate_timepoint_condition_scientist_m/z_Retentiontime-start__Retentiontime-stop_scannumber
# e.g. 191207_T0_C_inclusion_AvP_538.33220_1284.00__1344.00_8061

for folder in os.listdir(PATH):
    # to skip unrelated measurements (folder also included H4 measurements)
    if float(folder.split("_")[5]) > 700:
        continue

    df_raw = pd.read_csv(PATH + "/" + folder + "/fragments.csv")
    df_raw = df_raw.fillna(0)

    row = {
        "precursor_mz": df_raw.precursor_mz.values[0],
        "date": folder.split("_")[0],
        "Time in h": folder.split("_")[1],
        "Exp": folder.split("_")[2],        
        "fitted_area_b2_light": get_fitted_area(df_raw, "b2"),
        "fitted_area_b2_heavy": get_fitted_area(df_raw, "b2", low=False),
        "fitted_area_y7_light": get_fitted_area(df_raw, "y7"),
        "fitted_area_y7_heavy": get_fitted_area(df_raw, "y7", low=False),    
    }
    
    output = output.append(row, ignore_index=True)

output["Rep"] = output.date.map(REPLICATES)
output["Time in h"] = output["Time in h"].map(TIMEPOINTS)

# Abundance calculation of the different acetylation sites
output["K18"] = (get_percentage(output, "b2", "light") + get_percentage(output, "y7", "heavy")) / 2
output["K23"] = (get_percentage(output, "b2", "heavy") + get_percentage(output, "y7", "light")) / 2

output.to_csv("Results_H3.csv", index=False)