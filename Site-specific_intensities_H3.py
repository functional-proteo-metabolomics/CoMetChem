import pandas as pd

# PICor corrcted MS1 data as .csv
MS1 = "210712_H3_PICor_corrected_only_MS1.csv"

SPECIES = {
    536.8228 : {
        "K18" : "12-13",
        "K23" : "13-12"
    },
    538.3322 : {
        "K18" : "12-*",
        "K23" : "*-12"
    },
    539.3356 : {
        "K18" : "13-*",
        "K23" : "*-13"
    }
}

def get_MS1(time, rep, cond, species):
    """Returns the desired MS1 intensity of given experiment and desired species

    Parameters
    ----------
    time : int/float
        Timepoint of the experiment

    rep : int
        Replicate of the experiment

    cond : str
        Condition/Treatment of the experiment

    species : str
        Species of which the MS1 data should be returned

    Returns
    -------
    float
        MS1 Intensity of given experiment
    """
    return df_MS1[(df_MS1["Time in h"] == time) & (df_MS1.Rep == rep) & (df_MS1.Exp == cond)][species].values[0]

drop_cols = [
    "date",
    "fitted_area_b2_light", "fitted_area_b2_heavy",
    "fitted_area_y7_light", "fitted_area_y7_heavy"
    ]

df = pd.read_csv("Results_H3.csv")
df.dropna(inplace=True)
df.drop(columns=drop_cols, inplace=True)

df_MS1 = pd.read_csv(MS1)

df_536 = df.groupby("precursor_mz").get_group(536.8228)
df_538 = df.groupby("precursor_mz").get_group(538.3322)
df_539 = df.groupby("precursor_mz").get_group(539.3356)

df_536.rename(columns=SPECIES[536.8228], inplace= True)
df_536.drop(columns="precursor_mz", inplace=True)

df_538.rename(columns=SPECIES[538.3322], inplace= True)
df_538.drop(columns="precursor_mz", inplace=True)

df_539.rename(columns=SPECIES[539.3356], inplace= True)
df_539.drop(columns="precursor_mz", inplace=True)


df_site_specific = df_536.merge(df_538, on=["Time in h", "Exp", "Rep"])
df_site_specific = df_site_specific.merge(df_539, on=["Time in h", "Exp", "Rep"])

df_site_specific["Exp"] = df_site_specific.apply(lambda row: "MS275" if row["Exp"] == "MS" else row["Exp"], axis=1)
df_site_specific = df_site_specific[df_site_specific.Exp.isin(["C", "SAHA", "MS275"])]



df_site_specific["MS1 2C13"] = df_site_specific.apply(lambda row: get_MS1(
    row["Time in h"],
    row["Rep"],
    row["Exp"],
    "2C13"), axis = 1)


df_site_specific["MS1 2C13 3H02"] = df_site_specific.apply(lambda row: get_MS1(
    row["Time in h"],
    row["Rep"],
    row["Exp"],
    "2C13 3H02"), axis = 1)

df_site_specific["MS1 4C13 3H02"] = df_site_specific.apply(lambda row: get_MS1(
    row["Time in h"],
    row["Rep"],
    row["Exp"],
    "4C13 3H02"), axis = 1)

df_site_specific["*-*"] = df_site_specific.apply(lambda row: get_MS1(
    row["Time in h"],
    row["Rep"],
    row["Exp"],
    "4C13 6H02"), axis = 1)

df_site_specific["12-12"] = df_site_specific.apply(lambda row: get_MS1(
    row["Time in h"],
    row["Rep"],
    row["Exp"],
    "No label"), axis = 1)

df_site_specific["13-13"] = df_site_specific.apply(lambda row: get_MS1(
    row["Time in h"],
    row["Rep"],
    row["Exp"],
    "4C13"), axis = 1)

df_site_specific["12-13"] = df_site_specific["12-13"] * df_site_specific["MS1 2C13"] / 100
df_site_specific["13-12"] = df_site_specific["13-12"] * df_site_specific["MS1 2C13"] / 100

df_site_specific["12-*"] = df_site_specific["12-*"] * df_site_specific["MS1 2C13 3H02"] / 100
df_site_specific["*-12"] = df_site_specific["*-12"] * df_site_specific["MS1 2C13 3H02"] / 100

df_site_specific["13-*"] = df_site_specific["13-*"] * df_site_specific["MS1 4C13 3H02"] / 100
df_site_specific["*-13"] = df_site_specific["*-13"] * df_site_specific["MS1 4C13 3H02"] / 100

df_site_specific.drop(columns=["MS1 2C13", "MS1 2C13 3H02", "MS1 4C13 3H02"]).to_csv("H3_site_specific.csv", index=False)