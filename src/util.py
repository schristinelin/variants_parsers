import itertools
import os

import pandas as pd

from constants import classification_map, mutations_dict


def wrangle_clinvar_txt(df):
    # brand new data container
    df = df[
        ~df["Germline review status"].isin(
            ["no assertion criteria provided", "no classification provided"]
        )
    ]
    # df = df[df['Name'].str.contains('\(p.')] ## do not filter this out, put NA in occurences that are null
    df_new = pd.DataFrame()

    # built under the assumption that data retrieved from clinvar all have the same column names
    df_new["transcript_variant"] = df["Name"].str.split(":").str[1]
    df_new[["transcript_variant", "protein_variant"]] = df_new[
        "transcript_variant"
    ].str.split(" ", expand=True)
    df_new["protein_variant"] = df_new["protein_variant"].fillna("NA")
    df_new["protein_variant"] = df_new["protein_variant"].str.replace(
        r"\(|\)", "", regex=True
    )
    df_new["classification"] = df["Germline classification"]
    df_new["classification"] = (
        df_new["classification"]
        .apply(lambda i: [k for k, v in classification_map.items() if i in v])
        .str[0]
    )
    df_new = df_new[df_new["classification"].isin(["Pathogenic", "Benign"])]
    df_new[["suffix", "mid_string", "str1", "str2"]] = ""
    df_new.loc[(df_new["protein_variant"] != "NA"), "suffix"] = (
        df_new.loc[df_new["protein_variant"] != "NA", "protein_variant"]
        .str.split(".")
        .str[0]
        + "."
    )
    df_new["mid_string"] = df_new.loc[
        df_new["protein_variant"] != "NA", "protein_variant"
    ].str.extract("(\d+)")
    df_new[["str1", "str2"]] = (
        df_new.loc[df_new["protein_variant"] != "NA", "protein_variant"]
        .str.split("\d+", expand=True)
        .iloc[:, 0:2]
    )
    df_new["str1"] = df_new["str1"].str.replace("p.", "")
    df_new = df_new[df_new["protein_variant"] != "p.?"]  ## what's with these?
    df_new.loc[(df_new["protein_variant"] != "NA"), "protein_variant"] = (
        df_new.loc[(df_new["protein_variant"] != "NA"), "suffix"]
        + df_new.loc[(df_new["protein_variant"] != "NA"), "str1"]
        .apply(lambda i: [k for k, v in mutations_dict.items() if i in v])
        .str[0]
        + df_new.loc[(df_new["protein_variant"] != "NA"), "mid_string"]
        + df_new.loc[(df_new["protein_variant"] != "NA"), "str2"]
        .apply(lambda i: [k for k, v in mutations_dict.items() if i in v])
        .str[0]
    )

    df_new = df_new.drop(["suffix", "mid_string", "str1", "str2"], axis=1)
    return df_new


def map_genome_codes(str_val):
    import re

    suffix = str_val.split(".")[0] + "."
    mid_string = re.findall("\d+", str_val)[0]
    print(str_val.split("\d+"))
    str1, str2 = str_val.split("\d+")[0], str_val.split("\d+")[1]
    str1 = str1.replace("p.", "")
    new_str_val = (
        suffix
        + str1.apply(lambda i: [k for k, v in mutations_dict.items() if i in v]).str[0]
        + mid_string
        + str2.apply(lambda i: [k for k, v in mutations_dict.items() if i in v]).str[0]
    )

    return new_str_val


def calc_odds_path(df):
    ## p1
    sample_size = len(df)
    num_pathogenic = len(
        df[df["classification"].str.lower().str.contains("pathogenic")]
    )
    # num_benign = len(df[df["classification"].str.lower().str.contains("benign")])
    print(num_pathogenic)
    pathogenic_prop_p1 = num_pathogenic / sample_size
    # benign_prop_p1 = num_benign / sample_size

    ## p2
    ## FUNC
    func_df = df[df["func.class"] == "FUNC"]
    print(df)
    print(func_df)
    sample_size_1 = len(func_df)
    num_pathogenic_func = len(
        func_df[func_df["classification"].str.lower().str.contains("pathogenic")]
    )
    # num_benign_func = len(
    #    func_df[func_df["classification"].str.lower().str.contains("benign")]
    # )

    pathogenic_prop_func_p2 = num_pathogenic_func / sample_size_1
    # benign_prop_func_p2 = num_benign_func / sample_size_1
    print(pathogenic_prop_func_p2)
    ## LOF
    lof_df = df[df["func.class"] == "LOF"]
    print(lof_df)
    sample_size_2 = len(lof_df)
    num_pathogenic_lof = len(
        lof_df[lof_df["classification"].str.lower().str.contains("pathogenic")]
    )
    # num_benign_lof = len(
    #    lof_df[lof_df["classification"].str.lower().str.contains("benign")]
    # )

    pathogenic_prop_lof_p2 = num_pathogenic_lof / sample_size_2
    # benign_prop_lof_p2 = num_pathogenic_lof / sample_size_2

    ## func oddspath
    func_oddspath = (pathogenic_prop_func_p2 * (1 - pathogenic_prop_p1)) / (
        (1 - pathogenic_prop_func_p2) * pathogenic_prop_p1
    )
    func_oddspath_evidence = oddspath_strength_evidence(func_oddspath)

    ## lof oddspath
    lof_oddspath = (pathogenic_prop_lof_p2 * (1 - pathogenic_prop_p1)) / (
        (1 - pathogenic_prop_lof_p2) * pathogenic_prop_p1
    )
    lof_oddspath_evidence = oddspath_strength_evidence(lof_oddspath)

    row_list = [
        func_oddspath,
        func_oddspath_evidence,
        lof_oddspath,
        lof_oddspath_evidence,
    ]

    return row_list


def oddspath_strength_evidence(oddspath_val):
    if oddspath_val < 0.053:
        evidence_strength = "BS3"
    elif oddspath_val >= 0.053 and oddspath_val < 0.23:
        evidence_strength = "BS3_moderate"
    elif oddspath_val >= 0.23 and oddspath_val < 0.48:
        evidence_strength = "BS3_supporting"
    elif oddspath_val >= 0.48 and oddspath_val <= 2.1:
        evidence_strength = "Indetermine"
    elif oddspath_val > 2.1 and oddspath_val <= 4.3:
        evidence_strength = "PS3_supporting"
    elif oddspath_val > 4.3 and oddspath_val <= 18.7:
        evidence_strength = "PS3_moderate"
    elif oddspath_val > 18.7 and oddspath_val <= 350:
        evidence_strength = "PS3"
    elif oddspath_val > 350:
        evidence_strength = "PS3_very_strong"

    return evidence_strength


def plot_scatter(x, y, color, x_label, ylabel, fpath, title, fname):
    import matplotlib.patches
    import matplotlib.pyplot as plt

    levels, categories = pd.factorize(color)
    colors = [plt.cm.tab10(i) for i in levels]  # using the "tab10" colormap
    handles = [
        matplotlib.patches.Patch(color=plt.cm.tab10(i), label=c)
        for i, c in enumerate(categories)
    ]

    plt.scatter(x, y, c=colors, s=5)
    plt.legend(handles=handles)
    plt.gca().set(title=title, xlabel=x_label, ylabel=ylabel)
    plt.savefig(os.path.join(fpath, str(fname)))
    plt.clf()


def plot_hist_pathogenic(score1, score2, model_name, fpath, fname):
    import matplotlib.pyplot as plt

    plt.hist(score1, alpha=0.5, label="benign")
    plt.hist(score2, alpha=0.5, label="pathogenic")
    plt.legend(loc="upper right")
    plt.gca().set(title="Pathogenecity Histogram of " + model_name)
    plt.savefig(os.path.join(fpath, str(fname)))
    plt.clf()


def train_naive_bayes_model(functional_df, clinvar_df, save_plot_path, save_plot_name):
    from sklearn.metrics import accuracy_score
    from sklearn.model_selection import train_test_split
    from sklearn.naive_bayes import GaussianNB
    from sklearn.preprocessing import OneHotEncoder

    clinvar_df = clinvar_df.drop("transcript_variant", axis=1)
    clinvar_df["protein_variant"] = clinvar_df["protein_variant"].str.replace("p.", "")
    df = functional_df.merge(clinvar_df, on="protein_variant")
    df["classification"] = df["classification"].map({"Benign": 0, "Pathogenic": 1})
    X = df.iloc[:, 1:-1].fillna(0).values  # debatable way of filling NAs
    y = df.iloc[:, -1:].values.ravel()

    print(X)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=1
    )
    x_label, y_label = df.iloc[:, 1:-1].columns[0], df.iloc[:, 1:-1].columns[1]
    ## fit NB model on X and y
    classifier = GaussianNB(priors=[0.5, 0.5])
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)
    # accuracy score
    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy:", accuracy)

    all_predictions = classifier.predict(X)

    df_predict = pd.DataFrame(
        {
            "protein_variant": df.iloc[:, 0].values,
            x_label: df[x_label],
            y_label: df[y_label],
            "predicted_class": all_predictions,
        }
    )
    df_predict["func.class"] = df_predict["predicted_class"].map({0: "FUNC", 1: "LOF"})

    plot_scatter(
        df_predict[x_label],
        df_predict[y_label],
        df_predict["predicted_class"],
        x_label,
        y_label,
        save_plot_path,
        "Gaussian NB Clustering Results",
        save_plot_name,
    )

    return df_predict


def train_kmeans_model(functional_df, save_plot_path, save_plot_name):
    import matplotlib.pyplot as plt
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler

    functional_df = functional_df.drop("protein_variant", axis=1)
    x_label, y_label = functional_df.columns[0], functional_df.columns[1]
    X = functional_df.fillna(0).values
    # again, debatable way of filling NAs
    print(X)

    # normalize data
    standard_scaler = StandardScaler()
    X_standard = standard_scaler.fit_transform(X)

    ## elbow method
    inertias = []
    K = range(1, 10)
    for k in K:
        kmeanModel = KMeans(n_clusters=k).fit(X_standard)
        kmeanModel.fit(X_standard)
        inertias.append(kmeanModel.inertia_)

    print("Elbow Method plot saved to ", save_plot_path)
    plt.plot(K, inertias, "bx-")
    plt.xlabel("Values of K")
    plt.ylabel("Inertia")
    plt.title("The Elbow Method using Inertia")
    plt.savefig(os.path.join(save_plot_path, "KMeans_elbow_plot.png"))
    plt.clf()
    k_val = input("Please review the elbow method plot and enter an optimal k-value:")

    # initialize k means model
    model = KMeans(n_clusters=int(k_val), random_state=42)
    kmeans = model.fit(X_standard)

    plot_scatter(
        X[:, 0],
        X[:, 1],
        kmeans.labels_,
        x_label,
        y_label,
        save_plot_path,
        "K-Means Clustering Results",
        save_plot_name,
    )

    df_predict = pd.DataFrame(
        {
            "protein_variant": functional_df.iloc[:, 0].values,
            x_label: functional_df[x_label],
            y_label: functional_df[y_label],
            "predicted_class": kmeans.labels_,
        }
    )


def gaussian_mixture_model(functional_df, save_plot_path, save_plot_name):
    import matplotlib.pyplot as plt
    from sklearn.mixture import GaussianMixture
    from sklearn.preprocessing import StandardScaler

    # Create a Gaussian Mixture model
    gmm = GaussianMixture(n_components=4, random_state=1)
    X = functional_df.drop("protein_variant", axis=1)
    x_label, y_label = functional_df.columns[0], functional_df.columns[1]
    X = X.fillna(0).values

    # normalize data
    standard_scaler = StandardScaler()
    X_standard = standard_scaler.fit_transform(X)

    gmm.fit(X_standard)
    labels = gmm.predict(X_standard)

    plot_scatter(
        X[:, 0],
        X[:, 1],
        labels,
        x_label,
        y_label,
        save_plot_path,
        "Gaussian Mixture Results",
        save_plot_name,
    )

    functional_df["func.class"] = labels
    functional_df["func.class"] = functional_df["func.class"].map(
        {0: "LOF", 1: "FUNC", 2: "LOF", 3: "LOF"}
    )

    return functional_df


# print(df)
