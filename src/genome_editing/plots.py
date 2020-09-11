import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import pypipegraph as ppg
from pathlib import Path
import pandas as pd
import numpy as np
import scipy as sp
import cutadapt
import cutadapt.align



def plot_variants(output_file, sample_to_results, dependencies = []):
    """
    Displays all detected variants along the amplicon for multiple samples.
    SNPs Insertions and deletions are depicted as separate subplots.
    """
    if isinstance(output_file, str):
        outfile = Path(output_file)
    outfile.parent.mkdir(parents = True, exist_ok = True)
    def __dump():
        
        mycolorcycle = cycler(color=["r", "b", "c", "g", "m", "y", "k"])
        figure, axes = plt.subplots(3, figsize=(14, 8), sharex=True)

        for i, ydict in enumerate([ys, non_reference, ys_coverage]):
            axes[i].set_prop_cycle(cy)
            for name in ydict:
                y = ydict[name]
                if i == 2:
                    xi = x[y != 0]
                    yi = y[y != 0]
                else:
                    xi = x
                    yi = y
                axes[i].plot(
                    xi,
                    yi,
                    label=name.replace("_", " "),
                    linestyle=linestyles[i],
                    marker=markerstyles[i],
                )
                to_df["position"].extend(xi)
                to_df["value"].extend(yi)
                to_df["label"].extend([ylabels[i]] * len(xi))
            ylim = axes[i].get_ylim()
            axes[i].plot(
                [
                    pam_positions["R175H_1"] - xoffset,
                    pam_positions["R175H_1"] - xoffset,
                ],
                ylim,
                color="orange",
                label="PAM 1",
            )
            axes[i].plot(
                [
                    pam_positions["R175H_2"] - xoffset,
                    pam_positions["R175H_2"] - xoffset,
                ],
                ylim,
                color="olive",
                label="PAM 2",
            )
            axes[i].set(ylabel=ylabels[i])
        plt.xlim([-75, 75])
        plt.gca().invert_xaxis()
        ticks, labels = plt.xticks()
        newticks = np.sort(
            np.array(
                [int(t) for t in ticks]
                + [7675085 - xoffset, 7675086 - xoffset, 7675087 - xoffset]
            )
        )
        ticklabels = (
            [str(int(-t)) for t in ticks if t < (7675086 - xoffset)]
            + ["", "R175", ""]
            + [str(-t) for t in ticks if t > (7675086 - xoffset)]
        )
        plt.xticks(newticks, ticklabels)
        plt.legend(loc="right")
        plt.xlabel("Position")
        plt.tight_layout()
        figure.savefig(str(outfile1))
        dfout = pd.DataFrame(to_df)




    return ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)

def compare_mutation_distribution_nice(
    sample_to_results, outputfile_prefix, callertype, region, genome, dependencies=[]
    ):
    """
    plot the distribution of mutations per positon for sample and control, then perform a T-test for the mean and a KS-Test for the shape.
    """
    deps = [
        r.call_snps() for sublist in list(sample_to_results.values()) for r in sublist
    ] + dependencies
    prefix = Path(outputfile_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outfile1 = outputfile_prefix + "_variants.png"
    outfile4 = outputfile_prefix + "_variants_sliding.pdf"
    pam_positions = {"R175H_1": 7675100, "R175H_2": 7675103}

    def __plot():
        #first combine caller results
        columns = [
            "sample",
            "ID",
            "chr",
            "pos",
            "genes",
            "ref",
            "alt",
            "reference read count",
            "alternate read count",
            "reference frequency",
            "alternate frequency",
            "coverage",
            "p-value",
            "qual",
            "filter",
            "info",
            "format",
            "frequencies",
        ]
        dfs = {}
        for name, caller_results in sample_to_results.items():
            combine_dfs = []
            for caller_result in caller_results:
                print(caller_result.name)
                
                if hasattr(caller_result, "output_type") and (
                    caller_result.output_type == "varscan"
                ):
                    next_df = variant_calling.parse_varscan(
                        caller_result, genome, type=callertype
                    )
                else:
                    next_df = variant_calling.parse_vcf(
                        caller_result, genome, type=callertype
                    )
                next_df["ID"] = next_df["sample"]
                next_df["sample"] = caller_result.input_sample.name
                next_df = next_df[columns]
                # only use the positions in our reference region
                next_df = next_df[
                    (next_df["chr"] == str(region[0]))
                    & (region[1] <= next_df["pos"])
                    & (next_df["pos"] <= region[2])
                ]
                combine_dfs.append(next_df)
            dfs[name] = pd.concat(combine_dfs)
        #now we have the caller results in a df per sample
        
        # plot_histograms
        x = np.arange(region[1], region[2], 1)
        # offset the x
        xoffset = pam_positions["R175H_1"]
        x = x - xoffset
        codon = np.array([7675085, 7675087], dtype=int)
        codon = codon - xoffset
        ys = {}
        ys_coverage = {}
        non_reference = {}
        for name in sample_to_results:
            df = dfs[name]
            y = np.zeros(len(x), dtype=float)
            cov = np.zeros(len(x), dtype=int)
            non_ref = np.zeros(len(x), dtype=float)
            for pos, df_sub in df.groupby("pos"):
                index = int(pos) - region[1]
                y[index] = df_sub["alternate frequency"].sum()
                coverage = df_sub["coverage"].sum()
                reference = df_sub["reference read count"].sum()
                cov[index] = coverage
                non_ref[index] = (coverage - reference) / coverage
            ys[name] = y
            ys_coverage[name] = cov
            non_reference[name] = non_ref
        ylabels = ["variant frequency", "non reference frequency", "coverage"]
        linestyles = ["", "", "-"]
        markerstyles = [".", ".", ""]
        from cycler import cycler

        plt.rcParams["font.family"] = "Latin Modern Roman"
        #        import matplotlib.font_manager
        #        flist = matplotlib.font_manager.get_fontconfig_fonts()
        #        names = [matplotlib.font_manager.FontProperties(fname=fname).get_name() for fname in flist]
        to_df = {"position": [], "value": [], "label": []}
        cy = cycler(color=["r", "b", "c", "g", "m", "y", "k"])
        figure, axes = plt.subplots(3, figsize=(14, 8), sharex=True)
        for i, ydict in enumerate([ys, non_reference, ys_coverage]):
            axes[i].set_prop_cycle(cy)
            for name in ydict:
                y = ydict[name]
                if i == 2:
                    xi = x[y != 0]
                    yi = y[y != 0]
                else:
                    xi = x
                    yi = y
                axes[i].plot(
                    xi,
                    yi,
                    label=name.replace("_", " "),
                    linestyle=linestyles[i],
                    marker=markerstyles[i],
                )
                to_df["position"].extend(xi)
                to_df["value"].extend(yi)
                to_df["label"].extend([ylabels[i]] * len(xi))
            ylim = axes[i].get_ylim()
            axes[i].plot(
                [
                    pam_positions["R175H_1"] - xoffset,
                    pam_positions["R175H_1"] - xoffset,
                ],
                ylim,
                color="orange",
                label="PAM 1",
            )
            axes[i].plot(
                [
                    pam_positions["R175H_2"] - xoffset,
                    pam_positions["R175H_2"] - xoffset,
                ],
                ylim,
                color="olive",
                label="PAM 2",
            )
            axes[i].set(ylabel=ylabels[i])
        plt.xlim([-75, 75])
        plt.gca().invert_xaxis()
        ticks, labels = plt.xticks()
        newticks = np.sort(
            np.array(
                [int(t) for t in ticks]
                + [7675085 - xoffset, 7675086 - xoffset, 7675087 - xoffset]
            )
        )
        ticklabels = (
            [str(int(-t)) for t in ticks if t < (7675086 - xoffset)]
            + ["", "R175", ""]
            + [str(-t) for t in ticks if t > (7675086 - xoffset)]
        )
        plt.xticks(newticks, ticklabels)
        plt.legend(loc="right")
        plt.xlabel("Position")
        plt.tight_layout()
        figure.savefig(str(outfile1))
        dfout = pd.DataFrame(to_df)
        dfout.to_csv(str(outfile1) + ".tsv", index=False, sep="\t")
        figure, axes = plt.subplots(3, figsize=(14, 8), sharex=True)
        to_df = {"position": [], "value": [], "label": []}
        window_size = 10
        for i, ydict in enumerate([ys, non_reference, ys_coverage]):
            axes[i].set_prop_cycle(cy)
            for name in ydict:
                y = ydict[name]
                if i == 2:
                    xi = x[y != 0]
                    yi = y[y != 0]
                    axes[i].plot(
                        xi,
                        yi,
                        label=name.replace("_", " "),
                        linestyle=linestyles[i],
                        marker=markerstyles[i],
                    )
                    to_df["position"].extend(xi)
                    to_df["value"].extend(yi)
                    to_df["label"].extend([ylabels[i]] * len(xi))
                else:
                    xi = x
                    yi = y
                    sliding_y = []
                    for j in range(len(yi) + 1 - window_size):
                        sliding_y.append(np.mean(yi[j : min(len(yi), j + window_size)]))
                    sliding_y = np.array(sliding_y)
                    axes[i].plot(xi[: len(sliding_y)], sliding_y, label=name)
                    to_df["position"].extend(xi)
                    to_df["value"].extend(yi)
                    to_df["label"].extend([ylabels[i]] * len(xi))
            ylim = axes[i].get_ylim()
            axes[i].plot(
                [
                    pam_positions["R175H_1"] - xoffset,
                    pam_positions["R175H_1"] - xoffset,
                ],
                ylim,
                color="orange",
                label="PAM 1",
            )
            axes[i].plot(
                [
                    pam_positions["R175H_2"] - xoffset,
                    pam_positions["R175H_2"] - xoffset,
                ],
                ylim,
                color="olive",
                label="PAM 2",
            )
            axes[i].set_ylabel(ylabels[i])
        plt.xlim([-75, 75])
        plt.gca().invert_xaxis()
        ticks, labels = plt.xticks()
        newticks = np.sort(
            np.array(
                [int(t) for t in ticks]
                + [7675085 - xoffset, 7675086 - xoffset, 7675087 - xoffset]
            )
        )
        ticklabels = (
            [str(int(-t)) for t in ticks if t < (7675086 - xoffset)]
            + ["", "R175", ""]
            + [str(-t) for t in ticks if t > (7675086 - xoffset)]
        )
        plt.xticks(newticks, ticklabels)
        plt.legend(loc="right")
        plt.xlabel("Position")
        plt.tight_layout()
        figure.savefig(str(outfile4))
        dfout = pd.DataFrame(to_df)
        dfout.to_csv(str(outfile4) + ".tsv", index=False, sep="\t")

    return ppg.MultiFileGeneratingJob(
        [str(outfile1), str(outfile4)], __plot
    ).depends_on(deps)



    