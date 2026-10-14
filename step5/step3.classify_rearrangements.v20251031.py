import pandas as pd
import argparse
import matplotlib.pyplot as plt

# LIS function weighted by block length
def longest_increasing_subsequence_weighted_indices(values, weights, indices):
    n = len(values)
    if n == 0:
        return [], 0
    dp = [w for w in weights]
    prev = [-1]*n

    for i in range(n):
        for j in range(i):
            if values[j] < values[i] and dp[j] + weights[i] > dp[i]:
                dp[i] = dp[j] + weights[i]
                prev[i] = j

    max_idx = max(range(n), key=lambda i: dp[i])
    max_weight = dp[max_idx]

    lis_indices = []
    cur = max_idx
    while cur != -1:
        lis_indices.append(indices[cur])
        cur = prev[cur]
    lis_indices.reverse()
    return lis_indices, max_weight

# Classification and plotting function
def classify_and_plot(df, output_prefix):
    log_file = f"{output_prefix}.log"
    open(log_file, "w").close()  # clear log first

    # Initialize category column
    df["category"] = "NA"
#    df["note"] = ""

    for (sp1, sp2), sub in df.groupby(["sp1", "sp2"]):
        sub_sorted = sub.sort_values("start1").reset_index()
        start1_order = sub_sorted["start1"].tolist()
        start2_order = sub_sorted["start2"].tolist()
        len1_order = sub_sorted["len1"].tolist()
        indices = sub_sorted["index"].tolist()  # old index in df
        dirs = sub_sorted["dir"].tolist()
#        start1 = sub_sorted["start1"].tolist()
#        sub_indices = list(range(len(start1)))  # new index in the sub

        # LIS forward
        forward_idx = [i for i,d in enumerate(dirs) if d==1]
        print(forward_idx)
        print([start1_order[i] for i in forward_idx])
        print([start2_order[i] for i in forward_idx])
        lis_fwd, weight_fwd = longest_increasing_subsequence_weighted_indices(
            [start2_order[i] for i in forward_idx],
            [len1_order[i] for i in forward_idx],
            [i for i in forward_idx]
        )

        # LIS reverse
        reverse_idx = [i for i,d in enumerate(dirs) if d==-1]
        reverse_idx_sorted = sorted(reverse_idx, key=lambda i: start1_order[i])
        
        start2_rev_flipped = [-start2_order[i] for i in reverse_idx_sorted]
        len1_rev_flipped = [len1_order[i] for i in reverse_idx_sorted]
        print(reverse_idx_sorted)
        print([start1_order[i] for i in reverse_idx_sorted])
        print([start2_order[i] for i in reverse_idx_sorted])

        lis_rev, weight_rev = longest_increasing_subsequence_weighted_indices(
            start2_rev_flipped,
            len1_rev_flipped,
            [i for i in reverse_idx_sorted]
        )

        # Determine LIS orientation
        if weight_fwd >= weight_rev:
            orientation = "forward"
            lis_blocks = lis_fwd
            print(lis_blocks)
            print([start1_order[i] for i in lis_blocks])
            # Convert to global df indices:
            lis_global = sub_sorted.loc[lis_blocks, "index"].tolist()
            # Now safely update the original df:
            df.loc[lis_global, "category"] = "normal"  # use the global indices

#            df.loc[lis_blocks, "category"] = "normal"  # mark LIS immediately
#            idx = sub.index
#            sub_flipped = df.loc[idx]  # unflipped
            
        else:
            orientation = "reversed"
            lis_blocks = lis_rev
            # flip all blocks in this chr pair
            idx = sub.index
            df.loc[idx, "start2"] = -sub["start2"].values
            df.loc[idx, "end2"] = -sub["end2"].values
            df.loc[idx, "dir"] = -sub["dir"]

            # Ensure start2 < end2 for consistency
            start2_flipped = df.loc[idx, "start2"].copy()
            end2_flipped = df.loc[idx, "end2"].copy()

            # flip the start and end too so start is still smaller than end
            df.loc[idx, "start2"] = start2_flipped.where(start2_flipped < end2_flipped, end2_flipped)
            df.loc[idx, "end2"]   = end2_flipped.where(end2_flipped > start2_flipped, start2_flipped)

            # Convert to global df indices:
            lis_global = sub_sorted.loc[lis_blocks, "index"].tolist()

            # Now safely update the original df:
            df.loc[lis_global, "category"] = "normal"  # use the global indices
#            df.loc[lis_blocks, "category"] = "normal"
            sub_flipped = df.loc[idx]  # get flipped coordinates


        # change the non-normal blocks to trans or trans-invert based on the dir
        # If category == "NA" and dir == -1 → intra_trans_invert
        df.loc[(df["category"] == "NA") & (df["dir"] == -1), "category"] = "intra_trans_invert"

        # If category == "NA" and dir == 1 → intra_trans
        df.loc[(df["category"] == "NA") & (df["dir"] == 1), "category"] = "intra_trans"
        
        # Log info
        with open(log_file, "a") as log:
            lis_start2 = df.loc[lis_blocks, "start2"].tolist()
            lis_len1 = df.loc[lis_blocks, "len1"].tolist()
            log.write(f"{sp1}\t{sp2}\tLIS_total_weight={sum(lis_len1)}\tOrientation={orientation}\n")
            log.write(f"LIS_start2={lis_start2}\nLIS_len1={lis_len1}\n")

        # Plot
        plt.figure(figsize=(8,6))
        color_map = {
            "normal": "red",
            "intra_trans": "green",
            "intra_invert": "orange", # currently can't detect invert only
            "intra_trans_invert": "orange",
            "NA": "grey"
        }
        print(df)
        print(sub_sorted)
        print("sub_sorted index:", sub_sorted.index[:5].tolist())
        print("sub_sorted['index']:", sub_sorted["index"][:5].tolist())

#        for i, r in sub_sorted.iterrows():
#        for i, r in df.iterrows():
#        for i, r in sub_sorted.iterrows():
        sub_current = df[(df["sp1"] == sp1) & (df["sp2"] == sp2)]  # <- reselect from df
        for _, r in sub_current.iterrows():
            x0, x1 = r["start1"], r["end1"]
            y0, y1 = r["start2"], r["end2"]
            if r["dir"] == -1:
                y0, y1 = y1, y0
            cat_color = color_map[r["category"]]
            plt.plot([x0/1000000,x1/1000000], [y0/1000000,y1/1000000], color=cat_color, linewidth=2)

        plt.xlabel(f"{sp1} coordinates (Mb)")
        plt.ylabel(f"{sp2} coordinates (Mb)")
        plt.title(f"Synteny blocks: {sp1} vs {sp2}")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_{sp1}_{sp2}_synteny.pdf")
        plt.close()

# Main
def main():
    parser = argparse.ArgumentParser(description="Longest Increasing Subsequence-based synteny rearrangement classification and plotting\n"
                                                 "This pipeline classifies syntenic blocks into:\n"
                                                 " Normal: the block belongs to the longest chain of blocks between two chromosome pairs\n"
                                                 " Translocation: the block does not belong to the longest chain\n"
                                                 " Translocation and invertion: the block does not belong to the longest chain and inverted between two chromosomes\n"
                                                 "Updates:\n"
                                                 " v20251030 - Use LIS to orient chromosomes\n"
                                                 " v20251029 - Group by sp1 and sp2 chr pairs and output log\n",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help="Input synteny file")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output prefix for log and plots")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep=r"\s+", header=None,
                     names=["sp1","start1","end1","sp2","start2","end2","len1","len2","dir"])

    classify_and_plot(df, args.output_prefix)

    # Save classified df
    df.to_csv(f"{args.output_prefix}_classified.tsv", sep="\t", index=False)
    print("Classification complete. Counts:")
    print(df["category"].value_counts())

if __name__ == "__main__":
    main()

