from matplotlib import rc
from cycler import cycler

# Axes
rc("axes", linewidth=0.5)
rc("axes", edgecolor="gray")
rc("axes.spines", top=False)
rc("axes.spines", right=False)

# Lines
rc("lines", linewidth=1.25)
rc("lines", markersize=5)

# Ticks
rc("xtick", labelsize=8)
rc("ytick", labelsize=8)

# Legend
rc("legend", frameon=False)
rc("legend", loc="upper right")

# Grid
rc("grid", linestyle="solid")
rc("grid", linewidth=0.5)
rc("grid", alpha=0.5)

# TeX
rc("text.latex", preamble=r"\usepackage{txfonts}")
rc("text", usetex=True)

# Fonts
rc("font", family="sans-serif")
rc("font", serif="times")

# color and line style cycle
solarized_colors = cycler(
    color=["#2aa198", "#cb4b16", "#268bd2", "#859900", "#b58900",
           "#d33682", "#6c71c4"])
# https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=8
dark2 = cycler(
    color=['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e',
           '#e6ab02', '#a6761d', '#666666']
)
d3_catagory10 = cycler(
    color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
)
line_style_cycler = cycler(linestyle=['-', '--', ':', '-.',
                                      '-', '--', ':', '-.',
                                      '-', '--'])

rc("axes", prop_cycle=d3_catagory10 + line_style_cycler)

# Figure size
# onecol_width = 4.2519699737097
# onecol_height = 2.627861962896592
onecol_width = 3.4  # (Phys Rev)
twocol_width = 7 # (Phys Rev)
onecol_height = 3
# margins
left_margin = 0.15  # / onecol_width
right_margin = 0.95  # / onecol_width
bottom_margin = 0.15  # / onecol_height
top_margin = 0.95  # / onecol_height
rc("figure.subplot", bottom=bottom_margin)
rc("figure.subplot", right=right_margin)
rc("figure.subplot", left=left_margin)
rc("figure.subplot", top=top_margin)

def plot_for(journal="PRD", plot_type="onecol", nrows=1):
    if plot_type == "onecol":
        rc("figure", figsize=(onecol_width, onecol_height * nrows))
    else:
        rc("figure", figsize=(twocol_width, onecol_height * nrows))
    if journal=="PRD":
        labelsize = 10
        fontsize = 10
    elif journal=="APJ":
        labelsize = 8
        fontsize = 8
    else:
        labelsize = 9
        fontsize = 9
        
    rc("axes", labelsize=labelsize)
    rc("axes", titlesize=labelsize)
    rc("legend", fontsize=fontsize)

def customize_plot(ax):
    ax.legend()
    ax.grid()
