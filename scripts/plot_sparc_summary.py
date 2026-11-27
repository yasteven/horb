"""Display the checked SPARC comparisons; no fitting or sample selection here."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BG, PANEL, TEXT, MUTED, GRID = "#0b1320", "#101d2e", "#edf3ff", "#9dafc5", "#28384b"
COLORS = {"squared": "#54dece", "psi": "#bb9aff", "nfw": "#ffb35c", "burkert": "#75a8ff",
          "piso": "#a3cf87", "einasto": "#f582ad", "vms": "#ce9275", "mond": "#7ac9dd",
          "qi": "#f0cf6c", "baryons": "#a7b3c3"}


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "text.color": TEXT,
        "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.edgecolor": GRID, "figure.facecolor": BG, "axes.facecolor": PANEL,
        "savefig.facecolor": BG, "pdf.fonttype": 42})


def save(fig, name):
    (ROOT/"images").mkdir(exist_ok=True)
    (ROOT/"plots").mkdir(exist_ok=True)
    fig.savefig(ROOT/"images"/f"{name}.png", dpi=160)
    fig.savefig(ROOT/"plots"/f"{name}.pdf")
    plt.close(fig)


def plot_rms(summary):
    from make_sparc_summary import LABELS
    style()
    fig = plt.figure(figsize=(13.6, 9.5))
    fig.text(.055, .947, "HORB  /  SPARC ROTATION CURVES", color=COLORS["squared"], weight="bold", size=12)
    fig.text(.055, .878, "Average error across SPARC 175", size=29, weight="bold")
    fig.text(.055, .828, "Mean of each galaxy’s RMS  ·  fixed baryons  ·  lower is better", color=MUTED, size=13)
    ax = fig.add_axes([.275, .235, .555, .535])
    records = sorted(summary, key=lambda r:r["mean_rms_kms"])
    values = [r["mean_rms_kms"] for r in records]
    y = np.arange(len(records))
    ax.barh(y, values, color=[COLORS[r["model"]] for r in records], height=.52)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlim(0, max(values)*1.15)
    ax.set_xlabel("Mean RMS [km/s]", labelpad=11)
    ax.xaxis.grid(True, color=GRID, alpha=.65)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    # Same-galaxy mean provides a direct comparison for every model.
    ax.scatter([r["common_mean_rms_kms"] for r in records], y, marker="|", s=125,
               color=TEXT, linewidths=1.6, zorder=3)
    for i, row in enumerate(records):
        model = row["model"]
        label = LABELS[model]
        detail = ("up to 4 signed weights" if model in ("squared", "psi") else
                  "fixed prediction" if model in ("qi", "mond", "baryons") else
                  "1 fitted amplitude; fixed scale rule" if model == "vms" else
                  "3 fitted parameters" if model == "einasto" else "2 fitted parameters")
        ax.text(-.035, i-.08, label, transform=ax.get_yaxis_transform(), ha="right", va="center",
                color=COLORS[model], size=13, weight="bold" if model in ("squared", "psi") else "normal")
        ax.text(-.035, i+.23, detail, transform=ax.get_yaxis_transform(), ha="right", va="center", color=MUTED, size=9)
        ax.text(row["mean_rms_kms"]+.65, i, f'{row["mean_rms_kms"]:.2f}', va="center", size=13, weight="bold")
        ax.text(1.16, i, str(row["galaxies"])+("†" if row["galaxies"]<175 else ""),
                transform=ax.get_yaxis_transform(), ha="center", va="center", color=MUTED, size=12)
    fig.text(.92, .79, "GALAXIES", size=9, color=MUTED, ha="center", weight="bold")
    fig.text(.055, .148, "Equal weight per galaxy; all 3,391 observations retained for the fitted models.", size=11, color=MUTED)
    fig.text(.055, .111, "† MOND, QI and baryons: 174 valid galaxies. UGC01281 has outward baryonic force at two radii; no clipping.", size=10, color=MUTED)
    fig.text(.055, .075, "White ticks: mean on the same 174 galaxies for every model. Full values and fit settings are saved with the figure.", size=10, color=MUTED)
    fig.text(.055, .039, "In-sample comparison. HORB also searches orbital choices and scales; model flexibility differs.", size=10, color=MUTED)
    save(fig, "sparc_model_rms")


def plot_activity(rows, stats):
    style()
    fig = plt.figure(figsize=(14.5, 11.4))
    gold = "#ffc76a"
    fig.text(.055, .956, "HORB  /  RADIO ACTIVITY AND ORBITAL WEIGHTS", color=COLORS["squared"], weight="bold", size=12)
    fig.text(.055, .904, "Jets and orbital concentration", weight="bold", size=29)
    fig.text(.055, .867, "21 optically classified AGN  ·  7 with resolved jets/lobes  ·  both HORB density families", color=MUTED, size=12)
    handles = [Line2D([],[],marker="*",color=gold,ls="",ms=12,label="Resolved jet/lobes"),
               Line2D([],[],marker="o",color=TEXT,ls="",ms=5,label="Other optical AGN: detected core"),
               Line2D([],[],marker="<",color=MUTED,ls="",ms=6,label="Core luminosity upper limit")]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(.05,.829), frameon=False, ncol=3, fontsize=10, columnspacing=2)
    offsets = {
        "squared": {"NGC2841":(-45,9),"UGC11914":(9,-16),"NGC4051":(-51,-14),"NGC5985":(8,8),
                    "NGC5005":(8,-15),"UGC03546":(9,3),"UGC02487":(-70,11)},
        "psi": {"NGC2841":(-52,-28),"UGC11914":(-58,-14),"NGC4051":(8,8),"NGC5985":(8,-13),
                "NGC5005":(8,8),"UGC03546":(8,-13),"UGC02487":(-70,11)},
    }
    for col, model in enumerate(("squared", "psi")):
        left = .08+col*.48
        ax = fig.add_axes([left, .464, .395, .301])
        selected = [r for r in rows if r["density_model"]==model]
        agn = [r for r in selected if r["optical_class"] in ("Seyfert", "LINER")]
        for row in agn:
            x, y = row["log_Lnu_W_Hz"], 100*row["concentration"]
            jet = row["group"]=="resolved_jets"
            if row["radio_core_limit"]=="upper":
                ax.annotate("", xy=(x-.27,y), xytext=(x,y), arrowprops={"arrowstyle":"->", "color":MUTED, "lw":1.1}, zorder=2)
                ax.plot(x,y,marker="|",color=MUTED,ms=5)
            else:
                ax.scatter(x,y,c=gold if jet else COLORS[model],s=155 if jet else 40,
                           marker="*" if jet else "o",edgecolors=BG,linewidths=.5,zorder=4)
            if jet:
                ax.annotate(row["name"], (x,y),
                    xytext=offsets[model][row["galaxy"]], textcoords="offset points", color=gold, size=8.5,
                    arrowprops={"arrowstyle":"-", "color":gold,"lw":.5}, zorder=5)
        ax.set(xlim=(17.7,24.15),ylim=(24,108),xticks=np.arange(18,25),yticks=[25,50,75,100])
        ax.set_xlabel(r"Radio core luminosity: log$_{10}$ $L_\nu$ [W Hz$^{-1}$]",size=10,labelpad=10)
        if not col:
            ax.set_ylabel("Largest component’s absolute weight [%]",size=10,labelpad=10)
        ax.grid(color=GRID,alpha=.55)
        ax.set_axisbelow(True)
        ax.set_title("Squared density  |ψ|²" if model=="squared" else "Signed wavefunction lobes  ψ",
                     color=COLORS[model],size=15,weight="bold",pad=16)
        s = stats[model]
        fig.text(left,.388, f'Detected AGN cores: Spearman ρ = {s["spearman_rho_detections_only"]:+.2f}, '
                 f'p = {s["spearman_p_two_sided"]:.2f} (n = {s["detections"]})',size=10,color=MUTED)
        bx = fig.add_axes([left,.185,.395,.144])
        groups = ("resolved_jets","other_optical_AGN","HII_comparison")
        for i, group in enumerate(groups):
            members = sorted([r for r in selected if r["group"]==group],key=lambda r:r["galaxy"])
            values = [100*r["concentration"] for r in members]
            jitter = np.random.default_rng(42+i).uniform(-.19,.19,len(members))
            color = gold if i==0 else COLORS[model] if i==1 else MUTED
            bx.scatter(i+jitter,values,s=28,c=color,alpha=.8,edgecolors=BG,linewidths=.4)
            median = float(np.median(values))
            bx.plot([i-.25,i+.25],[median,median],color=TEXT,lw=2)
            bx.text(i+.30,median,f"{median:.0f}%",size=9,va="center",color=TEXT)
        bx.set(xlim=(-.5,2.6),ylim=(24,105),yticks=[25,50,75,100],xticks=[0,1,2],
               xticklabels=["Resolved jets\nn = 7","Other optical AGN\nn = 14","HII comparison\nn = 29"])
        bx.tick_params(axis="x",labelsize=9)
        bx.tick_params(axis="y",labelsize=9)
        bx.grid(axis="y",color=GRID,alpha=.55)
        bx.set_axisbelow(True)
        if not col:
            bx.set_ylabel("Concentration [%]",size=10)
    fig.text(.055,.353,"Context: all measured classifications  ·  white bars show group medians",color=MUTED,size=11)
    fig.text(.055,.111, r"Concentration = max$_i$ |w$_i$ M$_{\mathrm{ref},i}$| / Σ$_i$ |w$_i$ M$_{\mathrm{ref},i}$|.  "
             "Signed coefficients; absolute strength retains zero-net-mass ψ lobes.",size=10,color=MUTED)
    fig.text(.055,.079,"1.5-GHz LeMMINGs cores rescaled to SPARC distances: radio activity is a proxy, not mechanical jet power. Limits are not detections.",size=9.5,color=MUTED)
    fig.text(.055,.050,"124 unmatched galaxies and 1 uncertain ALG nucleus excluded. HII classification is not proof of no AGN or jet; no secure non-jet null sample.",size=9.5,color=MUTED)
    fig.text(.055,.022,"Rotation-only fits, up to four terms, n ≤ 3 and |w| ≤ 20. These descriptive comparisons do not identify unique orbital populations or a jet correlation.",size=9.5,color=MUTED)
    save(fig,"sparc_jet_concentration")
