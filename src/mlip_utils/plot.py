import matplotlib.pyplot as plt

import numpy as np

from ase_atom_tools import aseTools as at


def RMSE(db,dataset_name,labs,energy_label,forces_label,unit):
    maps = ['viridis','coolwarm']

    if unit == "meV":
        scale = 1000.0
    elif unit == "eV":
        scale = 1.0
    else:
        raise ValueError(f"Unknown unit: {unit}")

    # ==== 1. get RMSE energy====

    reference_energy = at.get_total_energy(db,energy_label[0],True).flatten() * scale
    model_energy = at.get_total_energy(db,energy_label[1],True).flatten() * scale

    rmse_energy = np.sqrt(np.mean((reference_energy-model_energy)**2))

    # ==== 2. get RMSE forces ====
    reference_forces = np.concatenate(at.get_forces(db,forces_label[0])).flatten() * scale
    model_forces = np.concatenate(at.get_forces(db,forces_label[1])).flatten() * scale

    rmse_forces = np.sqrt(np.mean((reference_forces - model_forces)**2))

    # ==== 3. Relative Forces ==========
    rel_rmse_forces = np.sqrt(np.sum((reference_forces - model_forces)**2)) / np.sqrt(np.sum((reference_forces - np.mean(reference_forces))**2)) * 100.0

    # ==== PLOT ======

    plt.figure(figsize=(9,6),dpi=100)
    plt.subplot(1,2,1)
    plt.title(f"Energy per atom ($\\mathrm{{{unit}/atom}}$) of {dataset_name}")

    # ---- ADD RMSE TEXT ----
    ax = plt.gca()
    ax.text(
        0.5, 0.1,                           # x,y in axes coordinates (0–1)
        f"RMSE  = {rmse_energy:.5f}",
        transform=ax.transAxes,             # use Axes coordinates
        ha='center', va='center', fontsize=12
    )
    # --------------------------------

    plt.xlabel(labs[0])
    plt.ylabel(labs[1])
    plt.scatter(reference_energy,model_energy,s=4,c=model_energy,cmap=maps[0],rasterized=True,edgecolors='none')

    # reference line
    vmin = min(reference_energy.min(),model_energy.min())
    vmax = max(reference_energy.max(),model_energy.max())

    plt.plot([vmin, vmax], [vmin, vmax], linestyle="-", alpha=0.6)

    plt.subplot(1,2,2)
    plt.title(f"Forces per atom ($\\mathrm{{{unit}/\\AA}}$) of {dataset_name}")

    ax = plt.gca()
    ax.text(
        0.5,0.1,
        f"RMSE = {rmse_forces:.5f}\nRRMSE = {rel_rmse_forces:.2f} %",
        transform=ax.transAxes,
        ha = 'center', va = 'center', fontsize=12
    )

    plt.xlabel(labs[0])
    plt.ylabel(labs[1])

    plt.scatter(reference_forces,model_forces,s=4,c=model_forces,cmap=maps[0],rasterized=True,edgecolors='none')

    # reference line
    vmin2 = min(reference_forces.min(),model_forces.min())
    vmax2 = max(reference_forces.max(),model_forces.max())

    plt.plot([vmin2, vmax2], [vmin2, vmax2], linestyle="-", alpha=0.6)

    plt.tight_layout()

    plt.show()