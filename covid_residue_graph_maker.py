import numpy as np
import matplotlib.pyplot as plt

def sum_over(arr):
    sum = 0
    summed_arr = []
    for i in range(len(arr)):
        sum += arr[i]
        summed_arr.append(sum)

    return summed_arr

def LinkedInsertionSort(arr1, arr2):

    for i in range(1, len(arr1)):
        for j in range(i, 0, -1):
            if arr1[j] < arr1[j-1]:
                temp = arr1[j]
                arr1[j] = arr1[j-1]
                arr1[j-1] = temp

                temp = arr2[j]
                arr2[j] = arr2[j-1]
                arr2[j-1] = temp

    return (arr1, arr2)

posAAs = ["ARG", "LYS"]
negAAs = ["ASP", "GLU"]
chargedAAs = ["ARG", "LYS", "ASP", "GLU"]

data = open("ace2_spike_cloc_energies_fsapt_all.csv", "r")

datalines = data.readlines()

energy = []
distance = []
residueACE2 = []
residueSPIKE = []

for i in range(1, len(datalines)):
    line = datalines[i].strip("\n").split(',')
    energy.append(float(line[7]))
    distance.append(float(line[4]))
    residueACE2.append(line[1])
    residueSPIKE.append(line[3])

LinkedInsertionSort(distance, energy)

distance_pos_neg = []
energy_pos_neg = []

for i in range(len(distance)):
    if (residueACE2[i] in posAAs and residueSPIKE[i] in negAAs) or (residueACE2[i] in negAAs and residueSPIKE[i] in posAAs):
        distance_pos_neg.append(distance[i])
        energy_pos_neg.append(energy[i])

distance_lc = []
energy_lc = []

for i in range(len(distance)):
    if (residueACE2[i] in posAAs and residueSPIKE[i] in posAAs) or (residueACE2[i] in negAAs and residueSPIKE[i] in negAAs):
        distance_lc.append(distance[i])
        energy_lc.append(energy[i])

distance_neutral = []
energy_neutral = []

for i in range(len(distance)):
    if (residueACE2[i] not in posAAs and residueACE2[i] not in negAAs and residueSPIKE[i] not in posAAs and residueSPIKE[i] not in negAAs):
        distance_neutral.append(distance[i])
        energy_neutral.append(energy[i])

distance_cn = []
energy_cn = []

for i in range(len(distance)):
    if (residueACE2[i] in chargedAAs and residueSPIKE[i] not in chargedAAs) or (residueACE2[i] not in chargedAAs and residueSPIKE[i] in chargedAAs):
        distance_cn.append(distance[i])
        energy_cn.append(energy[i])

energy_summed, energy_pos_neg_summed, energy_lc_summed, energy_neutral_summed, energy_cn_summed = sum_over(energy), sum_over(energy_pos_neg), sum_over(energy_lc), sum_over(energy_neutral), sum_over(energy_cn)

plt.plot(distance, energy_summed, label="All")
plt.plot(distance_pos_neg, energy_pos_neg_summed, label="Pos-Neg")
plt.plot(distance_lc, energy_lc_summed, label="Like Charges")
plt.plot(distance_neutral, energy_neutral_summed, label="Neutral-Neutral")
plt.plot(distance_cn, energy_cn_summed, label="Neutral-Charged")
plt.legend()

plt.title("Cumulative Interaction Energies")
plt.xlabel("Close Contact Distance (Angstroms)")
plt.ylabel("Cumulative SAPT0 Energy (kcal/mol)")

plt.savefig("ace2_spike_cloc_energies_fsapt_cumulative_ALL5.png")
plt.close("all")

def graph_all():

    distance = np.array(distance, dtype=float)
    energy = np.array(energy, dtype=float)

    plt.plot(distance, energy)

    plt.xlabel("Close Contact Distance (Angstroms)")
    plt.ylabel("SAPT0 Energy (Adjusted for Capping, kcal/mol)")

    plt.savefig("ace2_spike_cloc_energies_fsapt_all.png")

    plt.close("all")

def graph_pos_neg():

    distance_pos_neg = []
    energy_pos_neg = []

    for i in range(len(distance)):
        if (residueACE2[i] in posAAs and residueSPIKE[i] in negAAs) or (residueACE2[i] in negAAs and residueSPIKE[i] in posAAs):
            distance_pos_neg.append(distance[i])
            energy_pos_neg.append(energy[i])

    distance_pos_neg = np.array(distance_pos_neg, dtype=float)
    energy_pos_neg = np.array(energy_pos_neg, dtype=float)

    plt.plot(distance_pos_neg, energy_pos_neg)

    plt.xlabel("Close Contact Distance (Angstroms)")
    plt.ylabel("SAPT0 Energy (Adjusted for Capping, kcal/mol)")

    plt.savefig("ace2_spike_cloc_energies_fsapt_salt_bridge.png")

    plt.close("all")

def graph_like_charges():

    distance_lc = []
    energy_lc = []

    for i in range(len(distance)):
        if (residueACE2[i] in posAAs and residueSPIKE[i] in posAAs) or (residueACE2[i] in negAAs and residueSPIKE[i] in negAAs):
            distance_lc.append(distance[i])
            energy_lc.append(energy[i])

    print(len(distance_lc))

    plt.scatter(distance_lc, energy_lc)
    plt.autoscale()

    plt.title("Like Charges")
    plt.xlabel("Close Contact Distance (Angstroms)")
    plt.ylabel("SAPT0 Energy (Adjusted for Capping, kcal/mol)")

    plt.savefig("ace2_spike_cloc_energies_fsapt_like_charges.png")

    plt.close("all")

def graph_neutral_residues():

    distance_neutral = []
    energy_neutral = []

    for i in range(len(distance)):
        if (residueACE2[i] not in posAAs and residueACE2[i] not in negAAs and residueSPIKE[i] not in posAAs and residueSPIKE[i] not in negAAs):
            distance_neutral.append(distance[i])
            energy_neutral.append(energy[i])

    print(len(distance_neutral))

    plt.scatter(distance_neutral, energy_neutral)

    plt.title("Neutral Residue - Neutral Residue")

    plt.xlabel("Close Contact Distance (Angstroms)")
    plt.ylabel("SAPT0 Energy (Adjusted for Capping, kcal/mol)")

    plt.savefig("ace2_spike_cloc_energies_fsapt_neutral_neutral.png")

    plt.close("all")

def graph_charge_neutral():

    distance_cn = []
    energy_cn = []

    for i in range(len(distance)):
        if (residueACE2[i] in chargedAAs and residueSPIKE[i] not in chargedAAs) or (residueACE2[i] not in chargedAAs and residueSPIKE[i] in chargedAAs):
            distance_cn.append(distance[i])
            energy_cn.append(energy[i])

    print(len(distance_cn))

    plt.scatter(distance_cn, energy_cn)

    plt.title("Charged Residue - Neutral Residue")

    plt.xlabel("Close Contact Distance (Angstroms)")
    plt.ylabel("SAPT0 Energy (Adjusted for Capping, kcal/mol)")

    plt.savefig("ace2_spike_cloc_energies_fsapt_charged_neutral.png")

    plt.close("all")
