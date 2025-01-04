import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from collections import Counter


def aa_freq_counter(seqs, max_length=13):
    frequency_data = []

    for position in range(max_length):
        position_counter = Counter()

        for seq in seqs:
            if position < len(seq):
                position_counter[seq[position]] += 1

        frequency_data.append(position_counter)

    frequency_df = pd.DataFrame(frequency_data).fillna(0)

    return frequency_df


def logo_maker(frequency_df, label='', if_save=False, name_to_fig=''):
    logo = logomaker.Logo(frequency_df, fade_below=0.01, color_scheme='chemistry')  # hydrophobicity#charge#chemistry
    logo.ax.set_ylabel('Frequency')
    logo.ax.set_xlabel('Position')
    plt.title(f'CDR3 Sequence Logo for {label.upper()}')
    plt.tight_layout()
    plt.show()

    if if_save and name_to_fig:
        plt.savefig(name_to_fig)
        print(f"Logo saved as {name_to_fig}")


# Чтение данных из файла позитивной фракции
pos_file = r"C:\Users\Роман\Documents\Science\Moscow\Ph.D\Olya\YFV_B_Pos_cleanup.csv"

pos_seqs = pd.read_csv(pos_file, header=None)[0].tolist()

# Расчёт частот аминокислот для позитивной фракции
frequency_pos_B_filtered = aa_freq_counter(pos_seqs)

# Создание логотипа для позитивной фракции
logo_maker(frequency_pos_B_filtered, 'YFV POS', True, 'YFV_B_Pos_cleanup.png')
