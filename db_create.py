from pathlib import Path
import sqlite3
import pandas as pd
import re

# -------------------------------------------------------
# SETTINGS
# -------------------------------------------------------

ROOT = Path("D:/Andreas/Databases")
DB = "baringo_spectra.db"

# ==========================================================
# CONNECT
# ==========================================================

conn = sqlite3.connect(DB)
cur = conn.cursor()


def clean_name(name):
    """
    Convert folder/file names into valid SQLite table names.
    """
    name = re.sub(r"\W+", "_", name)
    if name[0].isdigit():
        name = "_" + name
    return name


def create_table(table_name):

    cur.execute(f'DROP TABLE IF EXISTS "{table_name}"')

    cur.execute(f'''
        CREATE TABLE "{table_name}"(
            Id INTEGER PRIMARY KEY AUTOINCREMENT,
            wavelength DOUBLE,
            reflectance DOUBLE,
            direct_transmittance DOUBLE,
            diffuse_transmittance DOUBLE
        );
    ''')


# ==========================================================
# IMPORT
# ==========================================================

for site_dir in ROOT.iterdir():

    if not site_dir.is_dir():
        continue

    site = site_dir.name

    for species_dir in site_dir.iterdir():

        if not species_dir.is_dir():
            continue

        species = species_dir.name

        for txt in species_dir.glob("*.txt"):

            sample = txt.stem

            table = clean_name(f"{site}_{species}_{sample}")

            print(f"Creating table {table}")

            create_table(table)

            df = pd.read_csv(
                txt,
                sep=r"\s+",
                engine="python"
            )

            rows = [
                (
                    r.wavelength,
                    r.reflectance,
                    r.direct_transmittance,
                    r.diffuse_transmittance
                )
                for r in df.itertuples(index=False)
            ]

            cur.executemany(
                f'''
                INSERT INTO "{table}"
                (
                    wavelength,
                    reflectance,
                    direct_transmittance,
                    diffuse_transmittance
                )
                VALUES (?,?,?,?)
                ''',
                rows
            )

conn.commit()
conn.close()

print("Finished.")